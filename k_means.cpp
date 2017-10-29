#include "k_means.h"

float euclidDistance(int dim,  	//no. dimensions 
               float *v1,   //[numdims] 
               float *v2)   //[numdims] 
{
    int i;
    float dist = 0.0;

    for (i=0; i < dim; ++i)
    {
        dist += (v1[i]-v2[i]) * (v1[i]-v2[i]);
    }

    return dist;
}

int getNearestClusterIndex(int   numClusters, 	//no. clusters 
                           int   numDims,   	//no. coordinates 
                           float  *vector,      //[numDims] 
                           float **clusters)    //[numClusters][numDims] 
{
    int index, i;
    float distance, minDistance;

    //init first cluster as closest to vector
    index    = 0;
	
    minDistance = euclidDistance(numDims, vector, clusters[0]);

    for (i = 1; i < numClusters; ++i) 
    {
        distance = euclidDistance(numDims, vector, clusters[i]);

        //square root is not computed - explained in doc

        if (distance < minDistance) 
        { 
            minDistance = distance;
            index = i;
        }
    }
    return index;
}

int k_means(float    **vectors,     			 //in: [numVectors][numDims]
            int        numDims,   				 //num of coordinates in vector 
            int        numVectors,   
            int        numClusters, 
            float      limit,   			 	 //max num of iterations	
            int       *vectorToClusterRelevance, //out: [numVectors] 
            float    **clusters,    			 //out: [numClusters][numDims] 
            MPI_Comm   comm)        
{
    int      i, j, rank, index, loop = 0, total_numVectors;
    int     *newClusterSize; //[numClusters]: no. vectors assigned in each new cluster                             
    int     *clusterSize;    //[numClusters]: temp buffer for reduction 
    float    delta;          //num of vectors that change their cluster
    float    delta_tmp;
    float  **newClusters;    //[numClusters][numDims] 

    //initialize vectorToClusterRelevance[]
    for (i=0; i < numVectors; ++i)
    {
    	vectorToClusterRelevance[i] = -1;
    }

    //initializing newClusterSize and newClusters[0] to 0
    newClusterSize = (int*) calloc(numClusters, sizeof(int));
    assert(newClusterSize != NULL);

    clusterSize    = (int*) calloc(numClusters, sizeof(int));
    assert(clusterSize != NULL);

    newClusters    = (float**) malloc(numClusters * sizeof(float*));
    assert(newClusters != NULL);

    newClusters[0] = (float*)  calloc(numClusters * numDims, sizeof(float));
    assert(newClusters[0] != NULL);

    for (i=1; i<numClusters; ++i)
    {
        newClusters[i] = newClusters[i-1] + numDims;
    }

    //not sure it's necessary
    //MPI_Allreduce(&numVectors, &total_numVectors, 1, MPI_INT, MPI_SUM, comm);

    do 
    {
        double timeStart = MPI_Wtime();

        delta = 0.0;

        for (i = 0; i < numVectors; ++i) 
        {
            //find the cluster nearest to the vector
            index = getNearestClusterIndex(numClusters, numDims, vectors[i], clusters);

            /*if the index found doesn't match the one in 
            vectorToClusterRelevance[i] - vector i needs to change cluster 
            and delta increase delta by 1 */
            if (vectorToClusterRelevance[i] != index) delta += 1.0;

            //update new cluster for vector i
            vectorToClusterRelevance[i] = index;

            //update new cluster center by summing the added vector i
            newClusterSize[index]++;
            for (j = 0; j < numDims; ++j)
            {
                newClusters[index][j] += vectors[i][j];
            }
        }

        /* sum all data vectors in newClusters */
        MPI_Allreduce(newClusters[0], clusters[0], numClusters*numDims,
                      MPI_FLOAT, MPI_SUM, comm);
        MPI_Allreduce(newClusterSize, clusterSize, numClusters, MPI_INT,
                      MPI_SUM, comm);

        /* average the sum and replace old cluster centers with newClusters */
        for (i=0; i<numClusters; i++) {
            for (j=0; j<numDims; j++) 
            {
                if (clusterSize[i] > 1)
                    clusters[i][j] /= clusterSize[i];
                newClusters[i][j] = 0.0;   /* set back to 0 */
            }
            newClusterSize[i] = 0;   /* set back to 0 */
        }
            
        MPI_Allreduce(&delta, &delta_tmp, 1, MPI_FLOAT, MPI_SUM, comm);
        delta = delta_tmp / total_numVectors;

        double maxTime;
        timeStart = MPI_Wtime() - timeStart;
        MPI_Reduce(&timeStart, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
        if (rank == 0) 
        {
        	printf("%2d: loop =%d time =%f sec\n",rank,loop,timeStart);
        }
        
    }while (delta > 0.0 || loop++ < limit);

    if (rank == 0) 
    {
    	printf("%2d: delta =%f limit =%f loop =%d\n",rank,delta,limit,loop);
    }

    //free all memory allocated
    free(newClusters[0]);
    free(newClusters);
    free(newClusterSize);
    free(clusterSize);

    return 0;
}







