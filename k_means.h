#ifndef __K_MEANS_H
#define __K_MEANS_H

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/*<>Measures the Euclidean distance between two vectors v1, v2 with dim components<>*/
float distance(int dim,  	//no. dimensions 
               float *v1,   //[numdims] 
               float *v2);  //[numdims]

/*<>Returns the index of the cluster that is closest to the given vector<>*/
int getNearestClusterIndex(int   numClusters, 	//no. clusters 
                           int   numCoords,   	//no. coordinates 
                           float  *vector,      //[numCoords] 
                           float **clusters);    //[numClusters][numCoords]

/*<>Computes K clusters and stores them in **clusters based on vectorial group 
**vectors and the limit - maximun number of algorithm iterations allowed<>*/
int k_means(float    **vectors,     			 //in: [numVectors][numDims]
            int        numDims,   				 //num of coordinates in vector 
            int        numVectors,   
            int        numClusters, 
            float      limit,   			 	 //max num of iterations	
            int       *vectorToClusterRelevance, //out: [numVectors] 
            float    **clusters,    			 //out: [numClusters][numDims] 
            MPI_Comm   comm);



#endif //__K_MEANS_H