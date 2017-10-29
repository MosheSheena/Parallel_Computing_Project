#pragma warning( disable : 4996 )
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

float** readVectorsFromFile(const char *fileName,				//file containing the data set     			
											int        *numDims,   				//number of dims in vectorial space
											int        *numVectors,   			//number of given vectors in file
											int        *maxNumOfClusters, 	//MAX number of cluster allowed according to file
											int        *iterationLimit,  				//limit of k-means iteration allowed
											float	   *qualityOfClusters) 		//quality of clusters to find according to file   		
{
	int i, j, counter = 0;
	float **vectors;
	FILE *f;

	f = fopen(fileName, "r");
	assert(f != NULL);

	fscanf(f, "%d\t%d\t%d\t%d\t%f\n", numVectors, numDims, maxNumOfClusters, iterationLimit, qualityOfClusters);

	//assiging the vectors matrix
	vectors = (float**)malloc(*numVectors * sizeof(float*));
	assert(vectors != NULL);

	//one call to calloc instead of calling it for each pointer in vectors[][]
	vectors[0] = (float*)calloc((*numVectors) * (*numDims), sizeof(float));
	assert(vectors[0] != NULL);

	for (i = 1; i < *numVectors; ++i)
	{
		vectors[i] = vectors[i - 1] + (*numDims);
	}

	for (i = 0; i < *numVectors; ++i)
	{
		for (j = 0; j < *numDims; ++j)
		{
			fscanf(f, "%f\t", &vectors[i][j]);
		}
		fscanf(f, "\n");
	}

	fclose(f);
	return vectors;
}

void printVectors(float** vectors, int numVectors, int numDims)
{
	int i, j;
	for (i = 0; i < numVectors; ++i)
	{
		for (j = 0; j < numDims; ++j)
			printf("%f\t", vectors[i][j]);

		printf("\n");
	}
}


