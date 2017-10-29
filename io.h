#ifndef __IO_H
#define __IO_H

#include <assert.h>

float** readVectorsFromFile(const char *fileName,				//file containing the data set     			
							int        *numDims,   				//number of dims in vectorial space
							int        *numVectors,   			//number of given vectors in file
							int        *maxNumOfClusters, 	//MAX number of cluster allowed according to file
							int        *iterationLimit,  				//limit of k-means iteration allowed
							float	   *qualityOfClusters) 		//quality of clusters to find according to file

void printVectors(float** vectors, int numVectors, int numDims);

#endif //__IO_H