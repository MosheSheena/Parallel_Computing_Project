#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//consts
static const int VECTOR_SIZE = 2;
static const int K = 2;
static const int NUM_VECTORS = 8;

struct Cluster;

//Structs
typedef struct Vector
{
    int numComponents;
    double *components = NULL;
    Cluster *cluster = NULL;
};

typedef struct Cluster
{
    int numVectors = 0;
    Vector *centroid = NULL;
    Vector **allVectors = NULL;
    double diameter;
};

//Prototypes
Cluster **kMeans(int k, Vector **vectorArray, int numVectors, int limit);
Vector* computeCentroid(const Vector *vectorArray, int numVectors);
double measureDistance(Vector *v1, Vector *v2);
int assignVectorToNearestCluster(Cluster **allClusters, int numClusters, Vector *v);
double computeClusterDiameter(Cluster *cluster);
int compareVectors(Vector *v1, Vector *v2);
void initClustersInKMeans(int k, Cluster **clusters, Vector **vectorArray);
void addVectorToCluster(Cluster *cluster, Vector *vector);
void setCentroid(Cluster *cluster, Vector *centroid);
void freeClusterMemory(Cluster *c);
int removeVectorFromCluster(Cluster *c, Vector *v);

int main()
{
    int limit = 100;
    int x = 1;
    Vector **vArr = (Vector**)calloc(NUM_VECTORS, sizeof(Vector*));
    for (int i = 0; i < NUM_VECTORS; ++i)
    {
        if(i == NUM_VECTORS/2) {x+=20;}

        vArr[i] = (Vector*)calloc(1, sizeof(Vector));
        vArr[i]->numComponents = VECTOR_SIZE;
        vArr[i]->components = (double*)malloc(VECTOR_SIZE*sizeof(double));
        vArr[i]->components[0] = x;
        vArr[i]->components[1] = x+1;
        x++;
    }

    for (int j = 0; j < NUM_VECTORS; ++j)
    {
        printf("(%lf, %lf)\n", vArr[j]->components[0], vArr[j]->components[1]);
    }

    Cluster **clusters = kMeans(K, vArr, NUM_VECTORS, limit);


    return 0;
}

//Measures the distance between two given vectors using euclidean distance
double measureDistance(Vector *v1, Vector *v2)
{
    double distance, sum = 0;

    //OMP
    for(int i = 0; i < VECTOR_SIZE; ++i)
    {
        sum += pow(v2->components[i] - v1->components[i], 2);
    }

    distance = sqrt(sum);

    return distance;
}

//computes the centroid of a given group of vectors
Vector* computeCentroid(Vector *vectorArray, int numVectors)
{
    double compSum;
    Vector *newCentroid = (Vector*)calloc(1, sizeof(Vector));

    //OMP
    for (int i = 0; i < VECTOR_SIZE; ++i)
    {
        compSum = 0;
        for(int j = 0; j < numVectors; ++j)
        {
            compSum += vectorArray[j].components[i];
        }

        newCentroid->components[i] = compSum / numVectors;
    }

    return newCentroid;
}

double computeClusterDiameter(Cluster *cluster)
{
    double diameter, distance;

    diameter = measureDistance(cluster->allVectors[0], cluster->allVectors[1]);

    for (int i = 0; i < cluster->numVectors; ++i)
    {
        for (int j = i + 1; j < cluster->numVectors; ++j)
        {
            distance = measureDistance(cluster->allVectors[i], cluster->allVectors[j]);
            if(distance > diameter)
            {
                diameter = distance;
            }
        }
    }
    return diameter;
}

int compareVectors(Vector *v1, Vector *v2) //Returns 1 if vectors are different otherwise returns 0
{
    for (int i = 0; i < VECTOR_SIZE; ++i)
    {
       if(v1->components[i] != v2->components[i])
           return 1;
    }

    return 0;
}

Cluster **kMeans(int k, Vector **vectorArray, int numVectors, int limit)
{
    //Create all clusters
    Cluster **clusters = (Cluster**)calloc(k, sizeof(Cluster*));
    for (int i = 0; i < k; ++i)
    {
        clusters[i] = (Cluster*)calloc(1, sizeof(Cluster));
    }

    initClustersInKMeans(k, clusters, vectorArray);

    int reposition = 0; //Keeps num of vectors that moved between clusters
    int numIterations = 0;

    do
    {
        for (int i = 0; i < numVectors; ++i)
        {
            reposition += assignVectorToNearestCluster(clusters, k, vectorArray[i]);
        }
        if(reposition == 0)
            break;//Each cluster keeps ALL of his vectors

        //Compute centroid in every cluster
        for (int j = 0; j < k; ++j)
        {
            if(clusters[j]->numVectors == 0)
                continue;//Skips 1 iteration

            clusters[j]->centroid = computeCentroid(*clusters[j]->allVectors, clusters[j]->numVectors);
        }

        numIterations++;

    }while(numIterations < limit);

    return clusters;
}

int assignVectorToNearestCluster(Cluster **allClusters, int numClusters, Vector *v)//returns 1 for change or 0 for no-change
{
    double minDistance, distance;
    int closestClusterIndex;

    minDistance  = measureDistance(allClusters[0]->centroid, v);
    closestClusterIndex = 0;

    for (int i = 1; i < numClusters; ++i)
    {
        distance = measureDistance(allClusters[i]->centroid, v);
        if(distance < minDistance)
        {
            minDistance = distance;
            closestClusterIndex = i;
        }
    }

    if(!v->cluster)
    {
        addVectorToCluster(allClusters[closestClusterIndex], v);
        return 1;
    }
    else if(compareVectors(v->cluster->centroid, allClusters[closestClusterIndex]->centroid) != 0)
    {
        removeVectorFromCluster(v->cluster, v);
        addVectorToCluster(allClusters[closestClusterIndex], v);
        return 1;
    }
    else
        return 0;
}

//determines the K first vectors as cluster centroids
void initClustersInKMeans(int k, Cluster **clusters, Vector **vectorArray)
{
    for (int i = 0; i < k; ++i)
    {
        clusters[i]->centroid = vectorArray[i];
    }
}

void addVectorToCluster(Cluster *cluster, Vector *vector)
{
    if(!cluster->allVectors)
    {
        cluster->allVectors = (Vector**)calloc(1, sizeof(Vector*));
    }
    else
    {
        cluster->allVectors = (Vector**)realloc(cluster->allVectors, sizeof(Vector*)*(cluster->numVectors + 1));
    }

    cluster->allVectors[cluster->numVectors++] = vector;

    vector->cluster = cluster;
}

void setCentroid(Cluster *cluster, Vector *centroid)
{
    cluster->centroid = centroid;
}

void freeClusterMemory(Cluster *c)
{
    free(c->allVectors);
}

int removeVectorFromCluster(Cluster *c, Vector *v)//Returns 1 if vector was successfully removed, otherwise 0
{
    for (int i = 0; i < c->numVectors; ++i)
    {
        if(compareVectors(c->allVectors[i], v) == 0)
        {
            for (int j = i; j < c->numVectors - 1; ++j)
            {
                c->allVectors[j] = c->allVectors[j+1];
            }

            c->numVectors--;
            c->allVectors = (Vector**)realloc(c->allVectors, sizeof(Vector*)*c->numVectors);

            return 1;
        }
    }
    return 0;
}





