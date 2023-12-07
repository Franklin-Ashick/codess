#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>  


int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);


double calculateInsertionCost(int vertex, int before, int after, double **distMatrix) {
    return distMatrix[before][vertex] + distMatrix[vertex][after] - distMatrix[before][after];
}

void insertVertexAtPosition(int *tour, int position, int vertex, int tourSize) {
    for (int i = tourSize; i > position; i--) {
        tour[i] = tour[i - 1];
    }
    tour[position] = vertex;
}

void cheapestInsertionTSP(double **distMatrix, int *tour, int startVertex, int N) {
    bool *visited = (bool *)calloc(N, sizeof(bool));

    int tourSize = 1;
    tour[0] = startVertex;

    visited[startVertex] = true;

    while (tourSize < N) {
        double minCost = DBL_MAX;
        int minPos = -1;
        int minVertex = -1;

        for (int i = 0; i < N; i++) {
            if (!visited[i]) {
                for (int j = 0; j < tourSize; j++) {
                    int k = (j + 1) % tourSize;
                 double cost = calculateInsertionCost(i, tour[j], tour[k], distMatrix);
                    if (cost < minCost) {
                        minCost = cost;
                        minVertex = i;
                        minPos = j + 1;
                    }
                }
            }
        }

        
        insertVertexAtPosition(tour, minPos, minVertex, tourSize);
        tourSize++;
        visited[minVertex] = true;
    }
}
double **calculate_distance(double **arr,int length)
{
        double **dMatrix = (double **)malloc(length* 100 * sizeof(double *));

        for(int i = 0; i < length; i++){
                dMatrix[i] = (double *) malloc(length*100 * sizeof(double));
                if (dMatrix[i] == NULL){
                       perror("Memory Allocation Failed");
                }
        }

   for (int i = 0; i < length; i++) {
                for(int j =0; j<length; j++){

            dMatrix[i][j] = sqrt(pow(arr[i][0]-arr[j][0], 2) +
                                 pow(arr[i][1]-arr[j][1],2));
                }
        }

   return dMatrix;
}

int main(int argc, char *argv[]){
    char *fileName = argv[1];
    char *outFileName = argv[2];

    int N = readNumOfCoords(fileName);

    double  **coords_twod_array = readCoords(fileName, N);
    double **distMatrix= calculate_distance(coords_twod_array, N);

    int *tour = (int *)malloc(N * sizeof(int));

    struct timeval wallStart, wallEnd;

    gettimeofday(&wallStart, NULL);


    cheapestInsertionTSP(distMatrix, tour, 0,N); 

    gettimeofday(&wallEnd, NULL); 
    double wallSecs = (wallEnd.tv_sec - wallStart.tv_sec);           
    double WALLtimeTaken = 1.0E-06 * ((wallSecs*1000000) + (wallEnd.tv_usec - wallStart.tv_usec)); 


    tour[N]=0;

    printf("Time taken for solving TSP using cheapest insertion, serially  =  %f seconds  \n", WALLtimeTaken);
    writeTourToFile(tour, N+1,outFileName);

    free(distMatrix);
    free(tour);
    free(coords_twod_array);

return 0;
}

