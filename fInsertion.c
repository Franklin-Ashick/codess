
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>  // for wallclock timing functions

int readNumOfCoords(char *fileName);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);


// Function to calculate the insertion cost
double calculateInsertionCost(int vertex, int before, int after, double **distMatrix) {
    return distMatrix[before][vertex] + distMatrix[vertex][after] - distMatrix[before][after];
}


// Function to insert the vertex into the tour at the specified position
void insertVertexAtPosition(int *tour, int position, int vertex, int tourSize) {
    for (int i = tourSize; i > position; i--) {
        tour[i] = tour[i - 1];
    }
    tour[position] = vertex;
}

// Function to calculate distance
double **calculate_distance(double **arr, int length) {
    double **dMatrix = (double **)malloc(length * sizeof(double *));
    for (int i = 0; i < length; i++) {
        dMatrix[i] = (double *)malloc(length * sizeof(double));
        for (int j = 0; j < length; j++) {
            dMatrix[i][j] = sqrt(pow(arr[i][0] - arr[j][0], 2) +
                                 pow(arr[i][1] - arr[j][1], 2));
        }
    }
    return dMatrix;
}

// Farthest Insertion TSP Algorithm
void farthestInsertionTSP(double **distMatrix, int *tour, int startVertex, int N) {
    bool *visited = (bool *)calloc(N, sizeof(bool));

    int tourSize = 1;
    tour[0] = startVertex;
    visited[startVertex] = true;

    while (tourSize < N) {
        double maxDistance = -1;
        int maxVertex = -1;
        // Find the farthest vertex from the tour
        for (int i = 0; i < N; i++) {
            if (!visited[i]) {
                for (int j = 0; j < tourSize; j++) {
                    if (distMatrix[tour[j]][i] > maxDistance) {
                        maxDistance = distMatrix[tour[j]][i];
                        maxVertex = i;
                    }
                }
            }
        }

        // Find the best position to insert maxVertex
        double minCost = DBL_MAX;
        int minPos = -1;
        for (int j = 0; j < tourSize; j++) {
            int k = (j + 1) % tourSize;
            double cost = calculateInsertionCost(maxVertex, tour[j], tour[k], distMatrix);

           if (cost < minCost) {
                minCost = cost;
                minPos = j + 1;
            }
        }

        // Insert maxVertex at minPos position
        insertVertexAtPosition(tour, minPos, maxVertex, tourSize);
        tourSize++;
        visited[maxVertex] = true;
    }
}

int main(char argc, char *argv[]) {
    char *fileName = argv[1];
    char *outFileName = argv[2];

    int N = readNumOfCoords(fileName);
    double **coords_twod_array = readCoords(fileName, N);
    double **distMatrix = calculate_distance(coords_twod_array, N);

    // The tour array
    int *tour = (int *)malloc(N * sizeof(int));
    /* for timing */
    struct timeval wallStart, wallEnd;

    gettimeofday(&wallStart, NULL); // save start time into variable 'wallStart'



    // Solve the TSP using Farthest Insertion
    farthestInsertionTSP(distMatrix, tour, 0, N); // Start from vertex 0

    gettimeofday(&wallEnd, NULL); // end time
    double wallSecs = (wallEnd.tv_sec - wallStart.tv_sec);           // just integral number of seconds
    double WALLtimeTaken = 1.0E-06 * ((wallSecs*1000000) + (wallEnd.tv_usec - wallStart.tv_usec)); // and now with any microseconds

    tour[N] = 0; // Closing the tour by returning to the start

    // Print the tour

//    printf("Tour: ");
//    for (int i = 0; i < N + 1; i++) {
//        printf("%d ", tour[i]);
//    }
//    printf("\n");
    printf("Time taken for solving TSP using farthest insertion, serially  =  %f seconds  \n", WALLtimeTaken);

    writeTourToFile(tour, N + 1, outFileName);

    free(distMatrix);
    free(tour);
    free(coords_twod_array);

    return 0;
}