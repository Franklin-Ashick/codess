#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <omp.h>
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
double **calculate_distance(double **arr,int length)
{
        double **dMatrix = (double **)malloc(length * sizeof(double *));

        for(int i = 0; i < length; i++){
                dMatrix[i] = (double *) malloc(length * sizeof(double));
                if (dMatrix[i] == NULL){
                        perror("Memory Allocation Failed");
                }
        }


    // Parallelize both loops with collapse
   #pragma omp parallel for collapse(2)
   for (int i = 0; i < length; i++) {
                for(int j =0; j<length; j++){

            dMatrix[i][j] = sqrt(pow(arr[i][0]-arr[j][0], 2) +
                                 pow(arr[i][1]-arr[j][1],2));
//            printf(" m[%d][%d] : (%lf)\n", i,j,dMatrix[i][j]);
               }
        }

   return dMatrix;
}

// Farthest Insertion TSP Algorithm with OpenMP
void farthestInsertionTSP(double **distMatrix, int *tour, int startVertex, int N) {
    bool *visited = (bool *)calloc(N, sizeof(bool));
    int tourSize = 1;
    tour[0] = startVertex;
    visited[startVertex] = true;

    // Main loop to insert all vertices into the tour
    while (tourSize < N) {
        double maxDistance = -1.0;
        int maxVertex = -1;
        // Declare structure to hold local results
        struct {
            double maxDistance;
            int maxVertex;

        } local;

        #pragma omp parallel private(local)
        {
            local.maxDistance = -1.0;
            local.maxVertex = -1;
            // Parallel for loop to find the farthest unvisited vertex from the tour
            #pragma omp for nowait
            for (int i = 0; i < N; i++) {
                if (!visited[i]) {
                    for (int j = 0; j < tourSize; j++) {
                        if (distMatrix[tour[j]][i] > local.maxDistance) {
                            local.maxDistance = distMatrix[tour[j]][i];
                            local.maxVertex = i;
                        }
                    }
                }
            }

            // Critical section to update the global maximum
            #pragma omp critical
            {
                if (local.maxDistance > maxDistance) {
                    maxDistance = local.maxDistance;
                    maxVertex = local.maxVertex;
                }
            }
        } // End of parallel region

        // Sequentially find the best position to insert the farthest vertex
        double minCost = DBL_MAX;
        int minPos = -1;
        for (int i = 0; i < tourSize; i++) {
            int next = (i + 1) % tourSize;
            double cost = calculateInsertionCost(maxVertex, tour[i], tour[next], distMatrix);
            if (cost < minCost) {
                minCost = cost;
                minPos = i + 1;
            }
        }

        // Sequentially insert the farthest vertex at the best position
        insertVertexAtPosition(tour, minPos, maxVertex, tourSize);
        tourSize++;
        visited[maxVertex] = true;
    }

    free(visited);
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

    gettimeofday(&wallStart, NULL); // save start time in to variable 'wallStart'

    // Solve the TSP using Farthest Insertion with OpenMP
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

    printf("Time taken for solving TSP using farthest insertion, parallely  =  %f seconds  \n", WALLtimeTaken);


    writeTourToFile(tour, N+1, outFileName);

    // Free allocated memory
    for (int i = 0; i < N; i++) {
        free(distMatrix[i]);
    }

    free(distMatrix);
    free(tour);
    free(coords_twod_array);

 return 0;
}