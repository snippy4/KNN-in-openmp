#include <mpi.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include"file-reader.c"
#include"knnomp.c"

void printArray(void *array, char array_name[500], int m, int n, char type);



int main(int argc, char *argv[]) {
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);


    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int max_procs = atoi(argv[4]);
    MPI_Comm new_comm = MPI_COMM_WORLD;
    int new_size = size;

    // processor amount control
    if (size > max_procs) {
        if (rank == 0) {
            printf("Reducing number of processes to %d.\n", max_procs);
        }
        if (rank < max_procs) {
            MPI_Comm_split(MPI_COMM_WORLD, 0, rank, &new_comm); 
        } else {
            MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, &new_comm);
        }

        if (rank >= max_procs) {
            MPI_Finalize();
            return 0;
        }

        MPI_Comm_size(new_comm, &new_size);
        MPI_Comm_rank(new_comm, &rank);
    }

    if (rank == 0) {
        printf("\n\n===============STARTING KNN WITH MPI===============\n\n");
    }

    // Variables
    char *inFile, *outFile;
    int k, numFolds;
    int totalNumPoints, numFeatures, numClasses;
    double *originalData;
    int *pointsInFold;
    
    if (rank == 0) {
        // Root process reads input and initializes variables
        inFile = argv[1];
        outFile = argv[2];
        k = atoi(argv[3]);
        numFolds = atoi(argv[4]);

        totalNumPoints = readNumOfPoints(inFile);
        numFeatures = readNumOfFeatures(inFile);
        numClasses = readNumOfClasses(inFile);
        originalData = readDataPoints(inFile, totalNumPoints, numFeatures);

        pointsInFold = (int *)malloc(numFolds * sizeof(int));
        int pointsPerFold = totalNumPoints / numFolds;
        int pointsPerFoldRemainder = totalNumPoints % numFolds;

        for (int fold = 0; fold < numFolds; fold++) {
            pointsInFold[fold] = pointsPerFold + ((fold == numFolds - 1) ? pointsPerFoldRemainder : 0);
        }
    }

    // Broadcast shared data to all processes
    MPI_Bcast(&totalNumPoints, 1, MPI_INT, 0, new_comm);
    MPI_Bcast(&numFeatures, 1, MPI_INT, 0, new_comm);
    MPI_Bcast(&k, 1, MPI_INT, 0, new_comm);
    MPI_Bcast(&numFolds, 1, MPI_INT, 0, new_comm);

    if (rank != 0) {
        pointsInFold = (int *)malloc(numFolds * sizeof(int));
        originalData = (double *)malloc(totalNumPoints * numFeatures * sizeof(double));
    }

    MPI_Bcast(pointsInFold, numFolds, MPI_INT, 0, new_comm);
    MPI_Bcast(originalData, totalNumPoints * numFeatures, MPI_DOUBLE, 0, new_comm);

    // Distribute folds among processes
    int distributed_folds = numFolds / new_size;
    int extra_folds = numFolds % new_size;
    int first_fold = rank * distributed_folds;
    int last_fold = first_fold + distributed_folds;
    double *results = (double *)malloc((last_fold - first_fold) * sizeof(double));
    size_t maxTrain_bytes = (totalNumPoints - pointsInFold[0]) * numFeatures * sizeof(double);
    size_t maxTest_bytes = pointsInFold[numFolds - 1] * numFeatures * sizeof(double);
    double *currTrain = (double *)malloc(maxTrain_bytes);
    double *currTest = (double *)malloc(maxTest_bytes);

   

    size_t testOffset = 0;
    for (int currFold = first_fold; currFold < last_fold; currFold++) {
        printf("fold: %d\n", currFold);
        testOffset = pointsInFold[0] * numFeatures * currFold;
        size_t currTestSize = pointsInFold[currFold] * numFeatures;
        size_t currTestSize_bytes = currTestSize * sizeof(double);
        size_t testOffset_bytes = testOffset * sizeof(double);

        memcpy(currTest, originalData + testOffset, currTestSize_bytes);

        size_t currTrainPoints = totalNumPoints - pointsInFold[currFold];
        size_t currTrainSize = currTrainPoints * numFeatures;
        size_t currTrainSize_bytes = currTrainSize * sizeof(double);

        memcpy(currTrain, originalData, testOffset_bytes);
        memcpy(currTrain + testOffset, originalData + testOffset + currTestSize, currTrainSize_bytes - testOffset_bytes);
        //printArray(currTest, "test", pointsInFold[currFold], numFeatures, 'd');
        //printArray(currTrain, "train", currTrainPoints, numFeatures, 'd');
        results[currFold - first_fold] = knn(currTrain, currTest, k, currTrainPoints, pointsInFold[currFold], numFeatures);
    }

    double *finalResults = NULL;
    if (rank == 0) {
        finalResults = (double *)malloc((numFolds + 1) * sizeof(double));
    }

    MPI_Gather(results, last_fold - first_fold, MPI_DOUBLE, finalResults, last_fold - first_fold, MPI_DOUBLE, 0, new_comm);

    if (rank == 0) {
        testOffset = distributed_folds * size * numFeatures * pointsInFold[0];
        for (int currFold = distributed_folds * size; currFold < numFolds; currFold++) {
            printf("fold: %d\n", currFold);
            size_t currTestSize = pointsInFold[currFold] * numFeatures;
            size_t currTestSize_bytes = currTestSize * sizeof(double);
            size_t testOffset_bytes = testOffset * sizeof(double);

            memcpy(currTest, originalData + testOffset, currTestSize_bytes);

            size_t currTrainPoints = totalNumPoints - pointsInFold[currFold];
            size_t currTrainSize = currTrainPoints * numFeatures;
            size_t currTrainSize_bytes = currTrainSize * sizeof(double);

            memcpy(currTrain, originalData, testOffset_bytes);
            memcpy(currTrain + testOffset, originalData + testOffset + currTestSize, currTrainSize_bytes - testOffset_bytes);

            finalResults[currFold - first_fold] = knn(currTrain, currTest, k, currTrainPoints, pointsInFold[currFold], numFeatures);
            testOffset += pointsInFold[currFold] * numFeatures;
        }
        double finalAccuracy = 0.0;
        for (int i = 0; i < numFolds; i++) {
            finalAccuracy += finalResults[i];
        }
        finalAccuracy /= numFolds;
        finalResults[numFolds] = finalAccuracy;
        writeResultsToFile(finalResults, 1, numFolds + 1, outFile);
        free(finalResults);
    }

    free(currTrain);
    free(currTest);
    free(results);
    free(pointsInFold);
    if (rank != 0) {
        free(originalData);
    }


    MPI_Finalize();
    return 0;
}

void printArray(void *array, char array_name[500], int m, int n, char type){
    printf("\n\n====Printing %s=====\n\n", array_name);

    int i, j;
    for(i = 0; i < m; i++){
        for(j = 0; j < n; j++){
            if(type == 'd') printf("%f, ", ((double*)array)[i * n + j]);
            else if (type == 'i') printf("%d, ", ((int*)array)[i * n + j]);

        }
        printf("\n");
    }
}