#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include"file-reader.c"
#include"knnomp.c"

void printArray(void *array, char array_name[500], int m, int n, char type);

int main(int argc, char *argv[]){

    printf("\n\n===============STARTING KNN===============\n\n");   

    //File reading
    char *inFile = argv[1];
    char *outFile = argv[2];
    int k = atoi(argv[3]);
    int numFolds = atoi(argv[4]);

    int totalNumPoints, numFeatures, numClasses;
    double *originalData;
    
    totalNumPoints = readNumOfPoints(inFile);
    numFeatures = readNumOfFeatures(inFile);
    numClasses = readNumOfClasses(inFile);
    originalData = readDataPoints(inFile, totalNumPoints, numFeatures);

    //Array contains the number of points in each fold
    int *pointsInFold = (int *)malloc(numFolds * sizeof(int));
    
    //Calculate the points in the fold and the remainder
    int pointsPerFold = totalNumPoints / numFolds;
    int pointsPerFoldRemainder = totalNumPoints % numFolds;

    //Allocates the number of points in the fold adding the remainder to the last fold. (Can handle unequally sized folds)
    int fold;
    for(fold = 0; fold < numFolds; fold++){
        pointsInFold[fold] = pointsPerFold + ((fold == numFolds - 1) ? pointsPerFoldRemainder : 0);
    }

    /**
     * Calculate the maximum size of the train and test folds
     * Since the last fold always has the amount, we use this to our advantage
    **/
    size_t maxTrain_bytes = (totalNumPoints - pointsInFold[0]) * numFeatures * sizeof(double);
    size_t maxTest_bytes = pointsInFold[numFolds - 1] * numFeatures * sizeof(double);

    //Allocate the memory for the train and test datasets
    double *currTrain = (double*)malloc(maxTrain_bytes);
    double *currTest = (double*)malloc(maxTest_bytes);
    double *results = (double*)malloc((numFolds+1) * sizeof(double));

    //Declare memcpy variables 
    size_t currTestSize, currTestSize_bytes, currTrainPoints, currTrainSize, currTrainSize_bytes, testOffset_bytes;

    //Test offset is used to track where we are in the originalData
    size_t testOffset = 0;
    
    int currFold;
    for(currFold = 0; currFold < numFolds; currFold++){
        
        currTestSize = pointsInFold[currFold] * numFeatures;
        currTestSize_bytes = currTestSize * sizeof(double);
        testOffset_bytes = testOffset * sizeof(double);

        //Copy the current fold into the currTest
        memcpy(currTest, originalData + testOffset, currTestSize_bytes);

        currTrainPoints = totalNumPoints - pointsInFold[currFold];
        currTrainSize = currTrainPoints * numFeatures;
        currTrainSize_bytes = currTrainSize * sizeof(double);

        //Copy the remaining folds into currTrain
        memcpy(currTrain, originalData, testOffset_bytes);
        memcpy(currTrain + testOffset, originalData + testOffset + currTestSize, currTrainSize_bytes - testOffset_bytes);

        //For debugging. (Don't print asteroids)
        //printArray(currTest, "test", pointsInFold[currFold], numFeatures, 'd');
       //printArray(currTrain, "train", currTrainPoints, numFeatures, 'd');
        
       double accuracy = knn(currTrain, currTest, k, currTrainPoints, pointsInFold[currFold], numFeatures);
       results[currFold] = accuracy;


        testOffset += pointsInFold[currFold] * numFeatures;
    }   
    double final = 0.0;
    for(int i = 0; i < numFolds; i++){
        final += results[i];
    }
    final = final / (double)numFolds;
    results[numFolds] = final;
    writeResultsToFile(results, 1, numFolds + 1, outFile);
    free(currTrain);
    free(currTest);
}

/**
 * Prints a given array
 * @param array: the array being given
 * @param array_name: the name of the array
 * @param m: the number of elements for the m dimension
 * @param n: the number of elements for the n dimension
 * @param type: the type of the array. e.g. d=double i=integer.
 **/
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