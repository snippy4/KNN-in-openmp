#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>


int readNumOfPoints(char*);
int readNumOfFeatures(char*);
int readNumOfClasses(char *filename);
double *readDataPoints(char*, int, int);
void *writeResultsToFile(double*, int, int, char*);

/*Gets the number of the coordinates in the file. Returns as a single integer*/
int readNumOfPoints(char *filename){
	FILE *file = fopen(filename, "r");
	int numOfPoints = 0;

	if(file == NULL){
		return -1;
	}

	char line[1000];

	printf("Reading num of points: ");

	while(fgets(line, sizeof(line), file) != NULL){
		numOfPoints++;
	}
	printf("%d\n", numOfPoints);

	fclose(file);

    return numOfPoints;
}

int readNumOfFeatures(char *filename){
	FILE *file = fopen(filename, "r");

	if(file == NULL){
		return -1;
	}

	char line[1000];

	char *p;

	printf("Reading num of features: ");
	if(fgets(line, sizeof(line), file) != NULL){
		int numOfFeatures = 1;

		for(p = line; *p != '\0'; p++){
			if(*p == ','){
				numOfFeatures++;
			}
		}

		fclose(file);
		printf("%d\n", numOfFeatures);
		return numOfFeatures;
	}

	fclose(file);
	return -1;
}

int readNumOfClasses(char *filename){
	FILE *file = fopen(filename, "r");

	if(file == NULL){
		return -1;
	}

	char line[1000];

	int numOfClasses = 0;

	printf("Reading num of classes: ");
	while(fgets(line, sizeof(line), file) != NULL){

		char *lastcomma = strrchr(line, ',');
		char *token = strtok(lastcomma, ",");

		int class = atoi(token);

		if(class > numOfClasses){
			numOfClasses = class;
		}

	}

	printf("%d\n", numOfClasses + 1);

	fclose(file);
	return numOfClasses + 1;
}
/*Gets the data from the file. Returns as a 2D array of doubles. 
The program returns all the features, and the final element is the class label*/
double *readDataPoints(char *filename, int numOfPoints, int numOfFeatures){
	FILE *file = fopen(filename,"r");
    int i;

	char line[1000];
    
    if(file == NULL) {
        printf("Unable to open file: %s\n", filename);
        return NULL;
    }

	double *dataPoints = (double *)malloc(numOfPoints * numOfFeatures * sizeof(double *));

	if(dataPoints == NULL){
		return NULL;
	}

	int lineNum = 0;

	int featureIndex;


	printf("Reading data points: ");
	while(fgets(line, sizeof(line), file) != NULL && lineNum < numOfPoints){
		char* token = strtok(line, ",");
		featureIndex = 0;

		while(token != NULL && featureIndex < numOfFeatures){
			
			dataPoints[lineNum * numOfFeatures + featureIndex] = atof(token);
			//if(lineNum<10) printf("LineNum: %d, NumFeatures: %d, Value = %f\n", lineNum, numOfFeatures, atof(token));
			token = strtok(NULL, ",");
			featureIndex++;
		}

		//printf("\n");

		lineNum++;

		/*
		if (sscanf(line, "%lf,%lf", &x, &y) == 2){
			coords[lineNum][0] = x;
			coords[lineNum][1] = y;
			lineNum++;
		}*/
	}

	printf("%d\n", lineNum + 1);

	fclose(file);

	return dataPoints;
}

void *writeResultsToFile(double *output, int numOfPoints, int numOfFeatures, char *filename){
	
	FILE *file = fopen(filename, "w");
	int i, j;	
	
	if(file == NULL){
		printf("Unable to open file: %s", filename);
		return NULL;
	}

	printf("Writing output data to file: %s:\nnumber of data points: %d\nnumber of features: %d\n", filename, numOfPoints, numOfFeatures);
    for(i=0; i < numOfPoints; i++) {
		for(j = 0; j < numOfFeatures; j++){
			if(j < numOfFeatures - 1) fprintf(file, "%lf,", output[i * numOfFeatures + j]);
			else fprintf(file, "%lf", output[i * numOfFeatures + j]);
		}
		fprintf(file, "\n");
        
    }

	fclose(file);
	return output;
}

typedef struct {
    double value;
    int index;
} ValueIndexPair;

int compare(const void *a, const void *b) {
    if (((ValueIndexPair*)a)->value > ((ValueIndexPair*)b)->value) return 1;
    else if (((ValueIndexPair*)a)->value < ((ValueIndexPair*)b)->value) return -1;
    else return 0;
}

void findKLowestInSubsection(double arr[], int total_size, int start_index, int n, int k, int result_indexes[]) {
    if (start_index + n > total_size) {
        n = total_size - start_index;
    }
    ValueIndexPair *subsection = (ValueIndexPair*) malloc(n * sizeof(ValueIndexPair));
    for (int i = 0; i < n; i++) {
        subsection[i].value = arr[start_index + i];
        subsection[i].index = start_index + i;
    }
    qsort(subsection, n, sizeof(ValueIndexPair), compare);
    for (int i = 0; i < k && i < n; i++) {
        result_indexes[i] = subsection[i].index;
    }
    free(subsection);
}

double findMostFrequent(double arr[], int n) {
    int max_count = 0;      
    double most_frequent = arr[0]; 
    for (int i = 0; i < n; i++) {
        int count = 0;  
        for (int j = 0; j < n; j++) {
            if (arr[j] == arr[i]) {
                count++;
            }
        }
        if (count > max_count) {
            max_count = count;
            most_frequent = arr[i];
        }
    }
    return most_frequent;
}
/**
 *this will containt the basic sequential solution to the problem
 *the goal of this program is to run and produce correct output, not to be optimised.
*/
int main(int argc, char *argv[]){
    clock_t time = clock();
    int train_rows = readNumOfPoints(argv[1]);
    int train_cols = readNumOfFeatures(argv[1]);
    double *train_data = readDataPoints(argv[1], train_rows, train_cols);
    int test_rows = readNumOfPoints(argv[2]);
    int test_cols = readNumOfFeatures(argv[2]);
    double *test_data = readDataPoints(argv[2], test_rows, test_cols);
    char *outfile = argv[3];
	int k = atoi(argv[4]);
	double dist = 0;
	double *point_distances;
	point_distances = malloc(train_rows*test_rows*sizeof(double));
	if(point_distances == NULL){
		printf("memory fail\n");
		exit(1);
	}
	printf("memory allocated\n");
	printf("max i %d\n", test_rows*test_cols);
	#pragma omp parallel for
	for (int i = 0; i < test_rows; i++) { 
        for (int j = 0; j < train_rows; j++) { 
            double dist = 0.0;
            for (int d = 0; d < test_cols-1; d++) {
                double diff = test_data[i * test_cols + d] - train_data[j * train_cols + d];
                dist += diff * diff;
            }
            point_distances[i * train_rows + j] = dist; 
        }
    }
	int k_lowest[k];
	double class[k];
	printf("distances calculated\n");
	//#pragma omp parallel for
	for(int i = 0; i<train_rows*test_rows; i+=train_rows){
		findKLowestInSubsection(point_distances, train_rows*test_rows, i, train_rows, k, k_lowest);
		for(int j = 0; j<k; j++){
			class[j] = train_data[(k_lowest[j]%train_rows)*train_cols + train_cols-1];
		}
		test_data[(i/train_rows+1)*test_cols-1] = findMostFrequent(class, k);
		for(int j = 0; j<k; j++){
		}
	}
	printf("classes calculated\n");

	writeResultsToFile(test_data, test_rows, test_cols, outfile);
    
    time = ( clock() - time);
    //printf("%f", (float) time / CLOCKS_PER_SEC );
	free(point_distances);
    free(train_data);
    free(test_data);
    return 0;
}