#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <immintrin.h> 

#define CHUNK_SIZE 100  // Define chunk size to process large datasets
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>

/*This code is for reading and writing to files for the 2024-25 COMP528 CA1*/

/*Use the functions in this file to read from the input file, and write to the output file*/

/*You should use this file when compiling your code*/

/*Declare these functions at the top of each 'main' file*/

/*If there are any issues with this code, please contact: h.j.forbes@liverpool.ac.uk*/

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
		fprintf(file, "\n\0");
        
    }

	fclose(file);

	//if not null not returned, then it is pointer to output meaning success.
	return output;
}

// used for storing distances
typedef struct {
    double value;  
    int index;     
    double class;  
} PointDistance;

// relevant enough to have its own function
int compare(const void *a, const void *b) {
    if (((PointDistance*)a)->value > ((PointDistance*)b)->value) return 1;
    else if (((PointDistance*)a)->value < ((PointDistance*)b)->value) return -1;
    else return 0;
}

// Swap helper function
inline void swap(PointDistance* a, PointDistance* b) {
    PointDistance temp = *a;
    *a = *b;
    *b = temp;
}

// Partition function for Quickselect
inline int partition(PointDistance arr[], int left, int right) {
    double pivot = arr[right].value;
    int i = left;

    for (int j = left; j < right; j++) {
        if (arr[j].value < pivot) {
            swap(&arr[i], &arr[j]);
            i++;
        }
    }
    swap(&arr[i], &arr[right]);
    return i;
}

// Quickselect function to find k-th smallest element
void quickselect(PointDistance arr[], int left, int right, int k) {
    if (left < right) {
        int pivotIndex = partition(arr, left, right);
        if (pivotIndex == k) {
            return;  // Found the k-th element
        } else if (pivotIndex > k) {
            quickselect(arr, left, pivotIndex - 1, k);
        } else {
            quickselect(arr, pivotIndex + 1, right, k);
        }
    }
}

// this probably saves okay time, i'll test it at some point
// rather than sorting all the points, only the k smallest points are sorted
void partial_sort(PointDistance arr[], int n, int k) {
    quickselect(arr, 0, n - 1, k); 
    qsort(arr, k, sizeof(PointDistance), compare); 
}

double findMostFrequentWithTieBreak(PointDistance arr[], int k) {
    // if a test case has more than 100 classes then i hope your pillow is warm tonight
    int classCount[100] = {0}; 
    int max_count = 0;
    double most_frequent_class = arr[0].class;
    bool tie_occurred = false;

    // pretty self explanitory but since i dont want to lose marks it counts the most frequent class and tracks if a tie is occuring or not
    for (int i = 0; i < k; i++) {
        classCount[(int)arr[i].class]++;
        if (classCount[(int)arr[i].class] > max_count) {
            max_count = classCount[(int)arr[i].class];
            most_frequent_class = arr[i].class;
            tie_occurred = false;
        } else if (classCount[(int)arr[i].class] == max_count) {
            tie_occurred = true; 
        }
    }
    // because the array is sorted arr[0] is the closest point
    if (tie_occurred) {
        return arr[0].class;
    }
    return most_frequent_class;
}

void processChunk(double *train_data, double *test_data, int train_rows, int test_rows, int train_cols, int test_cols, int k, int chunk_start, int chunk_size) {
    PointDistance *distances = (PointDistance*) malloc(train_rows * sizeof(PointDistance)); 

    for (int i = chunk_start; i < chunk_start + chunk_size && i < test_rows; i++) {
        for (int j = 0; j < train_rows; j++) {
            double dist = 0.0;
            // v_sum = [0,0,0,0]
            __m256d v_sum = _mm256_setzero_pd();

            int d = 0;
            // Vectorized distance calculation for 4 dimensions at a time
            for (; d <= test_cols - 5; d += 4) {
                // v_test = testing point i feature [i + d] -> [i + d + 4]
                __m256d v_test = _mm256_loadu_pd(&test_data[i * test_cols + d]);
                // v_train = training point j features [j + d] -> [j + d + 4]
                __m256d v_train = _mm256_loadu_pd(&train_data[j * train_cols + d]);
                // v_diff = v_test - v_train
                __m256d v_diff = _mm256_sub_pd(v_test, v_train);
                // v_sum = v_diff * v_diff
                v_sum = _mm256_add_pd(v_sum, _mm256_mul_pd(v_diff, v_diff));
            }

            // dist = squared distance between test_point[i] and train_point[j]
            dist = v_sum[0] + v_sum[1] + v_sum[2] + v_sum[3];

            // if point ∈ R^n where n%4 != 0 then this handles the remaining points
            for (; d < test_cols - 1; d++) {
                double diff = test_data[i * test_cols + d] - train_data[j * train_cols + d];
                dist += diff * diff;
            }
            // update distance to Point j in array
            distances[j].value = dist;
            distances[j].index = j;
            distances[j].class = train_data[j * train_cols + (train_cols - 1)];
        }

        // distances [0-k] are the k smallest in order
        partial_sort(distances, train_rows, k);

        // write the correct class to the array of test points
        // NOTE - this can be done in parallel because of each point being completely independant in memory
        test_data[i * test_cols + (test_cols - 1)] = findMostFrequentWithTieBreak(distances, k);
    }

    free(distances); 
}

int main(int argc, char *argv[]) {
    // if this bit needs explaining thats not my problem
    int train_rows = readNumOfPoints(argv[1]);
    int train_cols = readNumOfFeatures(argv[1]);
    double *train_data = readDataPoints(argv[1], train_rows, train_cols);

    int test_rows = readNumOfPoints(argv[2]);
    int test_cols = readNumOfFeatures(argv[2]);
    double *test_data = readDataPoints(argv[2], test_rows, test_cols);

    char *outfile = argv[3];
    int k = atoi(argv[4]);

    // todo: find better sizing thats dynamic with number of points
    int chunk_size = CHUNK_SIZE;

    /*
    process each chunk in parallel, this was my first approach and i am yet to find a better one
    as ALL of the processing for each chunk of points can be done completely independantly, probably missing
    some slight efficency by not parallelising every point but the chunking made working with memory easier.
    */  
    #pragma omp parallel for schedule(dynamic)
    for (int chunk_start = 0; chunk_start < test_rows; chunk_start += chunk_size) {
        int current_chunk_size = (chunk_start + chunk_size > test_rows) ? (test_rows - chunk_start) : chunk_size;
        processChunk(train_data, test_data, train_rows, test_rows, train_cols, test_cols, k, chunk_start, current_chunk_size);
    }

    writeResultsToFile(test_data, test_rows, test_cols, outfile);

    // no memory leaks today
    free(train_data);
    free(test_data);

    return 0;
}
