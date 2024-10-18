#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <immintrin.h> 

#define CHUNK_SIZE 10000  // Define chunk size to process large datasets
typedef struct {
    double value;  // Distance
    int index;     // Index of the training point
    double class;  // Class of the training point
} ValueIndexPair;

int compare(const void *a, const void *b) {
    if (((ValueIndexPair*)a)->value > ((ValueIndexPair*)b)->value) return 1;
    else if (((ValueIndexPair*)a)->value < ((ValueIndexPair*)b)->value) return -1;
    else return 0;
}

ValueIndexPair *pre_alloc_distances;

int readNumOfPoints(char*);
int readNumOfFeatures(char*);
double *readDataPoints(char*, int, int);
void *writeResultsToFile(double*, int, int, char*);
// Swap helper function
void swap(ValueIndexPair* a, ValueIndexPair* b) {
    ValueIndexPair temp = *a;
    *a = *b;
    *b = temp;
}

// Partition function for Quickselect
int partition(ValueIndexPair arr[], int left, int right) {
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
void quickselect(ValueIndexPair arr[], int left, int right, int k) {
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

// Partial sort to get the k smallest elements, followed by sorting those k elements
void partial_sort(ValueIndexPair arr[], int n, int k) {
    quickselect(arr, 0, n - 1, k);  // Rearrange the first k elements to be the smallest
    qsort(arr, k, sizeof(ValueIndexPair), compare);  // Sort only the k smallest elements
}
/*Gets the number of the coordinates in the file. Returns as a single integer*/
int readNumOfPoints(char *filename){
	FILE *file = fopen(filename, "r");
	int numOfPoints = 0;

	if(file == NULL){
		return -1;
	}

	char line[1000];

	while(fgets(line, sizeof(line), file) != NULL){
		numOfPoints++;
	}

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

	if(fgets(line, sizeof(line), file) != NULL){
		int numOfFeatures = 1;

		for(p = line; *p != '\0'; p++){
			if(*p == ','){
				numOfFeatures++;
			}
		}

		fclose(file);
		return numOfFeatures;
	}

	fclose(file);
	return -1;
}

double *readDataPoints(char *filename, int numOfPoints, int numOfFeatures){
	FILE *file = fopen(filename,"r");
	char line[1000];
    
    if(file == NULL) {
        printf("Unable to open file: %s\n", filename);
        return NULL;
    }

	double *dataPoints = (double *)malloc(numOfPoints * numOfFeatures * sizeof(double));
	if(dataPoints == NULL){
		return NULL;
	}

	int lineNum = 0, featureIndex;

	while(fgets(line, sizeof(line), file) != NULL && lineNum < numOfPoints){
		char* token = strtok(line, ",");
		featureIndex = 0;

		while(token != NULL && featureIndex < numOfFeatures){
			dataPoints[lineNum * numOfFeatures + featureIndex] = atof(token);
			token = strtok(NULL, ",");
			featureIndex++;
		}

		lineNum++;
	}

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

	for(i = 0; i < numOfPoints; i++) {
		for(j = 0; j < numOfFeatures; j++){
			if(j < numOfFeatures - 1) fprintf(file, "%lf,", output[i * numOfFeatures + j]);
			else fprintf(file, "%lf", output[i * numOfFeatures + j]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
	return output;
}


double findMostFrequentWithTieBreak(ValueIndexPair arr[], int k) {
    int classCount[100] = {0};  // Assuming class labels are between 0 and 99
    int max_count = 0;
    double most_frequent_class = arr[0].class;
    int tie_occurred = 0;

    for (int i = 0; i < k; i++) {
        classCount[(int)arr[i].class]++;
        if (classCount[(int)arr[i].class] > max_count) {
            max_count = classCount[(int)arr[i].class];
            most_frequent_class = arr[i].class;
            tie_occurred = 0;
        } else if (classCount[(int)arr[i].class] == max_count) {
            tie_occurred = 1;  // Indicate a tie occurred
        }
    }

    // If a tie occurs, we choose the class of the closest point (first one in the array after sorting)
    if (tie_occurred) {
        return arr[0].class;
    }
    return most_frequent_class;
}

void processChunk(double *train_data, double *test_data, int train_rows, int test_rows, int train_cols, int test_cols, double *point_distances, int k, int chunk_start, int chunk_size) {
    ValueIndexPair *distances = pre_alloc_distances;

    for (int i = chunk_start; i < chunk_start + chunk_size && i < test_rows; i++) {
        for (int j = 0; j < train_rows; j++) {
            double dist = 0.0;
            // Vectorized distance calculation (example with SIMD, depends on your CPU)
            __m256d v_sum = _mm256_setzero_pd();
            for (int d = 0; d < test_cols - 1; d += 4) {
                __m256d v_test = _mm256_loadu_pd(&test_data[i * test_cols + d]);
                __m256d v_train = _mm256_loadu_pd(&train_data[j * train_cols + d]);
                __m256d v_diff = _mm256_sub_pd(v_test, v_train);
                v_sum = _mm256_add_pd(v_sum, _mm256_mul_pd(v_diff, v_diff));
            }
            dist = v_sum[0] + v_sum[1] + v_sum[2] + v_sum[3];  // Sum the vector components

            // Store distance and class info
            distances[j].value = dist;
            distances[j].index = j;
            distances[j].class = train_data[j * train_cols + (train_cols - 1)];
        }

        // Partial sort (for k-nearest neighbors)
        partial_sort(distances, train_rows, k);  // Replace qsort with partial sort

        // Assign the most frequent class with tie-breaking
        test_data[i * test_cols + (test_cols - 1)] = findMostFrequentWithTieBreak(distances, k, i);
    }
}

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
    pre_alloc_distances = (ValueIndexPair*) malloc(train_rows * sizeof(ValueIndexPair));

	// Chunk processing
	int chunk_size = CHUNK_SIZE;
	double *point_distances = malloc(chunk_size * train_rows * sizeof(double));

	if (point_distances == NULL) {
		printf("Memory allocation failed for point_distances\n");
		exit(1);
	}
	#pragma omp parallel for schedule(dynamic)
	for (int chunk_start = 0; chunk_start < test_rows; chunk_start += chunk_size) {
		int current_chunk_size = (chunk_start + chunk_size > test_rows) ? (test_rows - chunk_start) : chunk_size;
		processChunk(train_data, test_data, train_rows, test_rows, train_cols, test_cols, point_distances, k, chunk_start, current_chunk_size);
	}

	writeResultsToFile(test_data, test_rows, test_cols, outfile);

	free(point_distances);
	free(train_data);
	free(test_data);
    free(pre_alloc_distances);

	time = (clock() - time);
	printf("Total time: %f seconds\n", (float)time / CLOCKS_PER_SEC);

	return 0;
}
