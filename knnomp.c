#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <immintrin.h> 

#define CHUNK_SIZE 40  // Define chunk size to process large datasets

/*This code is for reading and writing to files for the 2024-25 COMP528 CA1*/

/*Use the functions in this file to read from the input file, and write to the output file*/

/*You should use this file when compiling your code*/

/*Declare these functions at the top of each 'main' file*/

/*If there are any issues with this code, please contact: h.j.forbes@liverpool.ac.uk*/

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

void swap(PointDistance* a, PointDistance* b) {
    PointDistance temp = *a;
    *a = *b;
    *b = temp;
}


//online algorithm for inserting into max heap
static inline void maxHeapify(PointDistance arr[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    while (left < n) {
        if (arr[left].value > arr[largest].value) {
            largest = left;
        }
        if (right < n && arr[right].value > arr[largest].value) {
            largest = right;
        }
        //if the largest element has changed then swap arr[i] and largest
        if (largest != i) {
            swap(&arr[i], &arr[largest]);
            i = largest;
            left = 2 * i + 1;
            right = 2 * i + 2;
        } else {
            break;
        }
    }
}

//arr[0..k] is entered into the heap
void buildMaxHeap(PointDistance arr[], int k) {
    for (int i = (k / 2) - 1; i >= 0; i--) {
        maxHeapify(arr, k, i);
    }
}

// find k smallest
void findKSmallestElements(PointDistance arr[], int n, int k) {
    buildMaxHeap(arr, k);

    //if given element is smaller than the largest element in the max heap then insert the new element into the heap
    for (int i = k; i < n; i++) {
        if (arr[i].value < arr[0].value) {
            arr[0] = arr[i];
            maxHeapify(arr, k, 0);
        }
    }
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
    // select the closest point
    if (tie_occurred) {
        qsort(arr, k, sizeof(PointDistance), compare);
        return arr[0].class;
    }
    return most_frequent_class;
}

void processChunk(double *train_data, double *test_data, int train_rows, int test_rows, int train_cols, int test_cols, int k, int chunk_start, int chunk_size) {
    // create a max-heap of size k
    PointDistance *distances = (PointDistance*) malloc(k * sizeof(PointDistance)); 
    int end = (chunk_start + chunk_size > test_rows) ? test_rows : (chunk_start + chunk_size);

    for (int i = chunk_start; i < end; i++) {
        // initialize the heap
        for (int n = 0; n < k; n++) {
            distances[n].value = INFINITY; 
        }

        // calculate distances from the test point to each training point 
        for (int j = 0; j < train_rows; j++) {
            __builtin_prefetch(&train_data[(j + 1) * train_cols], 0, 1);

            // before calculating the full distance, check if pruning is possible
            double current_max_dist = distances[0].value; 
            double dist = 0.0;
            bool prune = false;
           

            // v_sum = [0,0,0,0]
            __m256d v_sum = _mm256_setzero_pd();
            // d = 1 to account for id
            int d = 1;

            // vectorized distance calculation for 4 dimensions at a time
            for (; d <= test_cols - 5; d += 4) {
                // v_test = testing point i feature [i + d] -> [i + d + 4]
                __m256d v_test = _mm256_loadu_pd(&test_data[i * test_cols + d]);
                // v_train = training point j features [j + d] -> [j + d + 4]
                __m256d v_train = _mm256_loadu_pd(&train_data[j * train_cols + d]);
                // v_diff = v_test - v_train
                __m256d v_diff = _mm256_sub_pd(v_test, v_train);
                // v_sum = v_diff * v_diff
                v_sum = _mm256_add_pd(v_sum, _mm256_mul_pd(v_diff, v_diff));
            // prune the point if the current distance is larger than the current k-th closest point
                if (v_sum[0] + v_sum[1] + v_sum[2] + v_sum[3] > current_max_dist) {
                    prune=true;
                    break; 
            }
            }
            if(prune){
                continue;
            }
            // dist = squared distance between test_point[i] and train_point[j]
            dist = v_sum[0] + v_sum[1] + v_sum[2] + v_sum[3];

            // If point ∈ R^n where n%4 != 0 then this handles the remaining points
            for (; d < test_cols - 1; d++) {
                double diff = test_data[i * test_cols + d] - train_data[j * train_cols + d];
                dist += diff * diff;
                // prune the point if the current distance is larger than the current k-th closest point
                if (dist > current_max_dist) {
                    prune=true;
                    break; 
            }
            }
            if(prune){
                continue;
            }
            
            //insert point into heap if it is closer than the k-th closest
            if (dist < current_max_dist) {
                distances[0].value = dist;
                distances[0].index = j;
                distances[0].class = train_data[j * train_cols + (train_cols - 1)];

                // heapify the heap to maintain the k smallest distances
                maxHeapify(distances, k, 0);
            }
        }

        // find the most frequent class
        test_data[i * test_cols + (test_cols - 1)] = findMostFrequentWithTieBreak(distances, k);
    }

    free(distances); 
}


double knn(double *train_data, double *test_data, int k, int train_rows, int test_rows, int features){
    int chunk_size = 10;
    double *points = (double *)malloc(sizeof(double) * train_rows);
    for(int i = 0; i<test_rows;i++){
        points[i] = test_data[i * features + (features - 1)];
    }

    #pragma omp parallel for schedule(dynamic)
    for (int chunk_start = 0; chunk_start < test_rows; chunk_start += chunk_size) {
        int current_chunk_size = (chunk_start + chunk_size > test_rows) ? (test_rows - chunk_start) : chunk_size;
        processChunk(train_data, test_data, train_rows, test_rows, features, features, k, chunk_start, current_chunk_size);
    }

    int differences = 0;

    #pragma omp parallel for reduction(+:differences)
    for (int i = 0; i < test_rows; i++) {
        if (points[i] != test_data[i * features + (features - 1)]) {
            differences++;
        }
    }

    double acc;
    acc = ((double)test_rows - (double)differences) / (double)test_rows;
    free(points);
    return acc;
}
