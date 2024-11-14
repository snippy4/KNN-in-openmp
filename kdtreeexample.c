#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
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
typedef struct KDNode {
    double *point;       // Pointer to point in space
    int index;           // Original index of point in dataset
    int class;           // Class label of point
    struct KDNode *left; // Left child
    struct KDNode *right;// Right child
} KDNode;
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


// Recursive function to build k-d tree
KDNode* buildKDTree(double *data, int start, int end, int depth, int num_features) {
    if (start > end) return NULL;

    int axis = depth % num_features;
    int median = (start + end) / 2;

    // Sort points by current axis and select the median
    qsort(data + start * num_features, end - start + 1, num_features * sizeof(double), 
          (int (*)(const void*, const void*))compareByAxis(axis));

    KDNode *node = malloc(sizeof(KDNode));
    node->point = data + median * num_features;
    node->class = node->point[num_features - 1];
    node->index = median;

    // Recursively build left and right subtrees
    node->left = buildKDTree(data, start, median - 1, depth + 1, num_features);
    node->right = buildKDTree(data, median + 1, end, depth + 1, num_features);

    return node;
}

// Global variable to hold the current axis for sorting
int current_axis;

int compareByAxis(const void *a, const void *b) {
    const double *pointA = *(const double **)a;
    const double *pointB = *(const double **)b;
    
    // Compare based on the specified axis
    if (pointA[current_axis] < pointB[current_axis]) {
        return -1;
    } else if (pointA[current_axis] > pointB[current_axis]) {
        return 1;
    } else {
        return 0;
    }
}

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

// Build a max-heap with the first k elements
void buildMaxHeap(PointDistance arr[], int k) {
    for (int i = (k / 2) - 1; i >= 0; i--) {
        maxHeapify(arr, k, i);
    }
}

// Function to find k smallest elements using a max-heap of size k
void findKSmallestElements(PointDistance arr[], int n, int k) {
    // Step 1: Build a max-heap with the first k elements
    buildMaxHeap(arr, k);

    // Step 2: Process the remaining elements
    for (int i = k; i < n; i++) {
        if (arr[i].value < arr[0].value) {
            // Replace the root (maximum element) with the current element and re-heapify
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
    // because the array is sorted arr[0] is the closest point
    if (tie_occurred) {
        qsort(arr, k, sizeof(PointDistance), compare);
        return arr[0].class;
    }
    return most_frequent_class;
}

void kNN_Search(KDNode *node, double *target, int depth, PointDistance *neighbors, int *neighbor_count, int k, int num_features) {
    if (node == NULL) return;

    double dist = 0.0;
    for (int i = 0; i < num_features - 1; i++) {
        dist += (target[i] - node->point[i]) * (target[i] - node->point[i]);
    }

    // Insert into max-heap of neighbors if distance is smaller than the current max
    if (*neighbor_count < k) {
        neighbors[*neighbor_count].value = dist;
        neighbors[*neighbor_count].index = node->index;
        neighbors[*neighbor_count].class = node->class;
        (*neighbor_count)++;
        if (*neighbor_count == k) buildMaxHeap(neighbors, k);
    } else if (dist < neighbors[0].value) {
        neighbors[0] = (PointDistance){ dist, node->index, node->class };
        maxHeapify(neighbors, k, 0);
    }

    int axis = depth % (num_features - 1);
    KDNode *nearBranch = (target[axis] < node->point[axis]) ? node->left : node->right;
    KDNode *farBranch = (nearBranch == node->left) ? node->right : node->left;

    // Traverse the near branch
    kNN_Search(nearBranch, target, depth + 1, neighbors, neighbor_count, k, num_features);

    // Backtrack and check the far branch if necessary
    if (fabs(target[axis] - node->point[axis]) < neighbors[0].value) {
        kNN_Search(farBranch, target, depth + 1, neighbors, neighbor_count, k, num_features);
    }
}


// Sample main function
int main(int argc, char *argv[]){
    // if this bit needs explaining thats not my problem
    int train_rows = readNumOfPoints(argv[1]);
    int train_cols = readNumOfFeatures(argv[1]);
    double *train_data = readDataPoints(argv[1], train_rows, train_cols);

    int test_rows = readNumOfPoints(argv[2]);
    int test_cols = readNumOfFeatures(argv[2]);
    double *test_data = readDataPoints(argv[2], test_rows, test_cols);

    char *outfile = argv[3];
    int k = atoi(argv[4]);

    // Build k-d tree from training data
    KDNode *kd_tree_root = buildKDTree(train_data, 0, train_rows - 1, 0, train_cols);

    // Array to store distances to nearest neighbors
    PointDistance *neighbors = (PointDistance *) malloc(k * sizeof(PointDistance));

    // For each test point, find the k-nearest neighbors
    for (int i = 0; i < test_rows; i++) {
        int neighbor_count = 0; // Track the number of neighbors found

        // Search for k-nearest neighbors using the k-d tree
        kNN_Search(kd_tree_root, &test_data[i * test_cols], 0, neighbors, &neighbor_count, k, test_cols);

        // Determine the most frequent class among the k-nearest neighbors
        test_data[i * test_cols + (test_cols - 1)] = findMostFrequentWithTieBreak(neighbors, k);
    }
    writeResultsToFile(test_data, test_rows, test_cols, outfile);


    return 0;
}
