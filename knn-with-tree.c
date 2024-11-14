#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include<string.h>
#include<stdbool.h>
#include<omp.h>

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

// Define the KD-Tree node structure
typedef struct KDNode {
    double* point;        // Point in n-dimensional space
    struct KDNode *left;  // Left subtree
    struct KDNode *right; // Right subtree
    double class;
} KDNode;

// Global variable to hold the number of dimensions (n)
int dimensions;

// Function to create a new KD-Tree node
KDNode* createNode(double point[], double class) {
    KDNode* node = (KDNode*)malloc(sizeof(KDNode));
    node->point = (double*)malloc(dimensions * sizeof(double));
    for (int i = 0; i < dimensions; i++)
        node->point[i] = point[i];
    node->left = node->right = NULL;
    node->class = class;
    return node;
}

// Insert a new point into the KD-Tree
KDNode* insert(KDNode* root, double point[], int depth, double class) {
    if (root == NULL) return createNode(point, class);
    // Determine current axis
    int axis = depth % dimensions;

    if (point[axis] < root->point[axis])
        root->left = insert(root->left, point, depth + 1, class);
    else
        root->right = insert(root->right, point, depth + 1, class);

    return root;
}

// Helper function to calculate Euclidean distance between two points
double distance(double point1[], double point2[]) {
    double dist = 0.0;
    for (int i = 0; i < dimensions; i++)
        dist += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    return dist;
}

// A helper struct to hold nearest neighbors
typedef struct {
    KDNode* node;
    double distance;
} Neighbor;
// Utility function to print a point

typedef struct {
    KDNode* closest_node;
    double closest_distance;
} ClosestNeighbor;

void printPoint(double point[]) {
    printf("(");
    for (int i = 0; i < dimensions; i++) {
        printf("%lf", point[i]);
        if (i < dimensions - 1) printf(", ");
    }
    printf(")");
}

// Quick helper to sort neighbors
int compareNeighbors(const void* a, const void* b) {
    Neighbor* n1 = (Neighbor*)a;
    Neighbor* n2 = (Neighbor*)b;
    if (n1->distance < n2->distance) return -1;
    else if (n1->distance > n2->distance) return 1;
    else return 0;
}


// Max-heapify function for neighbors array
void maxHeapify(Neighbor neighbors[], int k, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < k && neighbors[left].distance > neighbors[largest].distance)
        largest = left;

    if (right < k && neighbors[right].distance > neighbors[largest].distance)
        largest = right;

    if (largest != i) {
        Neighbor temp = neighbors[i];
        neighbors[i] = neighbors[largest];
        neighbors[largest] = temp;
        maxHeapify(neighbors, k, largest);
    }
}

// Function to insert a new neighbor in the max-heap if it's closer
void insertNeighbor(Neighbor neighbors[], KDNode* node, double dist, int* neighbor_count, int k) {
    if (*neighbor_count < k) {
        // Add directly if fewer than k neighbors
        neighbors[*neighbor_count].node = node;
        neighbors[*neighbor_count].distance = dist;
        (*neighbor_count)++;

        if (*neighbor_count == k) {
            // Build the heap once we've added k elements
            for (int i = k / 2 - 1; i >= 0; i--)
                maxHeapify(neighbors, k, i);
        }
    } else if (dist < neighbors[0].distance) {
        // Replace the farthest neighbor if the new one is closer
        neighbors[0].node = node;
        neighbors[0].distance = dist;
       maxHeapify(neighbors, k, 0);
    }
}

void kNearestNeighbors(KDNode* root, double target[], int depth, Neighbor* neighbors, int* neighbor_count, int k) {
    if (root == NULL) return;

    // Calculate the distance from the target to the current root node
    double dist = distance(target, root->point);

    // Insert the current node into the neighbors array if there is space
    if (*neighbor_count < k) {
        neighbors[*neighbor_count].node = root;
        neighbors[*neighbor_count].distance = dist;
        (*neighbor_count)++;
    } else {
        // If the neighbors array is full, only insert if the current distance is smaller
        if (dist < neighbors[k - 1].distance) {
            neighbors[k - 1].node = root;
            neighbors[k - 1].distance = dist;
        }
    }

    // Sort the neighbors array based on distance (ascending)
    qsort(neighbors, *neighbor_count, sizeof(Neighbor), compareNeighbors);

    // Determine the axis at the current depth of recursion
    int axis = depth % dimensions;
    KDNode* nearBranch = NULL;
    KDNode* farBranch = NULL;

    // Decide which subtree to explore first based on the current axis
    if (target[axis] < root->point[axis]) {
        nearBranch = root->left;
        farBranch = root->right;
    } else {
        nearBranch = root->right;
        farBranch = root->left;
    }

    // Explore the near branch first
    kNearestNeighbors(nearBranch, target, depth + 1, neighbors, neighbor_count, k);

    // After exploring the near branch, check if the far branch needs to be explored
    // If the far branch could contain closer points, check that too
    if (fabs(target[axis] - root->point[axis]) < neighbors[k - 1].distance) {
        kNearestNeighbors(farBranch, target, depth + 1, neighbors, neighbor_count, k);
    }
}




double findMostFrequentWithTieBreak(Neighbor arr[], int k) {
    // if a test case has more than 100 classes then i hope your pillow is warm tonight
    int classCount[100] = {0}; 
    int max_count = 0;
    double most_frequent_class = arr[0].node->class;
    bool tie_occurred = false;

    // pretty self explanitory but since i dont want to lose marks it counts the most frequent class and tracks if a tie is occuring or not
    for (int i = 0; i < k; i++) {
        //printPoint(arr[i].node->point);
        //printf("dist: %lf class: %lf\n", arr[i].distance, arr[i].node->class);
        classCount[(int)arr[i].node->class]++;
        if (classCount[(int)arr[i].node->class] > max_count) {
            max_count = classCount[(int)arr[i].node->class];
            most_frequent_class = arr[i].node->class;
            tie_occurred = false;
        } else if (classCount[(int)arr[i].node->class] == max_count) {
            tie_occurred = true; 
        }
    }
    // because the array is sorted arr[0] is the closest point
    if (tie_occurred) {
        return arr[0].node->class;
    }
    return most_frequent_class;
}

void findKNearest(KDNode* root, double target[], int k) {
    Neighbor* neighbors = (Neighbor*)malloc(k * sizeof(Neighbor));
    int neighbor_count = 0;

    // Perform k-Nearest Neighbors search
    kNearestNeighbors(root, target, 0, neighbors, &neighbor_count, k);

    // Find the most frequent class from the neighbors
    target[dimensions] = findMostFrequentWithTieBreak(neighbors, k);

    // Free the neighbors array
    free(neighbors);
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
    dimensions = train_cols-1;    

    KDNode* root = NULL;
    int num_points;
    num_points = train_rows;

    for (int i = 0; i < num_points; i++) {
        double* point = &train_data[i * train_cols];
        root = insert(root, point, 0, train_data[i * train_cols + dimensions]);   
    }

    #pragma omp parallel for
    for(int i = 0; i < test_rows; i++){
        findKNearest(root, &test_data[i * test_cols], k);
    }

    writeResultsToFile(test_data, test_rows, test_cols, outfile);


    return 0;
}
