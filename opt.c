#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<stdbool.h>
#include <float.h>
#include <omp.h>

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

int dimensions;

typedef struct Node{
    double *point;
    struct Node *left;
    struct Node *right;
    double class;
}Node;

Node* newnode(double *point, double class){
    Node* new = (Node*)malloc(sizeof(Node));
    new->point = (double*)malloc(dimensions * sizeof(double));
    for (int i = 0; i < dimensions; i++){
        new->point[i] = point[i];
    }
    new->class = class;
    new->left = NULL;
    new->right = NULL;
    return new;
}

Node* insert(Node* root, double *point, double class, int depth){

    if(root==NULL) return newnode(point, class);
    else if(point[depth%dimensions] <= root->point[depth%dimensions]){
        root->left = insert(root->left, point, class, (depth+1));
    }
    else{
        root->right = insert(root->right, point, class, (depth+1));
    }
    return root;
}
//proof of concept function
Node* findmin(Node* root, int dim, int depth){
    if(root == NULL) return NULL;
    if(depth % dimensions == dim){
        if(root->left == NULL) return root;
        else return findmin(root->left, dim, depth + 1);
    }

    else{
        Node* left = findmin(root->left, dim, depth + 1);
        Node* right = findmin(root->right, dim, depth + 1);
        if (left == NULL) return right;
        if (right == NULL) return left;
        if (left->point[dim] <= right->point[dim]) return left;
        else return right;
    }
}

double distance(double *point1, double *point2, int size){
    double dist = 0.0;

    for(int i = 0; i<size; i++){
        dist += (point1[i] - point2[i])*(point1[i] - point2[i]);
    }
    return dist;
}

typedef struct NodeDist{
    Node* node;
    double dist;
}NodeDist;

void swap(NodeDist* a, NodeDist* b) {
    NodeDist temp = *a;
    *a = *b;
    *b = temp;
}

static inline void maxHeapify(NodeDist arr[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    while (left < n) {
        if (arr[left].dist > arr[largest].dist) {
            largest = left;
        }
        if (right < n && arr[right].dist > arr[largest].dist) {
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
void buildMaxHeap(NodeDist arr[], int k) {
    for (int i = (k / 2) - 1; i >= 0; i--) {
        maxHeapify(arr, k, i);
    }
}

double findMostFrequentWithTieBreak(NodeDist arr[], int k) {
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

void KNNSearch(Node* root, double* target, int k, int depth, NodeDist heap[]){
    if (root == NULL) return;

    double dist = distance(target, root->point, dimensions);

    if(root->left == NULL && root->right == NULL){
        if(dist<heap[0].dist){
            NodeDist newnode;
            newnode.dist = dist;
            newnode.node = root;
            heap[0] = newnode;
            maxHeapify(heap, k, 0);

        }
    }else {
        int axis = depth % dimensions;
        int dir = 0;
        if (target[axis] < root->point[axis]){
            KNNSearch(root->left, target, k, depth+1, heap);
        }else{
            KNNSearch(root->right, target, k, depth+1, heap);
            dir = 1;
        }
        if(dist<heap[0].dist){
            NodeDist newnode;
            newnode.dist = dist;
            newnode.node = root;
            heap[0] = newnode;
            maxHeapify(heap, k, 0);
        }
        //determine if branch is pruned
        if(distance(target, heap[0].node->point, dimensions)*9 > target[axis] - root->point[axis]){
            if(dir == 1){
                KNNSearch(root->left, target, k, depth+1, heap);
            }else{
                KNNSearch(root->right, target, k, depth+1, heap);
            }
        }
    }


}

void KNN(Node* root, double* target, int k){
    NodeDist heap[k];
    // init heap
    for(int i = 0; i<k; i++){
        NodeDist next;
        next.node = root;
        next.dist = DBL_MAX;
        heap[i] = next;
    }
    buildMaxHeap(heap, k);

    KNNSearch(root, target, k, 0, heap);


   target[dimensions] = findMostFrequentWithTieBreak(heap, k);

        
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
   
    Node* root = NULL;
    for (int i = 0; i < train_rows; i++) {
        double* point = &train_data[i * train_cols];
        root = insert(root, point, train_data[i * train_cols + dimensions], 0);   
    }
    #pragma omp parallel for
    for(int i = 0; i < test_rows; i++){
        KNN(root, &test_data[i * test_cols], k);
    }
    /*
    printf("finding min in dim 0\n");
    Node* min = findmin(root, 1, 0);
    for(int i = 0; i < dimensions; i++){
        printf("%lf,", min->point[i]);
    }
    */


    writeResultsToFile(test_data, test_rows, test_cols, outfile);


    return 0;
}