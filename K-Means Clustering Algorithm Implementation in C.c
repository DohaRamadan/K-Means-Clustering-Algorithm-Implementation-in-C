#include <omp.h>
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <math.h>

struct Point {
	int x, y, newCluster, oldCluster;
	double minDist;
	double* distances;
};

struct Centroid {
	double x, y;
};

int maxim(int num1, int num2);
int minim(int num1, int num2);
double calcDistance(struct Point p1, struct Centroid p2);
bool AreSame(double a, double b);

int main()
{
	FILE* f;
	int numOfPoints = 0;
	// opening the file to know how many points are there.
	const char* filename = "data.txt";
	f = fopen(filename, "r");
	if (f == NULL) {
		perror("Error opening file");
		return(-1);
	}
	else {
		char c;
		for (c = getc(f); c != EOF; c = getc(f)) {
			if (c == '\n') // Increment count if this character is newline
				numOfPoints = numOfPoints + 1;
		}
		// Close the file
		fclose(f);
		numOfPoints++;
	}

	struct Point* data = (struct Point*)malloc(numOfPoints * sizeof(struct Point)); // array for storing the points

	// re-open file for storing the data points in the data array 
	int MAX_X = INT_MIN;
	int MIN_X = INT_MAX;
	int MAX_Y = INT_MIN;
	int MIN_Y = INT_MAX;
	f = fopen(filename, "r");
	if (f == NULL) {
		perror("Error opening file"); // print error
		return(-1);
	}
	else {
		int indx = 0;
		while (1) {
			int x, y;
			if (feof(f)) {
				break;
			}
			fscanf(f, "%d %d", &x, &y);
			struct Point p;
			p.x = x;
			p.y = y;
			p.minDist = INT_MAX;
			data[indx] = p;
			MAX_X = maxim(MAX_X, x);
			MIN_X = minim(MIN_X, x);
			MAX_Y = maxim(MAX_Y, y);
			MIN_Y = minim(MIN_X, y);
			indx++;
		}
		fclose(f); // close file
	}
	int i, j;
	int np;    // number of threads
	int nc;         // number of clusters
	#pragma omp parallel
	{
		np = omp_get_num_threads(); 
	}
	nc = np;
	for (i = 0; i < numOfPoints; i++) {
		data[i].distances = (double*)malloc(np * sizeof(double));
	}

	
	// Genertaing random intial centroids.
	struct Centroid* centroids = (struct Centroid*)malloc(nc * sizeof(struct Centroid));
	printf("Random Centroids: \n"); 
	for (i = 0; i < nc; i++) {
		struct Centroid c;
		c.x = (rand() % (MAX_X - MIN_X + 1)) + MIN_X;
		c.y = (rand() % (MAX_Y - MIN_Y + 1)) + MIN_Y;
		printf("Centroid %d (%lf, %lf)\n", i+1, c.x, c.y);
		centroids[i] = c;
	}
	
	int* clusterNumOfPoints = (int*)malloc(nc * sizeof(int));   // number of points in each cluster
	struct Point** clusters;
	clusters = (struct Point**)malloc(nc * sizeof(struct Point*)); // array of array of points represeting each cluster
	int count = 0; // counter for number of iterations
	while (true) {
		bool flag = true; //flag to indicate whether isClusterChanged is true for at least one point
		bool flag2 = true; // flag to indicate whether isCentroidChanged is true for at least one centroid
		// initializing isClusterChanged and isCentroidChanged with false
		bool* isClusterChanged = (bool*)malloc(numOfPoints * sizeof(bool));
		bool* isCentroidChanged = (bool*)malloc(nc * sizeof(bool));
		for (i = 0; i < numOfPoints; i++) {
			isClusterChanged[i] = false;
		}
		for (i = 0; i < nc; i++) {
			isCentroidChanged[i] = false;
		}

		// Calculating the distance between each point and each centroid
		#pragma omp parallel for private(i,j) num_threads(np) 
		for (i = 0; i < numOfPoints; i++) {
			for (j = 0; j < nc; j++) {
				double ds = calcDistance(data[i], centroids[j]);
				data[i].distances[j] = ds;
			}
		}

		for (i = 0; i < nc; i++) {
			clusterNumOfPoints[i] = 0;
		}

		// stroing the old cluster in the Point.oldCluster attribute
		// -only after the first iterarion because in the 1st iteration there is no oldCluster- . 
		if (count > 0) {
			for (i = 0; i < numOfPoints; i++) {
				data[i].oldCluster = data[i].newCluster; 
			}
		}
		// Assigning points to culstures based on minmum distance
		for (i = 0; i < numOfPoints; i++) {
			double min = 1.79769e+308;
			for (j = 0; j < nc; j++) {
				if (data[i].distances[j] < min) {
					min = data[i].distances[j];
					data[i].minDist = min;
					data[i].newCluster = j;
				}
			}
		}
		// Checking whether each point cluster has changed or not. 
		if (count > 0) {
			for (i = 0; i < numOfPoints; i++) {
				if (data[i].oldCluster != data[i].newCluster) {
					isClusterChanged[i] = true;
				}
				else {
					isClusterChanged[i] = false;
				}
			}
		}

		// Calculating the number of points in each cluster
		for (i = 0; i < numOfPoints; i++) {
			clusterNumOfPoints[data[i].newCluster]++;
		}

		// dynamically allocating meomry for each cluster based on the number of points that belongs to that cluster.
		for (i = 0; i < nc; i++) {
			clusters[i] = (struct Point*)malloc((clusterNumOfPoints[i]) * sizeof(struct Point));
		}
		// Storing the points in its cluster array 
		for (i = 0; i < nc; i++) {
			int count = 0;
			for (j = 0; j < numOfPoints; j++) {
				if (data[j].newCluster == i) {
					clusters[i][count] = data[j];
					count++;
				}
			}
		}

		//Calculate the mean for each cluster as new cluster centroid
		#pragma omp parallel for private(i,j) num_threads(np) 
		for (i = 0; i < nc; i++) {
			double sum_x = 0.0, sum_y = 0.0;
			for (j = 0; j < clusterNumOfPoints[i]; j++) {
				sum_x += clusters[i][j].x;
				sum_y += clusters[i][j].y;
			}
			struct Centroid oldCentroid = centroids[i];
			struct Centroid newCentroid;
			if (clusterNumOfPoints[i] > 0) {
				newCentroid.x = (double)(sum_x / clusterNumOfPoints[i]);
				newCentroid.y = (double)(sum_y / clusterNumOfPoints[i]);
				if (AreSame(newCentroid.x, oldCentroid.x) && AreSame(newCentroid.y, newCentroid.y)) {
					isCentroidChanged[i] = false;
				}
				else {
					isCentroidChanged[i] = true;
				}
				centroids[i] = newCentroid;
			}
		}

		// Checking whether the cluster has changed for at least one point 
		// and if the centroid is changed for at least one centroid
		if (count > 0) {
			for (i = 0; i < numOfPoints; i++) {
				if (isClusterChanged[i] == true) {
					flag = false;
					break;
				}
			}
			for (i = 0; i < nc; i++) {
				if (isCentroidChanged[i] == true) {
					flag2 = false;
					break;
				}
			}

		}
		// increamting the iteration counter
		count++;

		// if all points' cluster are the same and all centroids are the same -no change- and we are not in the first iteration, then stop. Otherwise, continue.
		if (flag && count > 1 && flag2) {
			break;
		}
		else {
			continue;
		}
	}

	// the final result 
	for (i = 0; i < nc; i++) {
		printf("Cluster %d (%lf, %lf): \n", i + 1, centroids[i].x, centroids[i].y);
		for (j = 0; j < clusterNumOfPoints[i]; j++)
		{
			printf("(%d, %d)\n", clusters[i][j].x, clusters[i][j].y);
		}
	}

}
int maxim(int num1, int num2)
{
	return (num1 > num2) ? num1 : num2;
}

int minim(int num1, int num2)
{
	return (num1 > num2) ? num2 : num1;
}
double calcDistance(struct Point p1, struct Centroid p2) {
	return (double)sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}
bool AreSame(double a, double b)
{
	return fabs(a - b) < 0.000001;
}