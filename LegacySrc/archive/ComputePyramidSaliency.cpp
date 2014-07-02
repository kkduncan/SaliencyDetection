/*
 * ComputePyramidSaliency.cpp
 *
 *  Created on: Apr 3, 2011
 *      Author: kesterduncan
 *
 *  This file computes the saliency of image pixels using
 *  an Iterative Kernel Density Estimation method.
 *
 *  Copyright (c) 2011 Kester Duncan
 *
 *  $Id: ComputePyramidSaliency$
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "image/RImage.h"
#include "image/SauvolasThreshold.h"
#include "gaussian/gaussian_filter.h"
#include "gaussian/GaussianPyramid.h"
#include "edge/canny.h"
#include "KDInfo.h"
#include "REM.h"

#define LOG2(X) (log(X) / 0.30103)

const int LEVELS = 4;

/** Global variables */
int 	width[LEVELS], height[LEVELS]; 	// Dimensions of pyramid images
float 	*gradDirections[LEVELS]; 		// Gradient directions map
float 	*magnitudes[LEVELS]; 			// Edge magnitudes map
KDInfo	*pixels[LEVELS];				// Stores relevant information for pixels
RImage 	maps[LEVELS]; 					// Saliency maps
int		*relevantPixels[LEVELS];		// Pixels that survive thresholding
RImage 	testImg; 						// Input image
GaussianPyramid pyr;					// Gaussian pyramid

/** Function Prototypes */
void 			magnitudeQuantization(int index);
void			computeIterativeSaliency(int index, int percentage);
void 			postProcessMap(float sigma);
static void 	allocateMemory();
static void 	deallocateMemory();
static char 	*getFileStem(const char *filename);
static void 	getFilename(char pathname[255], char output[255]);
static void 	convertValuesToChar(const int *srcVal, unsigned char **tgtVal, const int rows, const int cols);

void 			updateApplicableRegion(int index, Bounds bounds, KSumInfo ksum);
Bounds 			getApplicableBounds(Location locs[4], int neighborhood);
KSumInfo 		calculateKernelSum (int index, Location locs[4]);

void 			updateEntropy(KDInfo *p);
void 			updateSaliencyMap(int index);
int 			inValidImageBounds(int index, Location loc[4]);
void 			determineRelevance(int index);
int				getRelevantPixelCount(int index);

using namespace std;

int main(int argc, char *argv[]) {
	time_t	start, stop;
	char	imageFileName[100];
	char 	mapFileName[100];
	char	fileStem[100];
	int 	window; 			// neighborhood dimension
	int 	percentage; 		// sampling percentage
	int 	index;
	unsigned char *imgVals[LEVELS];

	if (argc != 3) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./ComputePyramidSaliency [IMAGE] [SAMPLING PERCENTAGE]\n");
		fprintf(stderr, "\te.g. ./ComputePyramidSaliency ads.pgm 25\n");
		exit(1);
	}

	/* Mark the start time for all processing */
	start = time(NULL);

	strcpy(imageFileName, argv[1]);
	strcpy(fileStem, getFileStem(imageFileName));
	printf("Evaluating saliency of image (using Gaussian pyramids): %s\n", imageFileName);

	percentage = atoi(argv[2]);

	testImg.read(imageFileName);
	maps[0].initialize(testImg.numRows, testImg.numCols);
	width[0] = testImg.numCols;
	height[0] = testImg.numRows;

	/* Perform Gaussian Pyramid Processing */
	pyr.setBaseLevel(testImg);
	pyr.reduceAll();
	for (int i = 1; i < LEVELS; i++) {
		width[i] = pyr.dimensions[i][1];
		height[i] = pyr.dimensions[i][0];
		maps[i].initialize(height[i], width[i]);

	}
	allocateMemory();

	/*
	 * Perform Relevance Processing
	 *
	 * -	we start with the highest pyramid level
	 * -	and compute the saliency of every pixel at this level
	 * -	we then perform Thresholding on the map returned and
	 * -	select the 25 neighbors of a pixel that survives thresholding for further processing
	 */
	for (index = LEVELS - 1; index >= 0; index--) {
		convertValuesToChar(pyr.levels[index].getChannel(RED), &imgVals[index], height[index], width[index]);
		canny(imgVals[index], height[index], width[index], 1.20, gradDirections[index], magnitudes[index]);
		magnitudeQuantization(index);
		initializeKDInfo(pixels[index], (width[index] * height[index]));

		computeIterativeSaliency(index, percentage);

		/*
		 * Perform Sauvola's Thresholding to determine the relevant pixels
		 */
		if (index > 0) {
//			SauvolasThreshold st; // Object to perform Sauvola's Thresholding
//			st.threshold(maps[index]);
			determineRelevance(index);
		}

		free(imgVals[index]);
	}

	sprintf(mapFileName, "results/pyramid_map_results/%s_map_it.jpg", fileStem);
	postProcessMap(16);

	maps[0].save(mapFileName);

	deallocateMemory();

	for (int i = 0; i < LEVELS; i++) {
		maps[i].flushMemory();
	}

//	free(gradDirections);
//	free(magnitudes);
//	free(pixels);


	/* Mark the end time for all processing and calculate */
	stop = time(NULL);
	printf("TOTAL PROCESSING TIME: %ld minute(s), %ld second(s)\n", (stop - start) / 60, (stop - start) % 60);

	return 0;
}


/**
 * Suppresses low magnitudes from the 'magnitudes' array
 * and increases magnitudes closer to the maximum magnitude
 * value
 */
void magnitudeQuantization(int index) {

	if (magnitudes[index] != NULL) {
		double maxMag = -999;

		for (int k = 0; k < (width[index] * height[index]); k++) {
			if (*(magnitudes[index] + k) > maxMag) {
				maxMag = magnitudes[index][k];
			}
		}

		/* Thresholds */
		double t1 = (double) maxMag / 3.0;
		double t2 = (double) maxMag / 1.5;	// was 2.0
		double t3 = (double) maxMag / 0.75; // was 1.5

		for (int k = 0; k < (width[index] * height[index]); k++) {

			if (magnitudes[index][k] >= 0 && magnitudes[index][k] < t1) {
				magnitudes[index][k] = 0.0005;

			} else if (magnitudes[index][k] >= t1 && magnitudes[index][k] < t2) {
				magnitudes[index][k] = 0.0025; // was 0.0025

			} else if (magnitudes[index][k] >= t2 && magnitudes[index][k] < t3) {
				magnitudes[index][k] = 0.0125; // was 0.0125

			} else if (magnitudes[index][k] >= t3){
				magnitudes[index][k] = 1.0;

			}
		}
	} else {
		fprintf(stderr, "ERROR: The magnitudes map is NULL.\n");
		exit(1);
	}
}


/**
 * Applies Gaussian filtering and image adjustment to a
 * saliency map
 */
void postProcessMap(float sigma) {
	int *filteredImage = NULL;
	float sig2 = 0.6 * sigma;

	/*
	 * First smoothing operation
	 */
	gaussianSmooth((unsigned int*) maps[0].getChannel(RED), maps[0].numRows, maps[0].numCols, sigma, &filteredImage);
	for (int i=0; i < height[0]; i++) {
		for (int j=0; j < width[0]; j++) {
			int ij = (i * width[0]) + j;
			maps[0](i, j, RED) = (int) (filteredImage[ij]);
			if (maps[0].isColorImage) {
				maps[0](i, j, GREEN) = (int) (filteredImage[ij]);
				maps[0](i, j, BLUE) = (int) (filteredImage[ij]);
			}
		}
	}

	free(filteredImage);
	filteredImage = NULL;

	/*
	 * Smoothing another time effectively gets rid of
	 * unsightly rings caused by the Gaussian
	 */
	gaussianSmooth((unsigned int*) maps[0].getChannel(RED), maps[0].numRows, maps[0].numCols, sig2, &filteredImage);
	for (int i=0; i < height[0]; i++) {
		for (int j=0; j < width[0]; j++) {
			int ij = (i * width[0]) + j;
			maps[0](i, j, RED) = (int) (filteredImage[ij]);
			if (maps[0].isColorImage) {
				maps[0](i, j, GREEN) = (int) (filteredImage[ij]);
				maps[0](i, j, BLUE) = (int) (filteredImage[ij]);
			}
		}
	}

	free(filteredImage);

	float maxVal = -999;
	float minVal = 999;

	for (int i = 0; i < height[0]; i++) {
		for (int j = 0; j < width[0]; j++) {
			if (maps[0](i, j, RED) > maxVal) {
				maxVal = maps[0](i, j);
			}
			if (maps[0](i, j, RED) < minVal) {
				minVal = maps[0](i, j);
			}
		}
	}

	/* Normalization */
	float sum = 0;
	for (int i = 0; i < height[0]; i++) {
		for (int j = 0; j < width[0]; j++) {
			maps[0](i, j, RED) = ((((maps[0](i, j, RED) - minVal) / (maxVal - minVal)) * 255.0));

			if (maps[0].isColorImage) {
				maps[0](i, j, GREEN) = ((((maps[0](i, j, GREEN) - minVal) / (maxVal - minVal)) * 255.0));
				maps[0](i, j, BLUE) = ((((maps[0](i, j, BLUE) - minVal) / (maxVal - minVal)) * 255.0));
			}
			sum += maps[0](i, j, RED);
		}
	}
}


/**
 * Process the saliency of each pixel in an image using an Iterative Approach
 * for calculating the Kernel Density Estimate of the local neighborhood
 * distribution, then calculating the entropy of that distribution for
 * determining the saliency.
 */
void computeIterativeSaliency(int index, int percentage) {
	const int M = 5;		// Dimension of the Local Pixel Neighborhood
	int pixelSamples;		// Maximum number of pixels sampled randomly
	int mainCount = 0;		// Counter for the main loop
	int mSquared = M * M;
	int halfM = M / 2;

	/* Ensure that the dimensions are valid before any processing */
	if (width[index] <= 0 || height[index] <= 0) {
		fprintf(stderr, "Error: Image dimensions are invalid.\n");
		exit(1);
	}

	/* Set the seed for randomization */
	srand(time(NULL));

	if (index == (LEVELS - 1)) {
		pixelSamples = (percentage / 100.0) * width[index] * height[index];

		while (mainCount < pixelSamples) {
			KSumInfo ksum;
			Bounds bounds;
			Location locs[4];

			/*
			 * Randomly select location of first pixel
			 */
			locs[0].row = rand() % height[index];
			locs[0].col = rand() % width[index];


			/*
			 * The other 3 pixels must be selected in the neighborhood
			 * of pixel 1
			 */
			locs[1].row = (rand() % M) + (locs[0].row - halfM);
			locs[1].col = (rand() % M) + (locs[0].col - halfM);

			locs[2].row = (rand() % M) + (locs[0].row - halfM);
			locs[2].col = (rand() % M) + (locs[0].col - halfM);

			locs[3].row = (rand() % M) + (locs[0].row - halfM);
			locs[3].col = (rand() % M) + (locs[0].col - halfM);

			ksum = calculateKernelSum(index, locs);
			bounds = getApplicableBounds(locs, M);
			updateApplicableRegion(index, bounds, ksum);
			mainCount++;
		}

	} else {
		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				KSumInfo ksum;
				Bounds bounds;
				Location locs[4];

				/*
				 * Randomly select location of first pixel
				 */
				locs[0].row = i;
				locs[0].col = j;

				int location = (locs[0].row * width[index]) + locs[0].col;

				if (relevantPixels[index][location] == 1) {
					/*
					 * The other 3 pixels must be selected in the neighborhood
					 * of pixel 1
					 */
					locs[1].row = (rand() % M) + (locs[0].row - halfM);
					locs[1].col = (rand() % M) + (locs[0].col - halfM);

					locs[2].row = (rand() % M) + (locs[0].row - halfM);
					locs[2].col = (rand() % M) + (locs[0].col - halfM);

					locs[3].row = (rand() % M) + (locs[0].row - halfM);
					locs[3].col = (rand() % M) + (locs[0].col - halfM);

					ksum = calculateKernelSum(index, locs);
					bounds = getApplicableBounds(locs, M);
					updateApplicableRegion(index, bounds, ksum);
					mainCount++;
				}
			}
		}
	}

	updateSaliencyMap(index);

}


/**
 * Calculate the intermediate kernel sum from the contribution of the
 * pixels in 'randomPixels'
 */
KSumInfo calculateKernelSum (int index, Location locs[4]) {
	KSumInfo ksum;
	double sampleDistance1, sampleAngle1, sampleMag1;
	double sampleDistance2, sampleAngle2, sampleMag2;
	double distanceKernel = 0, angleKernel = 0;
	double binDimension = 10.0; // Conceptual histogram's bin dimension
	double distanceBinWidth = sqrt(pow(width[index], 2) + pow(height[index], 2)) / binDimension;
	double angleBinWidth = 3.142 / binDimension;
	double twoPI = 6.283;
	double dNorm = (sqrt(twoPI) * distanceBinWidth);
	double aNorm = (sqrt(twoPI) * angleBinWidth);

	if (inValidImageBounds(index, locs)) {
		int loc1 = (locs[0].row * width[index]) + locs[0].col;
		int loc2 = (locs[1].row * width[index]) + locs[1].col;
		int loc3 = (locs[2].row * width[index]) + locs[2].col;
		int loc4 = (locs[3].row * width[index]) + locs[3].col;

		sampleDistance1 = sqrt(pow((locs[0].row - locs[1].row), 2) + pow((locs[0].col - locs[1].col), 2));
		sampleDistance2 = sqrt(pow((locs[2].row - locs[3].row), 2) + pow((locs[2].col - locs[3].col), 2));

		sampleAngle1 = fmod((fabs(gradDirections[index][loc1] - gradDirections[index][loc2])), 6.283);
		sampleAngle2 = fmod((fabs(gradDirections[index][loc3] - gradDirections[index][loc4])), 6.283);

		sampleMag1 = fabs(magnitudes[index][loc1] - magnitudes[index][loc2]);
		sampleMag2 = fabs(magnitudes[index][loc3] - magnitudes[index][loc4]);

		ksum.magSum1 = sampleMag1;
		ksum.magSum2 = sampleMag2;

		distanceKernel = (1/dNorm) * exp((pow(sampleDistance1 - sampleDistance2, 2) / (-2 * pow(distanceBinWidth, 2))));
		angleKernel = (1/aNorm) * exp((pow(sampleAngle1 - sampleAngle2, 2) / (-2 * pow(angleBinWidth, 2))));

		if (sampleMag1 > 0 && sampleMag2 > 0 && distanceKernel > 0 && angleKernel > 0) {
			ksum.kernelSum = (sampleMag1 * sampleMag2 * distanceKernel * angleKernel);
		} else {
			ksum.kernelSum = 0;
		}
	} else {
		ksum.kernelSum = 0;
		ksum.magSum1 = 0;
		ksum.magSum2 = 0;
	}

	return ksum;

}


/**
 * Determine if pixel locations are within the valid image
 * dimensions
 */
int inValidImageBounds(int index, Location loc[4]) {
	int valid = 0;
	int validCount = 0;

	int i = 0;
	while (i < 4) {
		if (loc[i].row >= 0 && loc[i].row < height[index] && loc[i].col >= 0 && loc[i].col < width[index]) {
			validCount++;
		}
		i++;
	}

	if (validCount == 4) {
		valid = 1;
	} else {
		valid = 0;
	}

	return (valid);
}


/**
 * Gets the applicable region of the pixels forming two samples
 * that can be updated
 */
Bounds getApplicableBounds(Location locs[4], int neighborhood) {
	Bounds bounds;
	int L = neighborhood;		// Local neighborhood dimension
	int maxX = locs[0].col;
	int minX = locs[0].col;
	int maxY = locs[0].row;
	int minY = locs[0].row;
	int xDiff;		// Differences between the max / min x values
	int yDiff;		// Differences between the max / min y values
	int xDisp;		// The disparity in the x axis for forming an M x M neighborhood
	int yDisp;		// The disparity in the y axis for forming an M x M neighborhood

	/* Get the maximum / minimum x and y values of the pixels */
	for (int i = 1; i < 4; i++) {
		if (locs[i].col > maxX) {
			maxX = locs[i].col;
		}

		if (locs[i].row > maxY) {
			maxY = locs[i].row;
		}

		if (locs[i].col < minX) {
			minX = locs[i].col;
		}

		if (locs[i].row < minY) {
			minY = locs[i].row;
		}
	}

	/* Calculate the differences between the max / min values */
	xDiff = maxX - minX;
	yDiff = maxY - minY;

	/* Get the x and y disparity */
	xDisp = (L - xDiff) - 1;
	yDisp = (L - yDiff) - 1;

	/* Calculate the applicable bounds */
	bounds.topLeft.col = minX - xDisp;
	bounds.topLeft.row = minY - yDisp;

	bounds.topRight.col = maxX + xDisp;
	bounds.topRight.row = minY - yDisp;

	bounds.botLeft.col = minX - xDisp;
	bounds.botLeft.row = maxY + yDisp;

	bounds.botRight.col = maxX + xDisp;
	bounds.botRight.row = maxY + yDisp;

	return bounds;
}


/**
 * Used for updating the kernel sums within the applicable region
 * specified
 */
void updateApplicableRegion(int index, Bounds bounds, KSumInfo ksum) {
	int sampleCountLimit = 128; // change back to 30 - 100

	for (int i = bounds.topLeft.row; i <= bounds.botLeft.row; i++) {
		for (int j = bounds.topLeft.col; j <= bounds.topRight.col; j++) {
			/*
			 * Must be within the image dimensions
			 */
			if (i >= 0 && i < height[index] && j >= 0 && j < width[index]) {
				int ij = (i * width[index]) + j;

				if (pixels[index][ij].sampleCount < sampleCountLimit) {
					pixels[index][ij].kernelSum += ksum.kernelSum;
					pixels[index][ij].magSum1 += ksum.magSum1;
					pixels[index][ij].magSum2 += ksum.magSum2;
					pixels[index][ij].sampleCount++;

					/*
					 * Update the pixel entropy every N iterations
					 */
					// was 32
					if (((pixels[index][ij].sampleCount + 1) % 16) == 0) {
						updateEntropy(&pixels[index][ij]);
					}
				}
			}
		}
	}
}


/**
 * This is used to update the entropy of a pixel in the image
 * after N updates. p is a ptr to a KDInfo object with the
 * required information for calculating the entropy.
 */
void updateEntropy(KDInfo *p) {
	if (p->sampleCount > 0) {
		double totalMagSum = p->magSum1 * p->magSum2;
		double estimation = 0;

		/* Special Case Handling */
		if (totalMagSum <= 0) {
			totalMagSum = p->sampleCount;
		}

		if (p->kernelSum < 0 || isnan(p->kernelSum)) {
			p->kernelSum = 0;
		}

		estimation = p->kernelSum / totalMagSum;

		if (estimation <= 0 || isnan(estimation)) {
			p->entropy = -1;
		} else {
			p->entropy = -1 * LOG2(estimation * estimation);
		}
	}
}


/**
 * Updates the saliency map to reflect the current calculations
 * This is called by computeIterativeSaliency after every
 * N iterations.
 */
void updateSaliencyMap(int index) {
	if (pixels != NULL) {
		/* For proper visualization */
		float maxEntropy = -999.0;
		float minEntropy = 999.0;

		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				int ij = (i * width[index]) + j;
				if (pixels[index][ij].entropy > maxEntropy) {
					maxEntropy = pixels[index][ij].entropy;
				}
			}
		}

		/*
		 * Areas with 0 magnitude difference are assigned an entropy
		 * value of -1. These areas are now given the value of the
		 * maximum entropy estimated.
		 *
		 * We find the minimum entropy in tandem
		 */
		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				int ij = (i * width[index]) + j;
				if (pixels[index][ij].entropy == -1) {
					pixels[index][ij].entropy = maxEntropy;
				}

				if (pixels[index][ij].entropy < minEntropy) {
					minEntropy = pixels[index][ij].entropy;
				}
			}
		}

		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				int ij = (i * width[index]) + j;
				pixels[index][ij].entropy -= minEntropy;
			}
		}

		/* Also adjust the maximum entropy */
		 maxEntropy -= minEntropy;

		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				int ij = (i * width[index]) + j;

				if (maxEntropy > 0) {
					maps[index](i, j) = (int) (255.0 - ((pixels[index][ij].entropy / maxEntropy) * 255.0));
				} else {
					maps[index](i, j) = 0;
				}
			}
		}

	} else {
		printf("Error: Pixels array is NULL.\n");
		exit(1);
	}
}


/**
 * Get the relevant pixels for a pyramid level specified by index
 * The relevant pixels in the successive level are then marked
 */
void determineRelevance(int index) {
	int window = 1;
	if (index > 0) {
		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				int ij = (i * width[index]) + j;

				if (maps[index](i, j) >= 217) {
					relevantPixels[index][ij] = 1;
				} else {
					relevantPixels[index][ij] = 0;
				}
			}
		}

		for (int i = 0; i < height[index - 1]; i++) {
			for (int j = 0; j < width[index - 1]; j++) {
				int ij = (i * width[index - 1]) + j;

				int row = i / 2;
				int col = j / 2;

				int mn = (row * width[index]) + col;

				relevantPixels[index - 1][ij] = relevantPixels[index][mn];

//				if (relevantPixels[index - 1][ij] == 1) {
//					for (int m = -window; m <= window; m++) {
//						for (int n = -window; n <= window; n++) {
//							int y = (i + m);
//							int x = (j + n);
//
//							int yx = (y * width[index - 1]) + x;
//
//							if (y < 0) y = 0;
//							if (x < 0) x = 0;
//							if (y >= height[index - 1]) y = height[index - 1] - 1;
//							if (x >= width[index - 1]) x = width[index - 1] - 1;
//
//							if (relevantPixels[index - 1][yx] == 0) {
//								relevantPixels[index - 1][yx] = 2;
//							}
//						}
//					}
//				}
			}
		}

	} else {
		fprintf(stderr, "ERROR: Maximum relevance level reached.\n");
		exit(1);
	}
}



int	getRelevantPixelCount(int index) {
	int count = 0;

	for (int i = 0; i < (height[index] * width[index]); i++) {
		if (relevantPixels[index][i] > 0) {
			count++;
		}
	}

	return count;
}

/**
 * Allocates memory for the global arrays
 * i.e. gradDirections, magnitudes, and entropies
 */
static void allocateMemory() {
	for (int i=0; i < LEVELS; i++) {
		pixels[i] = (KDInfo *) calloc(width[i] * height[i], sizeof(KDInfo));
		if (pixels[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the pixels map.\n");
			exit(1);
		}

		gradDirections[i] = (float *) calloc(width[i] * height[i], sizeof(float));
		if (gradDirections[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the gradDirections map.\n");
			exit(1);
		}

		magnitudes[i] = (float *) calloc(width[i] * height[i], sizeof(float));
		if (magnitudes[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the magnitudes map.\n");
			exit(1);
		}

		relevantPixels[i] = (int*) calloc(width[i] * height[i], sizeof(int));
		if (relevantPixels[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the relevance array.\n");
			exit(1);
		}
	}
}

static void deallocateMemory() {
	for (int i=0; i < LEVELS; i++) {
		free(pixels[i]);
		free(gradDirections[i]);
		free(magnitudes[i]);
		free(relevantPixels[i]);
	}

}


/**
 * Utility function to get a file name minus the extension
 */
static char *getFileStem(const char *filename) {
	char *str = new char [(strlen(filename) - 3)];
	char *chPtr = strrchr(filename, '.');
	int index = chPtr - filename;
	int i = 0;;

	while (i < index) {
		str[i] = filename[i];
		i++;
	}
	str[i] = '\0';

	char *output = new char [strlen(str)];
	getFilename(str, output);

	delete(str);

	return output;
}


static void getFilename(char pathname[255], char output[255]) {
  int i = 0, j = 0;

  while (pathname[i] != '\0') {
    if (pathname[i] == '/') j = 0;
    else output[j++] = pathname[i];
    i++;
  }
  output[j++] =  '\0';
}


/**
 * Utility function to convert an integer array into a character
 * array
 */
static void convertValuesToChar(const int *srcVal, unsigned char **tgtVal, const int rows, const int cols) {
	(*tgtVal) = (unsigned char*) calloc(rows * cols, sizeof(unsigned char));

	if (srcVal == NULL) {
		fprintf(stderr, "Empty source val at line: %d in %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if ((*tgtVal) == NULL) {
		fprintf(stderr, "Couldn't allocate memore for tgtVal in %s, line: %d\n", __FILE__, __LINE__);
		exit(1);
	}

	for (int i=0; i < rows; i++) {
		for (int j=0; j < cols; j++) {
			int ij = (i * cols) + j;
			(*tgtVal)[ij] = (unsigned char) srcVal[ij];
		}
	}
}
