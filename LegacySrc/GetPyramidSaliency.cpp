/*
 * GetPyramidSaliency.cpp
 *
 *  Created on: Aug 21, 2011
 *
 *  This file computes the saliency of image pixels using
 *  an Iterative Kernel Density Estimation method alongside
 *  kernel descriptors (gradient, local binary pattern, color).
 *
 *  Copyright (c) 2011 Kester Duncan
 *
 *  $Id: GetPyramidSaliency$
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "image/RImage.h"
#include "image/LBP.h"
#include "gaussian/gaussian_filter.h"
#include "gaussian/GaussianPyramid.h"
#include "edge/canny.h"
#include "KDInfo.h"

#define LOG2(X) (log(X) / 0.30103)

/** Global variables */
const int LEVELS = 4;					// Number of Pyramid Levels
int 	width[LEVELS], height[LEVELS]; 	// Dimensions of image
float	*intensities[LEVELS];			// Image intensity
float 	*gradDirections[LEVELS]; 		// Gradient directions map
float 	*magnitudes[LEVELS]; 			// Edge magnitudes map
int		*lbp[LEVELS];					// LBP map
float	*stdDev[LEVELS];				// Standard deviation of pixel values in a 3x3 neighborhood around a pixel
RImage 	src;		 					// Input image
RImage 	maps[LEVELS]; 					// Saliency map
KDInfo	*pixels[LEVELS];				// Stores relevant information for pixels
int 	*relevantPixels[LEVELS];		// Pixels that survice thresholding
const int M = 7;						// Dimension of the Local Pixel Neighborhood
GaussianPyramid pyr;					// Gaussian pyramid

/** Function Prototypes */
void 			magnitudeQuantization(int index);
void 			getStandardDeviations(int index);
void			computeKernelSaliency(int index, int percentage);
void 			postProcessMap(float sigma);
static void 	allocateMemory();
static void 	deallocateMemory();
static char 	*getFileStem(char *filename);
static void 	convertValuesToChar(const float *srcVal, unsigned char **tgtVal, int rows, int cols);
void			getIntensities(int index);
void 			updateApplicableRegion(int index, Bounds bounds, KSumInfo gsum, KSumInfo lsum, KSumInfo csum);
Bounds 			getApplicableBounds(Location locs[2], int neighborhood);
KSumInfo 		calculateGradKernelSum (int index, Location locs[2]);
KSumInfo 		calculateLBPKernelSum (int index, Location locs[2]);
KSumInfo 		calculateColorKernelSum (int index, Location locs[2]);
void 			updateEntropy(KDInfo *p);
void 			updateSaliencyMap(int index);
int 			inValidImageBounds(int index, Location loc[2]);
void 			determineRelevance(int index);
int				getRelevantPixelCount(int index);
void 			writeArrayToFile(const float* array, char *textFilename, int index);

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
	LBP localBinPat(32, 4);

	if (argc != 3) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./GetPyramidSaliency [IMAGE] [SAMPLING PERCENTAGE]\n");
		fprintf(stderr, "\te.g. ./GetPyramidSaliency ads.pgm 25\n");
		exit(-1);
	}

	/* Mark the start time for all processing */
	start = time(NULL);

	strcpy(imageFileName, argv[1]);
	strcpy(fileStem, getFileStem(imageFileName));
	printf("Evaluating (pyramid) saliency of image: %s\n", imageFileName);

	percentage = atoi(argv[2]);

	src.read(imageFileName);
	if (!src.isColorImage) {
		src.convertToColor();
	}

	maps[0].initialize(src.numRows, src.numCols);
	width[0] = src.numCols;
	height[0] = src.numRows;

	/* Perform Gaussian Pyramid Processing */
	pyr.setBaseLevel(src);
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
		getIntensities(index);
		convertValuesToChar(intensities[index], &imgVals[index], height[index], width[index]);
		canny(imgVals[index], height[index], width[index], 1.20, gradDirections[index], magnitudes[index]);
		localBinPat.getOriginalLBP(pyr.levels[index], lbp[index]);
		magnitudeQuantization(index);
		getStandardDeviations(index);
		initializeKDInfo(pixels[index], (width[index] * height[index]));

		computeKernelSaliency(index, percentage);

		/*
		 * Perform Sauvola's Thresholding to determine the relevant pixels
		 */
		if (index > 0) {
			determineRelevance(index);
		}

		free(imgVals[index]);

//		sprintf(mapFileName, "%s_inter_%d.jpeg", fileStem, index);
//		maps[index].save(mapFileName);

	}


	sprintf(mapFileName, "%s_map_it.jpeg", fileStem);
	postProcessMap(10);

	/*
	 * Make composite map
	 */
	RImage comp;
	comp.isColorImage = 1;
	comp.initialize(maps[0].numRows, (maps[0].numCols * 2));

	for (int m=0; m < maps[0].numRows; m++) {
		for (int n=0; n < maps[0].numCols; n++) {
			if (src.isColorImage) {
				comp(m, n, RED) = src(m, n, RED);
				comp(m, n, GREEN) = src(m, n, GREEN);
				comp(m, n, BLUE) = src(m, n, BLUE);
			} else {
				comp(m, n, RED) = src(m, n, RED);
				comp(m, n, GREEN) = src(m, n, RED);
				comp(m, n, BLUE) = src(m, n, RED);
			}
		}
	}

#define SHOW_MAP

#ifdef SHOW_MAP
	maps[0].makeThermalColorImage();
	for (int m=0; m < maps[0].numRows; m++) {
		for (int n=0; n < maps[0].numCols; n++) {
			int col = n + maps[0].numCols;
			comp(m, col, RED) = maps[0](m, n, RED);
			comp(m, col, GREEN) = maps[0](m, n, GREEN);
			comp(m, col, BLUE) = maps[0](m, n, BLUE);

		}
	}

#else
	for (int m=0; m < maps[0].numRows; m++) {
		for (int n=0; n < maps[0].numCols; n++) {
			int col = n + maps[0].numCols;
			double value = maps[0](m, n, RED);
			if (value > 128) {
				comp(m, col, RED) = src(m, n, RED);
				comp(m, col, GREEN) = src(m, n, RED);
				comp(m, col, BLUE) = src(m, n, RED);
			}
		}
	}
#endif

	deallocateMemory();

	for (int i = 0; i < LEVELS; i++) {
		maps[i].flushMemory();
	}

	comp.save(mapFileName);
	comp.flushMemory();
	src.flushMemory();

	/* Mark the end time for all processing and calculate */
	stop = time(NULL);
	printf("TOTAL PROCESSING TIME: %ld minute(s), %ld second(s)\n",
			(stop - start) / 60, (stop - start) % 60);

	return 0;
}


/**
 * Suppresses low magnitudes from the 'magnitudes' array
 * and increases magnitudes closer to the maximum magnitude
 * value
 */
void magnitudeQuantization(int index) {
	if (magnitudes[index] != NULL) {
		float magSum = 0;

		for (int k = 0; k < (width[index] * height[index]); k++) {
			magSum += pow((double) magnitudes[index][k], 2.0) + 0.05;
		}

		double denom = sqrt((double) magSum);
		for (int k = 0; k < (width[index] * height[index]); k++) {
			magnitudes[index][k] = magnitudes[index][k] / denom;
		}

	} else {
		fprintf(stderr, "ERROR: The magnitudes map is NULL.\n");
		exit(1);
	}
}


/**
 * Gets the standard deviation of pixel intensities in a 3x3
 * neighborhood around a central pixel.
 */
void getStandardDeviations(int index) {
	if (stdDev[index] != NULL) {
		double sdSum = 0;

		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				float sum = 0;
				float avg = 0;
				float sd = 0;
				int count = 0;
				int ij = (i * width[index]) + j;
				int neighbor = 1;

				for (int m = -neighbor; m <= neighbor; m++) {
					for (int n = -neighbor; n <= neighbor; n++) {
						int row = i + m;
						int col = j + n;

						if (row >= 0 && col >= 0 && row < height[index] && col < width[index]) {
							count++;
							int ij = (row * width[index]) + col;
							sum += intensities[index][ij];
						}
					}
				}

				avg = sum / count;
				sum = 0;
				for (int m = -neighbor; m <= neighbor; m++) {
					for (int n = -neighbor; n <= neighbor; n++) {
						int row = i + m;
						int col = j + n;

						if (row >= 0 && col >= 0 && row < height[index] && col < width[index]) {
							count++;
							int ij = (row * width[index]) + col;
							sum += pow((double)(intensities[index][ij] - avg), 2.0);
						}
					}
				}
				sd = sum / (count - 1) ;
				stdDev[index][ij] = sd;
				sdSum += pow((double) sd, 2.0) + 0.05;;
			}
		}

		double denom = sqrt((double) sdSum);
		for (int k = 0; k < (width[index] * height[index]); k++) {
			stdDev[index][k] = stdDev[index][k] / denom;
		}

	} else {
		fprintf(stderr, "ERROR: The lbp map is NULL.\n");
		exit(1);
	}
}


/**
 * Simply place the image intensities into an array
 */
void getIntensities(int index) {
	if (intensities[index] != NULL) {
		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				int ij = (i * width[index]) + j;
				intensities[index][ij] = (pyr.levels[index](i, j, RED) +
						pyr.levels[index](i, j, GREEN) +
						pyr.levels[index](i, j, BLUE)) / 3.0;
			}
		}

	} else {
		fprintf(stderr, "ERROR: The intensity map is NULL.\n");
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
 * Process the saliency of each pixel in an image using an Iterative Kernel
 * Approach for calculating the Kernel Density Estimate of the local neighborhood
 * distribution, then calculating the entropy of that distribution for
 * determining the saliency.
 */
void computeKernelSaliency(int index, int percentage) {
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
		pixelSamples = width[index] * height[index];

		while (mainCount < pixelSamples) {
			KSumInfo gsum, lsum, csum;
			Bounds bounds;
			Location locs[2];			// Array to store randomly selected pixel locations

			/*
			 * Randomly select location of first pixel
			 */
			locs[0].row = rand() % height[index];
			locs[0].col = rand() % width[index];

			/*
			 * The other pixel must be selected in the neighborhood
			 * of pixel 1
			 */
			locs[1].row = (rand() % M) + (locs[0].row - halfM);
			locs[1].col = (rand() % M) + (locs[0].col - halfM);

			gsum = calculateGradKernelSum(index, locs);
			lsum = calculateLBPKernelSum(index, locs);
			csum = calculateColorKernelSum(index, locs);
			bounds = getApplicableBounds(locs, M);
			updateApplicableRegion(index, bounds, gsum, lsum, csum);

			mainCount++;

		}
	} else {
		for (int i = 0; i < height[index]; i++) {
			for (int j = 0; j < width[index]; j++) {
				KSumInfo gsum, lsum, csum;
				Bounds bounds;
				Location locs[2];			// Array to store randomly selected pixel locations

				/*
				 * Randomly select location of first pixel
				 */
				locs[0].row = i;
				locs[0].col = j;

				int location = (locs[0].row * width[index]) + locs[0].col;

				if (relevantPixels[index][location] == 1) {

					locs[1].row = (rand() % M) + (locs[0].row - halfM);
					locs[1].col = (rand() % M) + (locs[0].col - halfM);

					if (locs[1].row >= 0 && locs[1].row < height[index] &&
						locs[1].col >= 0 && locs[1].col < width[index]) {
						gsum = calculateGradKernelSum(index, locs);
						lsum = calculateLBPKernelSum(index, locs);
						csum = calculateColorKernelSum(index, locs);
						bounds = getApplicableBounds(locs, M);
						updateApplicableRegion(index, bounds, gsum, lsum, csum);

					}
				}
			}
		}
	}

	updateSaliencyMap(index);

}


/**
 * Calculate the gradient kernel sum from the contribution of the
 * pixels in 'randomPixels'
 */
KSumInfo calculateGradKernelSum (int index, Location locs[2]) {
	KSumInfo ksum;
	double sampleDist1, sampleAngle1, sampleMag1, sampleMag2;
	double distKernel = 0, angleKernel = 0;
	double binDimension = 10.0;
	double distBinWidth = sqrt(pow((double) M, 2.0) + pow((double) M, 2.0)) / binDimension;
	double angleBinWidth = 3.142 / binDimension;
	double twoPI = 6.283;

	double pNorm = (2.5066 * distBinWidth);
	double gNorm = (2.5066 * angleBinWidth);

	if (inValidImageBounds(index, locs)) {
		int loc1 = (locs[0].row * width[index]) + locs[0].col;
		int loc2 = (locs[1].row * width[index]) + locs[1].col;

		sampleDist1 = sqrt(pow((double)(locs[0].row - locs[1].row), 2.0) + pow((locs[0].col - locs[1].col), 2.0)) / distBinWidth;

		sampleAngle1 = sqrt(pow(sin((double)gradDirections[index][loc1]) - sin((double)gradDirections[index][loc2]), 2.0) +
						   pow(cos((double)gradDirections[index][loc1]) - cos((double)gradDirections[index][loc2]), 2.0));

		sampleMag1 = magnitudes[index][loc1];
		sampleMag2 = magnitudes[index][loc2];

		distKernel = exp(sampleDist1 * -distBinWidth);
		angleKernel = exp(sampleAngle1 * -angleBinWidth);

		ksum.kernelSum = sampleMag1 * sampleMag2 * angleKernel * distKernel;

	} else {
		ksum.kernelSum = 0;
		ksum.magSum1 = 0;
		ksum.magSum2 = 0;
	}

	return ksum;
}


/**
 * Calculate the local binary pattern kernel sum
 */
KSumInfo calculateLBPKernelSum (int index, Location locs[2]) {
	KSumInfo ksum;
	double sampleDist1, sampleLBP1, sampleSD1, sampleSD2;
	double distKernel = 0, lbpKernel = 0;
	double binDimension = 10.0;
	double distBinWidth = sqrt(pow((double) M, 2.0) + pow((double) M, 2.0)) / binDimension;
	double lbpBinWidth = 255 / binDimension;
	double twoPI = 6.283;

	if (inValidImageBounds(index, locs)) {
		int loc1 = (locs[0].row * width[index]) + locs[0].col;
		int loc2 = (locs[1].row * width[index]) + locs[1].col;

		sampleDist1 = sqrt(pow((double)(locs[0].row - locs[1].row), 2.0) + pow((locs[0].col - locs[1].col), 2.0)) / distBinWidth;

		sampleLBP1 = sqrt(pow((double)(lbp[index][loc1] - lbp[index][loc2]), 2.0)) / 255;
		sampleSD1 = stdDev[index][loc1];
		sampleSD2 = stdDev[index][loc2];

		distKernel = exp(sampleDist1 * -distBinWidth);
		lbpKernel = exp(sampleLBP1 * -lbpBinWidth);

		ksum.kernelSum = sampleSD1 * sampleSD2 * lbpKernel * distKernel;

	} else {
		ksum.kernelSum = 0;
		ksum.magSum1 = 0;
		ksum.magSum2 = 0;
	}

	return ksum;
}


/**
 * Calculate the color kernel sum
 */
KSumInfo calculateColorKernelSum (int index, Location locs[2]) {
	KSumInfo ksum;
	double sampleDist1, sampleColor1;
	double distKernel = 0, colorKernel = 0;
	double binDimension = 10.0;
	double distBinWidth = sqrt(pow((double) M, 2.0) + pow((double) M, 2.0)) / binDimension;
	double colorBinWidth = 255 / binDimension;
	double twoPI = 6.283;

	if (inValidImageBounds(index, locs)) {
		int loc1 = (locs[0].row * width[index]) + locs[0].col;
		int loc2 = (locs[1].row * width[index]) + locs[1].col;

		sampleDist1 = sqrt(pow((double)(locs[0].row - locs[1].row), 2.0) + pow((locs[0].col - locs[1].col), 2.0)) / distBinWidth;
		sampleColor1 = sqrt(pow((double)(pyr.levels[index](locs[0].row, locs[0].col, RED) - pyr.levels[index](locs[1].row, locs[1].col, RED)), 2.0) +
							pow((double)(pyr.levels[index](locs[0].row, locs[0].col, GREEN) - pyr.levels[index](locs[1].row, locs[1].col, GREEN)), 2.0) +
							pow((double)(pyr.levels[index](locs[0].row, locs[0].col, BLUE) - pyr.levels[index](locs[1].row, locs[1].col, BLUE)), 2.0)) / 255;

		distKernel = exp(sampleDist1 * -distBinWidth);
		colorKernel = exp(sampleColor1 * colorBinWidth); // TODO: Changed from -ve

		ksum.kernelSum = colorKernel;

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
int inValidImageBounds(int index, Location loc[2]) {
	int valid = 0;
	int validCount = 0;

	int i = 0;
	while (i < 2) {
		if (loc[i].row >= 0 && loc[i].row < height[index] && loc[i].col >= 0 && loc[i].col < width[index]) {
			validCount++;
		}
		i++;
	}

	if (validCount == 2) {
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
Bounds getApplicableBounds(Location locs[2], int neighborhood) {
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
	for (int i = 1; i < 2; i++) {
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
void updateApplicableRegion(int index, Bounds bounds, KSumInfo gsum, KSumInfo lsum, KSumInfo csum) {
	int sampleCountLimit = 128; // TODO: changed from 64 (was 30-100 before)

	for (int i = bounds.topLeft.row; i <= bounds.botLeft.row; i++) {
		for (int j = bounds.topLeft.col; j <= bounds.topRight.col; j++) {
			/*
			 * Must be within the image dimensions
			 */
			if (i >= 0 && i < height[index] && j >= 0 && j < width[index]) {
				int ij = (i * width[index]) + j;

				float avg = (gsum.kernelSum + lsum.kernelSum + csum.kernelSum) / 3.0;
//				float avg = csum.kernelSum;

				if (pixels[index][ij].sampleCount < sampleCountLimit) {
					pixels[index][ij].kernelSum += avg;
					pixels[index][ij].sampleCount++;

					/*
					 * Update the pixel entropy every N iterations
					 */
					// Changed iteration count from 32
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
//		double totalMagSum = p->magSum1 * p->magSum2;
		double totalMagSum = 1; // TODO: changed to neglect magnitudes
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

				if (maps[index](i, j) >= 128) {	// TODO: Changed from 217
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

		intensities[i] = (float *) calloc(width[i] * height[i], sizeof(float));
		if (intensities[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the intensity map.\n");
			exit(1);
		}

		lbp[i] = (int *) calloc(width[i] * height[i], sizeof(int));
		if (lbp[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the LBP map.\n");
			exit(1);
		}

		stdDev[i] = (float *) calloc(width[i] * height[i], sizeof(float));
		if (stdDev[i] == NULL) {
			fprintf(stderr, "Cannot allocate memory for the LBP-Quantized map.\n");
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
		free(lbp[i]);
		free(stdDev[i]);
		free(relevantPixels[i]);
	}

}


/**
 * Utility function to get a file name minus the extension
 */
static char *getFileStem(char *filename) {
	char *str = (char*) malloc((strlen(filename) - 3) * sizeof(char));
	char *chPtr = strrchr(filename, '.');
	int index = chPtr - filename;
	int i = 0;;

	while (i < index) {
		str[i] = filename[i];
		i++;
	}
	str[i] = '\0';

	return str;
}


/**
 * Utility function to convert an integer array into a character
 * array
 */
static void convertValuesToChar(const float *srcVal, unsigned char **tgtVal, int rows, int cols) {
	(*tgtVal) = (unsigned char*) calloc(rows * cols, sizeof(unsigned char));

	if (srcVal == NULL) {
		fprintf(stderr, "Empty source val at line: %d in %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if ((*tgtVal) == NULL) {
		fprintf(stderr, "Couldn't allocate memory for tgtVal in %s, line: %d\n", __FILE__, __LINE__);
		exit(1);
	}

	for (int i=0; i < rows; i++) {
		for (int j=0; j < cols; j++) {
			int ij = (i * cols) + j;
			(*tgtVal)[ij] = (unsigned char) ((int)(srcVal[ij] + 0.5));
		}
	}
}


void writeArrayToFile(const float* array, char *textFilename, int index) {
	FILE *textFile;

	if ((textFile = fopen(textFilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing entropies\n",
				textFilename);
		exit(-1);
	}

	for (int i = 0; i < height[index]; i++) {
		for (int j = 0; j < width[index]; j++) {
			int ij = (i * width[index]) + j;
			fprintf(textFile, "%2.4f ", array[ij]);
		}
		fprintf(textFile, "\n");
	}
	fclose(textFile);
}

