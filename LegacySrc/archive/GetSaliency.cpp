/*
 * GetSaliency.cpp
 *
 *  Created on: Jul 17, 2011
 *      Author: kesterduncan
 *
 *  This file computes the saliency of image pixels using
 *  an Iterative Kernel Density Estimation method alongside
 *  kernel descriptors (gradient, local binary pattern, color).
 *
 *  Copyright (c) 2011 Kester Duncan
 *
 *  $Id: GetSaliency$
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "image/RImage.h"
#include "image/LBP.h"
#include "gaussian/gaussian_filter.h"
#include "edge/canny.h"
#include "KDInfo.h"

#define LOG2(X) (log(X) / 0.30103)

/** Global variables */
int 	width, height; 			// Dimensions of image
float	*intensities;			// Image intensity
float 	*gradDirections; 		// Gradient directions map
float 	*magnitudes; 			// Edge magnitudes map
int		*lbp;					// LBP map
float	*stdDev;				// Standard deviation of pixel values in a 3x3 neighborhood around a pixel
RImage 	src; 					// Input image
RImage 	map; 					// Saliency map
KDInfo	*pixels;				// Stores relevant information for pixels
const int M = 5;				// Dimension of the Local Pixel Neighborhood


/** Function Prototypes */
void 			magnitudeQuantization();
void 			getStandardDeviations();
void			computeKernelSaliency(int percentage);
void 			postProcessMap(float sigma);
static void 	allocateMemory();
static char 	*getFileStem(char *filename);
static void 	convertValuesToChar(const float *srcVal, unsigned char **tgtVal, int rows, int cols);
void			getIntensities();
void 			updateApplicableRegion(Bounds bounds, KSumInfo gsum, KSumInfo lsum, KSumInfo csum);
Bounds 			getApplicableBounds(Location locs[2], int neighborhood);
KSumInfo 		calculateGradKernelSum (Location locs[2]);
KSumInfo 		calculateLBPKernelSum (Location locs[2]);
KSumInfo 		calculateColorKernelSum (Location locs[2]);
void 			updateEntropy(KDInfo *p);
void 			updateSaliencyMap(int index);
int 			inValidImageBounds(Location loc[2]);

using namespace std;

int main(int argc, char *argv[]) {
	time_t	start, stop;
	char	imageFileName[100];
	char 	mapFileName[100];
	char 	entropyFilename[100];
	char	fileStem[100];
	int 	window; 			// neighborhood dimension
	int 	dimension; 			// bin dimension
	int 	percentage; 		// sampling percentage
	unsigned char *imgVals;
	LBP localBinPat(32, 4);

	if (argc != 3) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./GetSaliency [IMAGE] [SAMPLING PERCENTAGE]\n");
		fprintf(stderr, "\te.g. ./GetSaliency ads.pgm 25\n");
		exit(-1);
	}

	/* Mark the start time for all processing */
	start = time(NULL);

	strcpy(imageFileName, argv[1]);
	strcpy(fileStem, getFileStem(imageFileName));
	printf("Evaluating saliency of image: %s\n", imageFileName);

	percentage = atoi(argv[2]);

	src.read(imageFileName);
	if (!src.isColorImage) {
		src.convertToColor();
	}

	map.initialize(src.numRows, src.numCols);
	width = src.numCols;
	height = src.numRows;
	allocateMemory();

	getIntensities();
	convertValuesToChar(intensities, &imgVals, height, width);
	canny(imgVals, height, width, 1.20, gradDirections, magnitudes);
	localBinPat.getOriginalLBP(src, lbp);
	magnitudeQuantization();
	getStandardDeviations();

	/*
	 * Calculate iterative saliency
	 */
	initializeKDInfo(pixels, (width * height));
	computeKernelSaliency(percentage);

	sprintf(mapFileName, "%s_map_it.ppm", fileStem);

//	postProcessMap(10);
	/*
	 * Make composite map
	 */
	RImage comp;
	comp.isColorImage = 1;
	comp.initialize(map.numRows, (map.numCols * 2));

	for (int m=0; m < map.numRows; m++) {
		for (int n=0; n < map.numCols; n++) {
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
//	map.makeThermalColorImage();
	for (int m=0; m < map.numRows; m++) {
		for (int n=0; n < map.numCols; n++) {
			int col = n + map.numCols;
			comp(m, col, RED) = map(m, n, RED);
			comp(m, col, GREEN) = map(m, n, RED);
			comp(m, col, BLUE) = map(m, n, RED);

		}
	}

#else
	for (int m=0; m < map.numRows; m++) {
		for (int n=0; n < map.numCols; n++) {
			int col = n + map.numCols;
			double value = map(m, n, RED);
			if (value > 128) {
				comp(m, col, RED) = src(m, n, RED);
				comp(m, col, GREEN) = src(m, n, RED);
				comp(m, col, BLUE) = src(m, n, RED);
			}
		}
	}
#endif

	comp.save(mapFileName);
	comp.flushMemory();
	map.flushMemory();
	src.flushMemory();
	free(imgVals);
	free(gradDirections);
	free(magnitudes);
	free(pixels);
	free(intensities);
	free(lbp);
	free(stdDev);

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
void magnitudeQuantization() {
	if (magnitudes != NULL) {
		float magSum = 0;

		for (int k = 0; k < (width * height); k++) {
			magSum += pow((double) magnitudes[k], 2.0) + 0.05;
		}

		double denom = sqrt((double) magSum);
		for (int k = 0; k < (width * height); k++) {
			magnitudes[k] = magnitudes[k] / denom;
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
void getStandardDeviations() {
	if (stdDev != NULL) {
		double sdSum = 0;

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				float sum = 0;
				float avg = 0;
				float sd = 0;
				int count = 0;
				int ij = (i * width) + j;
				int neighbor = 1;

				for (int m = -neighbor; m <= neighbor; m++) {
					for (int n = -neighbor; n <= neighbor; n++) {
						int row = i + m;
						int col = j + n;

						if (row >= 0 && col >= 0 && row < height && col < width) {
							count++;
							sum += src(row, col, RED);
						}
					}
				}

				avg = sum / count;
				sum = 0;
				for (int m = -neighbor; m <= neighbor; m++) {
					for (int n = -neighbor; n <= neighbor; n++) {
						int row = i + m;
						int col = j + n;

						if (row >= 0 && col >= 0 && row < height && col < width) {
							count++;
							sum += pow((double)(src(row, col, RED) - avg), 2.0);
						}
					}
				}
				sd = sum / (count - 1) ;
				stdDev[ij] = sd;
				sdSum += pow((double) sd, 2.0) + 0.05;;
			}
		}

		double denom = sqrt((double) sdSum);
		for (int k = 0; k < (width * height); k++) {
			stdDev[k] = stdDev[k] / denom;
		}

	} else {
		fprintf(stderr, "ERROR: The lbp map is NULL.\n");
		exit(1);
	}
}


/**
 * Simply place the image intensities into an array
 */
void getIntensities() {
	if (intensities != NULL) {
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				intensities[ij] = (src(i, j, RED) + src(i, j, GREEN) + src(i, j, BLUE)) / 3.0;
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
	gaussianSmooth((unsigned int*) map.getChannel(RED), map.numRows, map.numCols, sigma, &filteredImage);
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			int ij = (i * width) + j;
			map(i, j, RED) = (int) (filteredImage[ij]);
			if (map.isColorImage) {
				map(i, j, GREEN) = (int) (filteredImage[ij]);
				map(i, j, BLUE) = (int) (filteredImage[ij]);
			}
		}
	}

	free(filteredImage);
	filteredImage = NULL;

	/*
	 * Smoothing another time effectively gets rid of
	 * unsightly rings caused by the Gaussian
	 */
	gaussianSmooth((unsigned int*) map.getChannel(RED), map.numRows, map.numCols, sig2, &filteredImage);
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			int ij = (i * width) + j;
			map(i, j, RED) = (int) (filteredImage[ij]);
			if (map.isColorImage) {
				map(i, j, GREEN) = (int) (filteredImage[ij]);
				map(i, j, BLUE) = (int) (filteredImage[ij]);
			}
		}
	}

	free(filteredImage);

	float maxVal = -999;
	float minVal = 999;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (map(i, j, RED) > maxVal) {
				maxVal = map(i, j);
			}
			if (map(i, j, RED) < minVal) {
				minVal = map(i, j);
			}
		}
	}

	/* Normalization */
	float sum = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			map(i, j, RED) = ((((map(i, j, RED) - minVal) / (maxVal - minVal)) * 255.0));

			if (map.isColorImage) {
				map(i, j, GREEN) = ((((map(i, j, GREEN) - minVal) / (maxVal - minVal)) * 255.0));
				map(i, j, BLUE) = ((((map(i, j, BLUE) - minVal) / (maxVal - minVal)) * 255.0));
			}
			sum += map(i, j, RED);
		}
	}
}


/**
 * Process the saliency of each pixel in an image using an Iterative Kernel
 * Approach for calculating the Kernel Density Estimate of the local neighborhood
 * distribution, then calculating the entropy of that distribution for
 * determining the saliency.
 */
void computeKernelSaliency(int percentage) {
	int pixelSamples;		// Maximum number of pixels sampled randomly
	int index = 1;			// Stores the index used for the intermediate maps
	int mainCount = 0;		// Counter for the main loop
	int mSquared = M * M;
	int halfM = 2;

	pixelSamples = (percentage / 100.0) * width * height * mSquared;

	/* Ensure that the dimensions are valid before any processing */
	if (width <= 0 || height <= 0) {
		fprintf(stderr, "Error: Image dimensions are invalid.\n");
		exit(1);
	}

	/* Set the seed for randomization */
	srand(time(NULL));

	while (mainCount < pixelSamples) {
		KSumInfo gsum, lsum, csum;
		Bounds bounds;
		Location locs[2];			// Array to store randomly selected pixel locations

		/*
		 * Randomly select location of first pixel
		 */
		locs[0].row = rand() % height;
		locs[0].col = rand() % width;

		/*
		 * The other pixel must be selected in the neighborhood
		 * of pixel 1
		 */
		locs[1].row = (rand() % M) + (locs[0].row - halfM);
		locs[1].col = (rand() % M) + (locs[0].col - halfM);

		if (locs[1].row >= 0 && locs[1].row < height &&
			locs[1].col >= 0 && locs[1].col < width) {
			gsum = calculateGradKernelSum(locs);
			lsum = calculateLBPKernelSum(locs);
			csum = calculateColorKernelSum(locs);
			bounds = getApplicableBounds(locs, M);
			updateApplicableRegion(bounds, gsum, lsum, csum);

			mainCount++;
		}
	}

	updateSaliencyMap(index);

}


/**
 * Calculate the gradient kernel sum from the contribution of the
 * pixels in 'randomPixels'
 */
KSumInfo calculateGradKernelSum (Location locs[2]) {
	KSumInfo ksum;
	double sampleDist1, sampleAngle1, sampleMag1, sampleMag2;
	double distKernel = 0, angleKernel = 0;
	double binDimension = 10.0;
	double distBinWidth = sqrt(pow((double) M, 2.0) + pow((double) M, 2.0)) / binDimension;
	double angleBinWidth = 3.142 / binDimension;
	double twoPI = 6.283;

	double pNorm = (2.5066 * distBinWidth);
	double gNorm = (2.5066 * angleBinWidth);

	if (inValidImageBounds(locs)) {
		int loc1 = (locs[0].row * width) + locs[0].col;
		int loc2 = (locs[1].row * width) + locs[1].col;

		sampleDist1 = sqrt(pow((double)(locs[0].row - locs[1].row), 2.0) + pow((locs[0].col - locs[1].col), 2.0)) / distBinWidth;

		sampleAngle1 = sqrt(pow(sin((double)gradDirections[loc1]) - sin((double)gradDirections[loc2]), 2.0) +
						   pow(cos((double)gradDirections[loc1]) - cos((double)gradDirections[loc2]), 2.0));

		sampleMag1 = magnitudes[loc1];
		sampleMag2 = magnitudes[loc2];

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
KSumInfo calculateLBPKernelSum (Location locs[2]) {
	KSumInfo ksum;
	double sampleDist1, sampleLBP1, sampleSD1, sampleSD2;
	double distKernel = 0, lbpKernel = 0;
	double binDimension = 10.0;
	double distBinWidth = sqrt(pow((double) M, 2.0) + pow((double) M, 2.0)) / binDimension;
	double lbpBinWidth = 255 / binDimension;
	double twoPI = 6.283;

	if (inValidImageBounds(locs)) {
		int loc1 = (locs[0].row * width) + locs[0].col;
		int loc2 = (locs[1].row * width) + locs[1].col;

		sampleDist1 = sqrt(pow((double)(locs[0].row - locs[1].row), 2.0) + pow((locs[0].col - locs[1].col), 2.0)) / distBinWidth;

		sampleLBP1 = sqrt(pow((double)(lbp[loc1] - lbp[loc2]), 2.0)) / 255;
		sampleSD1 = stdDev[loc1];
		sampleSD2 = stdDev[loc2];

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
KSumInfo calculateColorKernelSum (Location locs[2]) {
	KSumInfo ksum;
	double sampleDist1, sampleColor1;
	double distKernel = 0, colorKernel = 0;
	double binDimension = 10.0;
	double distBinWidth = sqrt(pow((double) M, 2.0) + pow((double) M, 2.0)) / binDimension;
	double colorBinWidth = 255 / binDimension;
	double twoPI = 6.283;

	if (inValidImageBounds(locs)) {
		int loc1 = (locs[0].row * width) + locs[0].col;
		int loc2 = (locs[1].row * width) + locs[1].col;

		sampleDist1 = sqrt(pow((double)(locs[0].row - locs[1].row), 2.0) + pow((locs[0].col - locs[1].col), 2.0)) / distBinWidth;
		sampleColor1 = sqrt(pow((double)(src(locs[0].row, locs[0].col, RED) - src(locs[1].row, locs[1].col, RED)), 2.0) +
							pow((double)(src(locs[0].row, locs[0].col, GREEN) - src(locs[1].row, locs[1].col, GREEN)), 2.0) +
							pow((double)(src(locs[0].row, locs[0].col, BLUE) - src(locs[1].row, locs[1].col, BLUE)), 2.0)) / 255;

		distKernel = exp(sampleDist1 * -distBinWidth);
		colorKernel = exp(sampleColor1 / colorBinWidth); // TODO: Changed from -ve

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
int inValidImageBounds(Location loc[2]) {
	int valid = 0;
	int validCount = 0;

	int i = 0;
	while (i < 2) {
		if (loc[i].row >= 0 && loc[i].row < height && loc[i].col >= 0 && loc[i].col < width) {
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
void updateApplicableRegion(Bounds bounds, KSumInfo gsum, KSumInfo lsum, KSumInfo csum) {
	int sampleCountLimit = 64; // TODO: change back to 30 - 100

	for (int i = bounds.topLeft.row; i <= bounds.botLeft.row; i++) {
		for (int j = bounds.topLeft.col; j <= bounds.topRight.col; j++) {
			/*
			 * Must be within the image dimensions
			 */
			if (i >= 0 && i < height && j >= 0 && j < width) {
				int ij = (i * width) + j;

				float avg = csum.kernelSum;

				if (pixels[ij].sampleCount < sampleCountLimit) {
					pixels[ij].kernelSum += avg;
					pixels[ij].sampleCount++;

					/*
					 * Update the pixel entropy every N iterations
					 */
					// Changed iteration count from 32
					if (((pixels[ij].sampleCount + 1) % 32) == 0) {
						updateEntropy(&pixels[ij]);
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

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				if (pixels[ij].entropy > maxEntropy) {
					maxEntropy = pixels[ij].entropy;
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
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				if (pixels[ij].entropy == -1) {
					pixels[ij].entropy = maxEntropy;
				}

				if (pixels[ij].entropy < minEntropy) {
					minEntropy = pixels[ij].entropy;
				}
			}
		}

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				pixels[ij].entropy -= minEntropy;
			}
		}

		/* Also adjust the maximum entropy */
		 maxEntropy -= minEntropy;

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;

				if (maxEntropy > 0) {
					map(i, j) = (int) (255.0 - ((pixels[ij].entropy / maxEntropy) * 255.0));
				} else {
					map(i, j) = 0;
				}
			}
		}

	} else {
		printf("Error: Pixels array is NULL.\n");
		exit(1);
	}
}


/**
 * Allocates memory for the global arrays
 * i.e. gradDirections, magnitudes, and entropies
 */
static void allocateMemory() {
	pixels = (KDInfo *) calloc(width * height, sizeof(KDInfo));
	if (pixels == NULL) {
		fprintf(stderr, "Cannot allocate memory for the pixels map.\n");
		exit(1);
	}

	gradDirections = (float *) calloc(width * height, sizeof(float));
	if (gradDirections == NULL) {
		fprintf(stderr, "Cannot allocate memory for the gradDirections map.\n");
		exit(1);
	}

	magnitudes = (float *) calloc(width * height, sizeof(float));
	if (magnitudes == NULL) {
		fprintf(stderr, "Cannot allocate memory for the magnitudes map.\n");
		exit(1);
	}

	intensities = (float *) calloc(width * height, sizeof(float));
	if (intensities == NULL) {
		fprintf(stderr, "Cannot allocate memory for the intensity map.\n");
		exit(1);
	}

	lbp = (int *) calloc(width * height, sizeof(int));
	if (lbp == NULL) {
		fprintf(stderr, "Cannot allocate memory for the LBP map.\n");
		exit(1);
	}

	stdDev = (float *) calloc(width * height, sizeof(float));
	if (stdDev == NULL) {
		fprintf(stderr, "Cannot allocate memory for the LBP-Quantized map.\n");
		exit(1);
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
		fprintf(stderr, "Couldn't allocate memore for tgtVal in %s, line: %d\n", __FILE__, __LINE__);
		exit(1);
	}

	for (int i=0; i < rows; i++) {
		for (int j=0; j < cols; j++) {
			int ij = (i * cols) + j;
			(*tgtVal)[ij] = (unsigned char) ((int)(srcVal[ij] + 0.5));
		}
	}
}
