/*
 * ComputeSaliency.cpp
 *
 *  Created on: Jan 30, 2010
 *      Author: kesterduncan
 *
 *  This file takes an image as input, computes the saliency of
 *  image pixels, and outputs a saliency map of the image.
 *
 *  Copyright (c) 2010 Kester Duncan
 *
 *  $Id: ComputeSaliency$
 */
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "image/RImage.h"
#include "gaussian/gaussian_filter.h"
#include "edge/canny.h"
#include "gaussian/GaussianPyramid.h"
#include "REM.h"


/** Global variables */
int 	width, height; 			// Dimensions of image
float 	*gradDirections; 		// Gradient directions map
float 	*magnitudes; 			// Edge magnitudes map
double 	*entropies; 			// Stores entropy for each pixel
RImage 	src; 					// Input image
RImage 	map; 					// Saliency map
GaussianPyramid pyr;			// Gaussian pyramid of the input image

//#define AVG

/** Function Prototypes */
void 			magnitudeQuantization();
void 			computeSaliency(int window, int dimension, int percentage);
void 			computeSaliencyPyramid(int level, int window, int dimension, int percentage);
void			computeAverageSaliency();
void 			postProcessMap(float sigma);
static void 	allocateMemory();
void 			writeEntropiesToFile(char *textFile);
void 			writeGradientsToFile(char *textFilename);
void 			writeMagnitudesToFile(char *textFilename);
static char 	*getFileStem(const char *filename);
static void 	convertValuesToChar(const int *srcVal, unsigned char **tgtVal, int rows, int cols);

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

	if (argc != 5) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./ComputeSaliency [IMAGE] [WINDOW] [HISTOGRAM_DIMENSION] [SAMPLING PERCENTAGE]\n");
		fprintf(stderr, "\te.g. ./ComputeSaliency ads.pgm 11 10 5\n");
		exit(-1);
	}

	/* Mark the start time for all processing */
	start = time(NULL);

	strcpy(imageFileName, argv[1]);
	strcpy(fileStem, getFileStem(imageFileName));
	printf("Evaluating saliency of image: %s\n", imageFileName);

	window = 		atoi(argv[2]);
	dimension = 	atoi(argv[3]);
	percentage = 	atoi(argv[4]);

	src.read(imageFileName);
	map.initialize(src.NR, src.NC);

	/*
	 * If AVG is defined above, we calculate salieny by the averaging method,
	 * which performs Gaussian Pyramid processing on the input, calculates
	 * the saliency of each image, and computes a weighted-avearge in order
	 * to create the saliency map.
	 *
	 * If AVG is not defined above, we simply process the saliency of the
	 * input image.
	 *
	 */
#ifdef AVG

	pyr.setBaseLevel(src);
	pyr.reduceAll();
	printf("\tCreated the pyramid...\n");

	for (int i=4; i >= 0; i--) {
		width = pyr.levels[i].NC;
		height = pyr.levels[i].NR;

		allocateMemory();
		convertValuesToChar(pyr.levels[i].getChannel(RED), &imgVals, height, width);
		canny(imgVals, height, width, 1.20, gradDirections, magnitudes);

		magnitudeQuantization();
		computeSaliencyPyramid(i, window, dimension, percentage);

		free(imgVals);
		free(gradDirections);
		free(magnitudes);
		free(entropies);

	}

	printf("\tComputed the Pyramid saliency...\n");

	pyr.expandAll();
	printf("\tExpanded the Saliency Pyramid...\n");

	width = src.NC;
	height = src.NR;

	computeAverageSaliency();
	postProcessMap(15);
	printf("\tAveraged the saliencies... outputting saliency map...\n");

#else

	width = src.NC;
	height = src.NR;
	allocateMemory();

	convertValuesToChar(src.getChannel(RED), &imgVals, height, width);
	canny(imgVals, height, width, 1.20, gradDirections, magnitudes);

	magnitudeQuantization();
	computeSaliency(window, dimension, percentage);
	postProcessMap(30);

#endif

	sprintf(mapFileName, "%s_map_%d_%d.pgm", fileStem, window, dimension);
	sprintf(entropyFilename, "%s_entropy_%d_%d.txt", fileStem, window, dimension);
	map.save(mapFileName);


#ifndef AVG
	free(imgVals);
	free(gradDirections);
	free(magnitudes);
	free(entropies);

#endif

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
		double maxMag = -999;

		for (int k = 0; k < (width * height); k++) {
			if (*(magnitudes + k) > maxMag) {
				maxMag = magnitudes[k];
			}
		}

		/* Thresholds */
		double t1 = (double) maxMag / 3.0;
		double t2 = (double) maxMag / 2.0;
		double t3 = (double) maxMag / 1.5;

		for (int k = 0; k < (width * height); k++) {

			if (magnitudes[k] >= 0 && magnitudes[k] < t1) {
				magnitudes[k] = 0.0005;

			} else if (magnitudes[k] >= t1 && magnitudes[k] < t2) {
				magnitudes[k] = 0.0025;

			} else if (magnitudes[k] >= t2 && magnitudes[k] < t3) {
				magnitudes[k] = 0.0125;

			} else if (magnitudes[k] >= t3){
				magnitudes[k] = 1.0;

			}
		}
	} else {
		fprintf(stderr, "ERROR: The magnitudes map is NULL.\n");
		exit(1);
	}
}


/**
 * Process the saliency of each pixel in an image
 */
void computeSaliency(int window, int dimension, int percentage) {

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int neighbor = (int) window / 2;
			int ij = (i * width) + j;
			double entropy = 0;
			double ratio = 0;
			double old_entropy = 1;

			REM group(window, window, dimension, percentage);

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int row = i + m;
					int col = j + n;

					if (row > 0 && col > 0 && row < map.NR && col < map.NC) {
						int mn = (row * width) + col;
						group.addElement(row, col, gradDirections[mn], magnitudes[mn]);
					} else {
						group.addElement(i, j, (M_PI / (m + 1)), 0);
					}
				}
			}

			if (percentage != 100) {
				group.copyElementsToArray();
				do {
					entropy = group.estimateEntropyFromSample();

					/* Hack to handle constant intensity regions */
					if (entropy != -1) {
						ratio = fabs(old_entropy - entropy) / old_entropy;
						old_entropy = entropy;

					} else {
						ratio = 0.00001;
					}

				} while (ratio > 0.01);

			} else {
				entropy = group.estimateEntropy();
			}

			entropies[ij] = entropy;
		}
	}

	/*
	 * Post Processing for proper visualization
	 */
	double maxEntropy = -999.0;
	double minEntropy = 999.0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			if (entropies[ij] > maxEntropy) {
				maxEntropy = entropies[ij];
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
			if (entropies[ij] == -1) {
				entropies[ij] = maxEntropy;
			}
			if (entropies[ij] < minEntropy) {
				minEntropy = entropies[ij];
			}
		}
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			if (minEntropy < 0) {
				entropies[ij] += fabs(minEntropy);
			} else {
				entropies[ij] -= minEntropy;
			}
		}
	}

	/*
	 * Also adjust the maximum entropy
	 */
	(minEntropy < 0) ? maxEntropy += fabs(minEntropy) : maxEntropy -= minEntropy;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			if (maxEntropy > 0 && minEntropy > 0) {
				map(i, j) = (int) (255.0 - ((entropies[ij] / maxEntropy) * 255.0));
			} else {
				map(i, j) = 0;
			}
		}
	}
}


/**
 * Process the saliency of each pixel in a pyramid image
 */
void computeSaliencyPyramid(int level, int window, int dimension, int percentage) {

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int neighbor = (int) window / 2;
			int ij = (i * width) + j;
			double entropy = 0;
			double ratio = 0;
			double old_entropy = 1;

			REM group(window, window, dimension, percentage);

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int row = i + m;
					int col = j + n;

					if (row > 0 && col > 0 && row < map.NR && col < map.NC) {
						int mn = (row * width) + col;
						group.addElement(row, col, gradDirections[mn], magnitudes[mn]);
					} else {
						group.addElement(i, j, (M_PI / (m + 1)), 0);
					}
				}
			}

			if (percentage != 100) {
				group.copyElementsToArray();
				do {
					entropy = group.estimateEntropyFromSample();

					/* Hack to handle constant intensity regions */
					if (entropy != -1) {
						ratio = fabs(old_entropy - entropy) / old_entropy;
						old_entropy = entropy;

					} else {
						ratio = 0.0001;
					}

				} while (ratio > 0.01);

			} else {
				entropy = group.estimateEntropy();
			}

			entropies[ij] = entropy;
		}
	}

	/*
	 * Post Processing for proper visualization
	 */
	double maxEntropy = -999.0;
	double minEntropy = 999.0;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			if (entropies[ij] > maxEntropy) {
				maxEntropy = entropies[ij];
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
			if (entropies[ij] == -1) {
				entropies[ij] = maxEntropy;
			}
			if (entropies[ij] < minEntropy) {
				minEntropy = entropies[ij];
			}
		}
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			if (minEntropy < 0) {
				entropies[ij] += fabs(minEntropy);
			} else {
				entropies[ij] -= minEntropy;
			}
		}
	}

	/*
	 * Also adjust the maximum entropy
	 */
	(minEntropy < 0) ? maxEntropy += fabs(minEntropy) : maxEntropy -= minEntropy;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			if (maxEntropy > 0 && minEntropy > 0) {
				pyr.levels[level](i, j) = (int) (255.0 - ((entropies[ij] / maxEntropy) * 255.0));
			} else {
				pyr.levels[level](i, j) = 0;
			}
		}
	}
}


/**
 * Calculate the weighted-avergage saliency of the pyramid maps
 */
void computeAverageSaliency() {
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			double sum = 0;
//			double weights[] = {0.35, 0.25, 0.17, 0.13, 0.10};

			for (int k=0; k < 5; k++) {
				sum += (pyr.levels[k](i, j));
			}

			sum /= 5.0;
			map(i, j, RED) = (int) ((sum > 255) ? 255 : sum);
		}
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
	gaussianSmooth((unsigned int*) map.getChannel(RED), map.NR, map.NC, sigma, &filteredImage);
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			int ij = (i * width) + j;
			map(i, j, RED) = (int) (filteredImage[ij]);
		}
	}

	free(filteredImage);
	filteredImage = NULL;

	/*
	 * Smoothing another time effectively gets rid of
	 * unsightly rings caused by the Gaussian
	 */
	gaussianSmooth((unsigned int*) map.getChannel(RED), map.NR, map.NC, sig2, &filteredImage);
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			int ij = (i * width) + j;
			map(i, j, RED) = (int) (filteredImage[ij]);
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

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			map(i, j, RED) = ((((map(i, j, RED) - minVal) / (maxVal - minVal)) * 255.0));
		}
	}
}


/**
 * Allocates memory for the global arrays
 * i.e. gradDirections, magnitudes, and entropies
 */
static void allocateMemory() {
	gradDirections = (float *) malloc(width * height * sizeof(float));
	if (gradDirections == NULL) {
		fprintf(stderr, "Cannot allocate memory for the gradDirections map.\n");
		exit(1);
	}

	magnitudes = (float *) malloc(width * height * sizeof(float));
	if (magnitudes == NULL) {
		fprintf(stderr, "Cannot allocate memory for the magnitudes map.\n");
		exit(1);
	}

	entropies = (double *) malloc(width * height * sizeof(double));
	if (entropies == NULL) {
		fprintf(stderr, "Cannot allocate memory for entropies array.\n");
		exit(1);
	}
}


/**
 * Utility function to get a file name minus the extension
 */
static char *getFileStem(const char *filename) {
	char *str = (char*) malloc((strlen(filename)-2) * sizeof(char));
	char *chPtr = strrchr(filename, '.');
	int index = chPtr - filename;
	int i = 0;;

	while (i < index) {
		str[i] = filename[i];
		i++;
	}
	str[i+1] = '\0';

	return str;
}


/**
 * Utility function to convert an integer array into a character
 * array
 */
static void convertValuesToChar(const int *srcVal, unsigned char **tgtVal, int rows, int cols) {
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


/**
 * Utility function to write entropy values for the image to a text file
 */
void writeEntropiesToFile(char *textFilename) {
	FILE *textFile;

	if ((textFile = fopen(textFilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing entropies\n",
				textFilename);
		exit(-1);
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			fprintf(textFile, "%2.4f ", pixels[ij].entropy);
		}
		fprintf(textFile, "\n");
	}
	fclose(textFile);
}


/**
 * Utility function to write gradient values for the image to a text file
 */
void writeGradientsToFile(char *textFilename) {
	FILE *textFile;

	if ((textFile = fopen(textFilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing entropies\n",
				textFilename);
		exit(-1);
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			fprintf(textFile, "%2.4f ", gradDirections[ij]);
		}
		fprintf(textFile, "\n");
	}
	fclose(textFile);
}


/**
 * Utility function to write magnitude values for the image to a text file
 */
void writeMagnitudesToFile(char *textFilename) {
	FILE *textFile;

	if ((textFile = fopen(textFilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing entropies\n",
				textFilename);
		exit(-1);
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			fprintf(textFile, "%2.5f ", magnitudes[ij]);
		}
		fprintf(textFile, "\n");
	}
	fclose(textFile);
}
