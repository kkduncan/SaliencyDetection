/*
 * ComputeRelativeSaliency.cpp
 *
 *  Created on: March 23, 2010
 *      Author: kesterduncan
 */

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "REM.h"
#include "image.h"
#include "GroupParallel.h"

/* Globals */
/** Width and height of image */
int width, height;

/** Stores the gradient directions for pixels */
float *gradDirections;

/** Stores the edge magnitudes for pixels */
int *magnitudes;

/** Stores the entropies */
double *entropies;

/** Items related to group windows */
Image src, rel, shan;
REM avg;

/* Prototypes */
void getWidthAndHeight(char *filename);
void findEdges(char imageFile[100], char outFileStem[100]);
int readGradientFile(char filename[100]);
int readMagnitudesFile(char filename[100]);
void eliminateLowMagnitudes();
void processPixelsForEntropy(int window, int dimension,	int percentage, char *relFilename);
int getint(FILE *fp);

using namespace std;

int main(int argc, char *argv[]) {
	char imageFileName[100];
	char outFileStem[100];
	char resultsFileName[100];
	char shanFileName[100];
	char renyiFilename[100];
	char shannonFilename[100];
	char str[50];
	int window;
	int dimension;
	int percentage;

	if (argc != 6) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./ComputeSaliency [IMAGE] [IMAGE_NAME] [WINDOW] [SQUARE_HISTOGRAM_DIMENSION] [SAMPLING]\n");
		fprintf(stderr, "\te.g. ./ComputeSaliency ads.pgm ads 11 10 5\n");
		exit(-1);
	}

	strcpy(imageFileName, argv[1]);
	strcpy(outFileStem, argv[2]);
	window = atoi(argv[3]);
	dimension = atoi(argv[4]);
	percentage = atoi(argv[5]);

	getWidthAndHeight(imageFileName);
	findEdges(imageFileName, outFileStem);

	sprintf(str, "texts/%s_GD.txt", outFileStem);
	readGradientFile(str);

	sprintf(str, "texts/%s_MAG.txt", outFileStem);
	readMagnitudesFile(str);
	eliminateLowMagnitudes();

	src.read(imageFileName);
	rel.initialize(src.NR, src.NC);
	shan.initialize(src.NR, src.NC);

	sprintf(resultsFileName, "maps/%s_rel_%d.pgm", outFileStem, window);
	sprintf(renyiFilename, "texts/%s_rel_%d.txt", outFileStem, window);

	processPixelsForEntropy(window, dimension, percentage, renyiFilename);
	rel.save(resultsFileName);

	free(gradDirections);
	free(magnitudes);
	free(entropies);

	return 0;
}

void getWidthAndHeight(char *filename) {
	FILE *imgFile;
	char str[50];

	if ((imgFile = fopen(filename, "r")) == NULL) {
		printf("Cannot open image file %s\n", filename);
		exit(-1);
	}

	fscanf(imgFile, "%s", str); // Get the P5 string
	width = getint(imgFile);
	height = getint(imgFile);

	if (width <= 0 || height <= 0) {
		printf("Image file %s is in an invalid format\n", filename);
		exit(-1);
	}

	fclose(imgFile);
}

void findEdges(char imageFile[100], char outFileStem[100]) {
	char command[500];
	sprintf(command, "./CannyEdge %s 1.20 0.20 0.70 %s", imageFile, outFileStem);
	system(command);

}

int readGradientFile(char filename[100]) {
	FILE *gradFile;

	if ((gradFile = fopen(filename, "r")) == NULL) {
		printf("Cannot open gradient text file %s\n", filename);
		exit(-1);
	}
	gradDirections = (float *) malloc(width * height * sizeof(float));

	for (int k = 0; k < (width * height); k++) {
		*(gradDirections + k) = 0;
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			float value;
			fscanf(gradFile, "%f ", &value);
			*(gradDirections + ij) = value;
		}
	}

	fclose(gradFile);
	return 0;
}

int readMagnitudesFile(char filename[100]) {
	FILE *magFile;

	if ((magFile = fopen(filename, "r")) == NULL) {
		printf("Cannot open gradient text file %s\n", filename);
		exit(-1);
	}
	magnitudes = (int *) malloc(width * height * sizeof(int));

	for (int k = 0; k < (width * height); k++) {
		*(magnitudes + k) = 0;
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			int value;
			fscanf(magFile, "%d ", &value);
			*(magnitudes + ij) = value;
		}
	}

	fclose(magFile);
	return 0;
}

void eliminateLowMagnitudes() {
	int maxMag = -999;

	for (int k = 0; k < (width * height); k++) {
		if (magnitudes[k] > maxMag) {
			int temp = maxMag;
			maxMag = max(temp, magnitudes[k]);
		}
	}

	for (int k = 0; k < (width * height); k++) {
		if (magnitudes[k] < (maxMag / 3)) {
			magnitudes[k] = 1;
		} else if (magnitudes[k] < (maxMag / 2)) {
			magnitudes[k] = 5;
		} else if (magnitudes[k] < (maxMag / 1.5)) {
			magnitudes[k] = 50;
		} else {
			magnitudes[k] = 1000;
		}
	}
}

void processPixelsForEntropy(int window, int dimension,	int percentage, char *relFilename) {
	FILE *relTextFile;

	if ((relTextFile = fopen(relFilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing renyi entropies\n",
				relFilename);
		exit(-1);
	}

	avg.init(1, window, window, dimension, percentage);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int neighbor = (int) window / 2;
			float entropy = 0;

			REM group;
			group.init(1, window, window, dimension, percentage);

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int row = i + m;
					int col = j + n;

					/* Check if neighborhood pixels are within image boundary */
					if (row > 0 && col > 0 && row < rel.NR && col < rel.NC) {
						int mn = (row * width) + col;
						group.addElement(row, col, gradDirections[mn], magnitudes[mn]);
					} else {
						group.addElement(i, j, M_PI, 0);
					}
				}
			}

			group.buildRelativeHistogram(avg.hist);
		}
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int neighbor = (int) window / 2;
			float entropy = 0;

			REM group;
			group.init(1, window, window, dimension, percentage);

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int row = i + m;
					int col = j + n;

					/* Check if neighborhood pixels are within image boundary */
					if (row > 0 && col > 0 && row < rel.NR && col < rel.NC) {
						int mn = (row * width) + col;
						group.addElement(row, col, gradDirections[mn], magnitudes[mn]);
					} else {
						group.addElement(i, j, M_PI, 0);
					}
				}
			}

			entropy = group.calculateRelativeEntropy(avg.hist);
			(entropy < 0) ? entropy = 0.999 : entropy = entropy;
			rel(i, j) = 255 - (int) (entropy * 255);
			fprintf(relTextFile, "%2.4f ", entropy);

		}
		fprintf(relTextFile, "\n");
	}

	fclose(relTextFile);
	relTextFile = NULL;

}

/**
 * Reads an integer value from the file pointed to by fp
 */
int getint(FILE *fp) {
	int item, i, flag;

	/* skip forward to start of next number */
	item = getc(fp);
	flag = 1;
	do {

		if (item == '#') { /* comment */
			while (item != '\n' && item != EOF)
				item = getc(fp);
		}

		if (item == EOF)
			return 0;
		if (item >= '0' && item <= '9') {
			flag = 0;
			break;
		}

		/* illegal values */
		if (item != '\t' && item != '\r' && item != '\n' && item != ',')
			return (-1);

		item = getc(fp);
	} while (flag == 1);

	/* get the number */
	i = 0;
	flag = 1;
	do {
		i = (i * 10) + (item - '0');
		item = getc(fp);
		if (item < '0' || item > '9') {
			flag = 0;
			break;
		}
		if (item == EOF)
			break;
	} while (flag == 1);

	return i;
}

