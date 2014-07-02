/*
 * GetScaleSpace.cpp
 *
 *  Created on: Feb 13, 2010
 *      Author: kesterduncan
 */

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "REM.h"
#include "image.h"


/********************************* GLOBALS **********************************/
int width, height;
float *gradDirections;
int *magnitudes;
Image src;


/******************************** PROTOTYPES *********************************/
void getWidthAndHeight(char *filename);
void findEdges(char imageFile[100], char outFileStem[100]);
int readGradientFile(char filename[100]);
int readMagnitudesFile(char filename[100]);
void eliminateLowMagnitudes();
int getint(FILE *fp);

/********************************** MAIN *************************************/
using namespace std;

int main(int argc, char *argv[]) {
	char imageFileName[100];
	char outFileStem[100];
	char str[50];
	FILE *textFile;

	/*
	 * Check command line arguments
	 */
	if (argc != 3) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./ScaleRelations [IMAGE] [IMAGE_NAME]\n");
		fprintf(stderr, "\te.g. ./ScaleRelations ads.pgm ads\n");
		exit(-1);
	}

	/*
	 * Miscellaneous steps that must be taken
	 */
	strcpy(imageFileName, argv[1]);
	strcpy(outFileStem, argv[2]);
	printf("Starting scale relations processing of image %s\n", imageFileName);

	getWidthAndHeight(imageFileName);
	findEdges(imageFileName, outFileStem);

	src.read(imageFileName);

	sprintf(str, "texts/%s_GD.txt", outFileStem);
	readGradientFile(str);
	sprintf(str, "texts/%s_MAG.txt", outFileStem);
	readMagnitudesFile(str);


	/*
	 * Process relational distributions at the different scales
	 */
	const int X_LOC = 210;
	const int Y_LOC = 150;
	int exponent = 1;
	int window = pow(2, exponent) + 1;

	while (window < 256) {
		REM group(1, window, window, 51, 100);
		int neighbor = (int) window / 2;

		for (int m = -neighbor; m <= neighbor; m++) {
			for (int n = -neighbor; n <= neighbor; n++) {
				int row = Y_LOC + m;
				int col = X_LOC + n;

				/* Check if neighborhood pixels are within image boundary */
				if (row > 0 && col > 0 && row < src.NR && col < src.NC) {
					int mn = (row * width) + col;
					group.addElement(row, col, gradDirections[mn], magnitudes[mn]);
				}
			}
		}
		group.buildHistogram();

		sprintf(str, "histograms/%s_H_%d.pgm", outFileStem, window);
		group.drawHistogram(str);

		sprintf(str, "scale_relations.txt");

		if ((textFile = fopen(str, "a")) == NULL) {
			printf("Cannot open text file %s\n", str);
			exit(-1);
		}

		fprintf(textFile, "Entropy for scale %d: %2.5f\n", window, group.calculateRenyiEntropy());

		exponent++;
		window = pow(2, exponent) + 1;
	}

	fclose(textFile);
	free(gradDirections);
	free(magnitudes);

	printf("Finished processing\n");

	return 0;
}

/*********************************** DEFINITIONS ****************************/



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


/**
 * Simply get the width and height of the image specified by 'filename'
 */
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

/**
 * Perform Canny edge detection on the image specified by 'filename'
 */
void findEdges(char imageFile[100], char outFileStem[100]) {
	char command[500];
	sprintf(command, "./CannyEdge %s 0.80 0.20 0.70 %s", imageFile, outFileStem);
	system(command);

}

/**
 * Read the gradient file specified by 'filename'
 */
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

/**
 * Read the magnitude file specified by 'filename'
 */
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


