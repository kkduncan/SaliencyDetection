/*
 * GetStructure.c
 *
 *  Created on: Jan 21, 2010
 *      Author: kesterduncan
 */
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Group.h"

/* Globals */
/*
 * Could probably define an array that stores the values of
 * each parameter
 */
int 	OPERATOR = 2; // canny
float 	scale = 0.1; // nice
float	low_thresh = 0.5; // nice
float	high_thresh = 0.5; // nice
float 	LENGTH_THRESH = 20.0; // good at 20
float 	SEG_INDEX = 400; // higher -> less junk segments
float 	MIN_PIXEL_DIST = 0.25; // 0.25
float 	MIN_SEG_LEN = 40; // 5, 15, 40


/** Output image */
int *outputImage;

/** Width and height of image */
int width, height;

/** Result from grouping */
int *PO;

/** Stores the gradient directions for pixels */
float *gradDirections;

/** Stores the magnitudes for each pixel */
int *magnitudes;

/** Group that corresponds to the whole image */
Group edgeImage;

/** Stores the entropy results */
FILE *resultsFile;


/* Prototypes */
void findEdges(char imageFile[100], char outFileStem[100]);
int readGradientFile(char filename[100]);
int readMagnitudesFile(char filename[100]);
void initializeGroup(int binDimension);
double calculatePerformance();
int getint(FILE *fp);

using namespace std;

int main(int argc, char *argv[]) {
	FILE *imgFile;
	char imageFileName[100];
	char outFileStem[100];
	char resultsFileName[100];
	char str[50];
	int binDimension;

	if (argc != 4) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		exit(-1);
	}

	strcpy(imageFileName, argv[1]);
	strcpy(outFileStem, argv[2]);
	binDimension = atoi(argv[3]);

	/* Get width & height of image */
	if ((imgFile = fopen(argv[1], "r")) == NULL) {
		printf("Cannot open image file %s\n", argv[1]);
		exit(-1);
	}

	fscanf(imgFile, "%s", str); // Get the P5 string
	width = getint(imgFile);
	height = getint(imgFile);

	if (width <= 0 || height <= 0) {
		printf("Image file %s is in an invalid format\n", argv[1]);
		exit(-1);
	}

	findEdges(imageFileName, outFileStem);

	sprintf(str, "texts/%s_GD.txt", outFileStem);
	readGradientFile(str);

	sprintf(str, "texts/%s_MAG.txt", outFileStem);
	readMagnitudesFile(str);

	initializeGroup(binDimension);

	sprintf(resultsFileName, "Results/%s_results/%s_results_whole.txt", outFileStem, outFileStem);
	if ((resultsFile = fopen(resultsFileName, "w")) == NULL) {
		fprintf(stderr, "Could not initialize results file\n");
		cerr << resultsFile;
		exit(-1);
	}

	fprintf(resultsFile, "The entropy value for the image %s  = %2.5f\n", imageFileName, edgeImage.calculateRenyiEntropy());
	fclose(resultsFile);

	sprintf(str, "histograms/%s_histogram.pgm", outFileStem);
	edgeImage.drawHistogram(str);
	sprintf(str, "texts/%s_histogram.txt", outFileStem);
	edgeImage.printHistogramMatrix(str);


	return 0;
}

void findEdges(char imageFile[100], char outFileStem[100]) {
//	int OPERATOR = 2; // canny
//	float scale = 0.1; // nice
//	float low_thresh = 0.5; // nice
//	float high_thresh = 0.5; // nice
//	float LENGTH_THRESH = 20.0; // good at 20
//	float SEG_INDEX = 400; // higher -> less junk segments
//	float MIN_PIXEL_DIST = 0.25; // 0.25
//	float MIN_SEG_LEN = 40; // 5, 15, 40
//	char temp[100];
	char command[500];
//
//	sprintf(temp, "%s.atts", outFileStem);
//	sprintf(command,
//			"./EdgeDetection %s %f -edge %d -lh %f -hh %f -lt %f -si %f -msl %f -md %f -seg 1 -o %s",
//			imageFile, scale, OPERATOR, low_thresh, high_thresh,
//			LENGTH_THRESH, SEG_INDEX, MIN_SEG_LEN, MIN_PIXEL_DIST,
//			outFileStem);
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

	fprintf(stderr, "Initialized gradient directions array\n");

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			float value;
			fscanf(gradFile, "%f ", &value);
			*(gradDirections + ij) = value;
		}
	}

	fprintf(stderr, "Read gradient directions from file\n");
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

	fprintf(stderr, "Initialized magnitudes array\n");

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			int value;
			fscanf(magFile, "%d ", &value);
			*(magnitudes + ij) = value;
		}
	}

	fprintf(stderr, "Read magnitudes directions from file\n");
	fclose(magFile);
	return 0;
}

void initializeGroup(int binDimension) {
	int i, j;

	edgeImage.init(1, width, height, binDimension);

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			int ij = (i * width) + j;
			edgeImage.addElement(i, j, gradDirections[ij], magnitudes[ij]);

		}
	}

	fprintf(stderr, "Finished adding %d elements\n", (int) edgeImage.members.size());

	edgeImage.buildSampledHistogram();
	fprintf(stderr, "Finished building histogram\n");
}


/**
 * Calculates the performance of the grouping as an entropy measure
 */
double calculatePerformance() {
	double measure;

	edgeImage.calculateRenyiEntropy();

	if(edgeImage.entropy != 0) {
		measure = 1 - edgeImage.entropy;
	} else {
		measure = 0;
	}

    return (measure);

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

