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
float 	LENGTH_THRESH = 3.0; // good at 20
float 	SEG_INDEX = 200; // higher -> less junk segments
float 	MIN_PIXEL_DIST = 0.15; // 0.25
float 	MIN_SEG_LEN = 40; // 5, 15, 40


/**
 * Output image
 */
int *outputImage;

/**
 * Width and height of image
 */
int width, height;

/**
 * Stores the segment ids
 */
int *PO;

/**
 * Stores the gradient directions for pixels
 */
float *gradDirections;

Group edgeImage;

FILE *resultsFile;


/* Prototypes */
void findEdges(char imageFile[100], char outFileStem[100]);
int readAttributesFile(char filename[100]);
void initializeGroup(int dimension);
double calculatePerformance();
int getint(FILE *fp);


int main(int argc, char *argv[]) {

	char imageFileName[100];
	char outFileStem[100];
	char resultsFileName[100];
	char str[50];
	int dimension = 0;
	FILE *imgFile;

	if (argc != 4) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		exit(-1);
	}

	strcpy(imageFileName, argv[1]);
	strcpy(outFileStem, argv[2]);
	dimension = atoi(argv[3]);

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
	sprintf(str, "%s.atts", outFileStem);
	readAttributesFile(str);

	initializeGroup(dimension);

	sprintf(resultsFileName, "%s_results_whole_%d.txt", outFileStem, dimension);
	if ((resultsFile = fopen(resultsFileName, "w")) == NULL) {
		fprintf(stderr, "Could not initialize results file\n");
	}
	fprintf(resultsFile, "The Renyi entropy value for the image %s  = %2.5f\n", imageFileName, edgeImage.calculateRenyiEntropy());
	fprintf(resultsFile, "The Shannon entropy value for the image %s  = %2.5f\n", imageFileName, edgeImage.calculateShannonEntropy());
	fprintf(stderr, "The entropy value for the image %s  = %2.5f\n", imageFileName, edgeImage.entropy);
	fclose(resultsFile);

	sprintf(resultsFileName, "histograms/%s_histogram.pgm", outFileStem);
	edgeImage.drawHistogram(resultsFileName);

	sprintf(resultsFileName, "texts/%s_histogram.txt", outFileStem);
	edgeImage.printHistogramMatrix(resultsFileName);

	return 0;
}

void findEdges(char imageFile[100], char outFileStem[100]) {
	char command[500];
	char temp[100];

	sprintf(temp, "%s.atts", outFileStem);

	sprintf(command,
			"./EdgeDetection %s %f -edge %d -lh %f -hh %f -lt %f -si %f -msl %f -md %f -seg 1 -o %s",
			imageFile, scale, OPERATOR, low_thresh, high_thresh,
			LENGTH_THRESH, SEG_INDEX, MIN_SEG_LEN, MIN_PIXEL_DIST,
			outFileStem);

	if (system(command)) printf("Edge detection successfully done\n");

}

int readAttributesFile(char filename[100]) {
	FILE *infile;
	char paren;
	int flag, flag1, x, y, i, len;
	float dummy;
	int dummy1, id;
	float gradDir;

	PO = (int *) malloc(width * height * sizeof(int));
	outputImage = (int *) malloc(width * height * sizeof(int));
	gradDirections = (float *) malloc(width * height * sizeof(float));

	// Initialize arrays
	for (i = 0; i < width * height; i++) {
		*(PO + i) = 0;
		*(outputImage + i) = 0;
		*(gradDirections + i) = 0.0;
	}

	if ((infile = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open file, %s for reading.\n", filename);
		return (0);
	}

	flag = 0;
	fscanf(infile, " %c", &paren);
	if ((int) paren == 40) {
		while (flag == 0) {
			fscanf(infile, " %c", &paren);
			if ((int) paren == 41) {
				flag = 1;
			} else {
				fscanf(infile, " %d", &(id));
				fscanf(infile, " %c", &paren);/** reads in ( */
				fscanf(infile, " %f", &(dummy)); // x center
				fscanf(infile, " %f", &(dummy)); // y center
				fscanf(infile, " %f", &(dummy)); // radius
				fscanf(infile, " %f", &(dummy)); // fit_error
				fscanf(infile, " %c", &paren);/** reads in ( */
				fscanf(infile, " %d", &dummy1); // start_pt
				fscanf(infile, " %d", &dummy1);
				fscanf(infile, " %c", &paren);/** reads in ) */
				fscanf(infile, " %c", &paren);/** reads in ( */
				fscanf(infile, " %d", &dummy1); // end_pt
				fscanf(infile, " %d", &dummy1);
				fscanf(infile, " %c", &paren);/** reads in ) */
				fscanf(infile, " %f", &(dummy)); // gray_slope
				fscanf(infile, " %f", &dummy); // dummy value
				fscanf(infile, " %f", &(dummy)); // mag+
				fscanf(infile, " %f", &(dummy)); // width+
				fscanf(infile, " %f", &(dummy)); // mag-
				fscanf(infile, " %f", &(dummy)); // width-

				fscanf(infile, " %f", &(gradDir)); /* reads in the grad dir */

				fscanf(infile, " %f", &dummy); //curvature
				fscanf(infile, " %c", &paren);/** reads in ( */
				fscanf(infile, " %f", &dummy); // centroid
				fscanf(infile, " %f", &dummy);
				fscanf(infile, " %c", &paren);/** reads in ) */
				fscanf(infile, " %f", &dummy); // start_slope
				fscanf(infile, " %f", &dummy); // end_slope
				fscanf(infile, " %c", &paren);/** reads in ) */
				fscanf(infile, " %c", &paren);/** reads in ( */
				flag1 = 0;
				len = 0;
				while (flag1 == 0) {
					fscanf(infile, " %c", &paren);
					if ((int) paren == 41) {
						flag1 = 1;
					} else {
						fscanf(infile, "%d ", &y);
						fscanf(infile, "%d ", &x);

						/* Contiguous segments might share end points */
						if (*(outputImage + (y) * width + x) == 0) {
							len++;

							/* Pixels that belong to the segment indicated by 'id'
							 * are labelled with 'id' in PO
							 *
							 * PO assigned id+1
							 */
							*(PO + ((y * width) + x)) = id + 1;
							*(outputImage + ((y * width) + x)) = id;
							*(gradDirections + ((y * width) + x)) = gradDir;

						}

						fscanf(infile, " %c", &paren);/** reads in ) */
					}
				}

				fscanf(infile, " %c", &paren);/** reads in ) */
			}
		}
	}

	fclose(infile);

	return 1;
}

void initializeGroup(int dimension) {
	int i, j;

	edgeImage.init(1, width, height, dimension);

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			int ij = (i * width) + j;

			if (PO == NULL) {
				printf("PO is NULL\n");
				exit(-1);
			}
			if (PO[ij] != 0) {
				edgeImage.addElement(i, j, *(gradDirections + ((i * width) + j)));
			}
		}
	}

	edgeImage.buildHistogram();
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

