/*
 * GetDifferentScales.cpp
 *
 *  Created on: March 23, 2010
 *      Author: kesterduncan
 */

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <list>
#include "ScaleSpace.h"
#include "image.h"
#include <pthread.h>


/********************************* GLOBALS **********************************/
int width, height;
float *gradDirections;
int *magnitudes;
Image src;
Image scaleImages[7];


/******************************** PROTOTYPES *********************************/
void getWidthAndHeight(char *filename);
void findEdges(char imageFile[100], char outFileStem[100]);
int readGradientFile(char filename[100]);
int readMagnitudesFile(char filename[100]);
void eliminateLowMagnitudes();
void processScaleSpace(int x, int y);
int getint(FILE *fp);

/********************************** MAIN *************************************/
using namespace std;

int main(int argc, char *argv[]) {
	char imageFileName[100];
	char outFileStem[100];
	char str[50];


	/*
	 * Check command line arguments
	 */
	if (argc != 3) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		fprintf(stderr, "\tUSAGE: ./ScaleSpace [IMAGE] [IMAGE_NAME]\n");
		fprintf(stderr, "\te.g. ./ScaleSpace ads.pgm ads\n");
		exit(-1);
	}

	/*
	 * Miscellaneous steps that must be taken
	 */
	strcpy(imageFileName, argv[1]);
	strcpy(outFileStem, argv[2]);
	printf("Starting scale space processing of image %s\n", imageFileName);

	getWidthAndHeight(imageFileName);
	findEdges(imageFileName, outFileStem);

	src.read(imageFileName);
	for (int i = 0; i < 7; i++) {
		scaleImages[i].initialize(src.NR, src.NC);
	}


	sprintf(str, "texts/%s_GD.txt", outFileStem);
	readGradientFile(str);
	sprintf(str, "texts/%s_MAG.txt", outFileStem);
	readMagnitudesFile(str);


	srand(time(NULL));
	eliminateLowMagnitudes();
	for (int i = 0; i < src.NR; i++) {
		for (int j = 0; j < src.NC; j++) {
			processScaleSpace(i, j);
		}
	}

	/*
	 * Output information
	 */
	for (int i = 0; i < 7; i++) {
		sprintf(str, "maps/%s_scale_%d_map.pgm", outFileStem, i);
		scaleImages[i].save(str);
	}

	/*
	 * Clean up and close
	 */
	free(gradDirections);
	free(magnitudes);
	printf("Finished processing\n");

	return 0;
}

/*********************************** DEFINITIONS ****************************/

/**
 * Processes the scale space with a greedy algorithm
 */
void processScaleSpace(int x, int y) {
	ScaleSpace scaleSpace;
	int pixelsChecked = 0;
	const int OFFSET = scaleSpace.MAX_SCALE;
	const int PERCENTAGE = (int) ceil(0.05 * (pow(scaleSpace.MAX_SCALE, 2)));
	double distance, angle, magnitude;
	list<GroupElement>::iterator xElem, yElem;

	/* Initialize groups for each scale */
	for (int i = 0; i < scaleSpace.NUM_SCALES; i++) {
		int n = i^2;
		int scale = (2 * n) + 1;
		scaleSpace.scaleGroups[i].init(i + 1, scale, scale, 10);
	}

	/*
	 * Build up the pixel group for the max scale by randomly
	 * picking pixel points up to a certain percentage
	 * of pixels in the image
	 * All the other scales are incorporated in this max scale
	 */
	while (pixelsChecked != PERCENTAGE) {
		int row = x + (rand() % (OFFSET * 2)) - 128;
		int col = y + (rand() % (OFFSET * 2)) - 128;
		int mn = (row * width) + col;

		if (row > 0 && row < scaleImages[0].NR && col > 0 && col < scaleImages[0].NC) {
			scaleSpace.scaleGroups[6].addElement(row, col, gradDirections[mn], magnitudes[mn]);

		} else {
			/* This compensates for pixels that are not within the
			 * image bounds. Doesn't contribute anything basically.
			 */
			scaleSpace.scaleGroups[6].addElement(x, y, M_PI, 0);
		}
		pixelsChecked++;
	}

	/*
	 * Build up the scale histograms
	 */
	for (xElem = scaleSpace.scaleGroups[6].members.begin(); xElem != scaleSpace.scaleGroups[6].members.end(); xElem++) {
		for (yElem = scaleSpace.scaleGroups[6].members.begin(); yElem != scaleSpace.scaleGroups[6].members.end(); yElem++) {
			distance = abs(xElem->x - yElem->x) + abs(xElem->y - yElem->y);
			angle = fmod(((2 * M_PI) + (xElem->angle - yElem->angle)), (2 * M_PI));
			if (angle > M_PI) {
				angle = (2 * M_PI) - angle;
			}
			magnitude = fabs(xElem->magnitude * yElem->magnitude);

			for (int i = 0; i < scaleSpace.NUM_SCALES; i++) {
				int n = 2^i;
				int scale = (2 * n) + 1;
				int x1_offset = abs(xElem->x - x);
				int y1_offset = abs(xElem->y - y);
				int x2_offset = abs(xElem->x - x);
				int y2_offset = abs(xElem->y - y);

				if (((x1_offset / scale) <= 1) && ((y1_offset / scale) <= 1) &&
					((x2_offset / scale) <= 1) && ((y2_offset / scale) <= 1)){
					scaleSpace.scaleGroups[i].updateHistogram(distance, angle, magnitude);
				}
			}
		}
	}

	/*
	 * Find entropies for each scale
	 */
	for (int i = 0; i < scaleSpace.NUM_SCALES; i++) {
		scaleSpace.scaleEntropies[i] = scaleSpace.scaleGroups[i].calculateShannonEntropy();
		scaleImages[i](x, y) = scaleSpace.scaleEntropies[i];
	}
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


