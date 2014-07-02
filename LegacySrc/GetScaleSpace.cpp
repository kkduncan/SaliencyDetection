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
#include <time.h>
#include <list>
#include "ScaleSpace.h"
#include "image.h"
#include <pthread.h>


/********************************* GLOBALS **********************************/
int width, height;
float *gradDirections;
float *maxScaleEntropies;
int *magnitudes;
int *scaleSpaceValues;
Image src;
Image scale1, scale2, scale3, scale4, scale5, scale6, scale7;



/******************************** PROTOTYPES *********************************/
void getWidthAndHeight(char *filename);
void findEdges(char imageFile[100], char outFileStem[100]);
int readGradientFile(char filename[100]);
int readMagnitudesFile(char filename[100]);
void eliminateLowMagnitudes();
void processPixelForEntropy(Image& tgt, int x, int y, char *textFilename);
void processScaleSpace(int x, int y);
void drawScales(char *filename);
void * threadProcess(void * params);
int getint(FILE *fp);

/********************************** MAIN *************************************/
using namespace std;

int main(int argc, char *argv[]) {
	char imageFileName[100];
	char outFileStem[100];
	char textFilename[100];
	char entropyFilename[100];
	char str[50];
	FILE *textFile;
	FILE *entropyFile;
	bool draw = false;

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
	if (atoi(argv[3]) == 1) draw = true;
	printf("Starting scale space processing of image %s\n", imageFileName);

	getWidthAndHeight(imageFileName);
	findEdges(imageFileName, outFileStem);

	src.read(imageFileName);
	tgt.initialize(src.NR, src.NC);

	/*
	 * Either perform the scale processing or draw the scales of a
	 * few select pixels.
	 */
	if (!draw) {
		sprintf(str, "texts/%s_GD.txt", outFileStem);
		readGradientFile(str);
		sprintf(str, "texts/%s_MAG.txt", outFileStem);
		readMagnitudesFile(str);

		/*
		 * Initialize scale space and max scale entropies arrays
		 */
		scaleSpaceValues = (int *) malloc(width * height * sizeof(int));
		maxScaleEntropies = (float *) malloc(width * height * sizeof(float));

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				scaleSpaceValues[(i * width) + j] = 0;
				maxScaleEntropies[(i * width) + j] = 0;
			}
		}

		/*
		 * Process Pixels in Parallel
		 */
		srand(time(NULL));
		int rowsPerThread = (int) ceil(height / NUM_THREADS);

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				processScaleSpace(i, j);
			}
		}

		/*
		 * Output information
		 */
		sprintf(str, "maps/%s_scale_map.pgm", outFileStem);
		tgt.save(str);

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				fprintf(textFile, "%d ", scaleSpaceValues[(i * width) + j]);
			}
			fprintf(textFile, "\n");
		}

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				fprintf(entropyFile, "%2.5f ", maxScaleEntropies[(i * width) + j]);
			}
			fprintf(entropyFile, "\n");
		}

		/*
		 * Clean up and close
		 */
		fclose(textFile);
		fclose(entropyFile);
		free(gradDirections);
		free(magnitudes);
		free(scaleSpaceValues);
		pthread_mutex_destroy(&writeMutex);

	} else {
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				tgt(i, j, RED) = src(i, j, RED);
				tgt(i, j, GREEN) = src(i, j, RED);
				tgt(i, j, BLUE) = src(i, j, RED);
			}
		}

		sprintf(textFilename, "texts/%s_scale_space.txt", outFileStem);
		srand(time(NULL));
		drawScales(textFilename);

		sprintf(str, "maps/%s_scales.ppm", outFileStem);
		tgt.save(str);
	}

	printf("Finished processing\n");
	pthread_exit(0);

	return 0;
}

/*********************************** DEFINITIONS ****************************/

void * threadProcess(void * params) {
	ThreadInfo *tInfo = (ThreadInfo *) params;

	for (int i = tInfo->rowStart; i < height && i < tInfo->rowEnd; i++) {
		for (int j = 0; j < width; j++) {
			processScaleSpace(i, j);
		}
	}
	return NULL;
}

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
		scaleSpace.scaleGroups[i].init(i + 1, scale, scale, 10, 5);
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

		if (row > 0 && row < tgt.NR && col > 0 && col < tgt.NC) {
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
	}

	/*
	 * Get the maximum scale. This is performed synchronously
	 */
	pthread_mutex_lock(&writeMutex);
	tgt(x, y) = (int) (((float)scaleSpace.getMaxScale() / scaleSpace.MAX_SCALE) * 255);
	scaleSpaceValues[(x * width) + y] = scaleSpace.maxScale;
	int index = ((int)(log(scaleSpace.maxScale - 1) / log(2))) - 1;
	maxScaleEntropies[(x * width) + y] = scaleSpace.scaleEntropies[index];
	pthread_mutex_unlock(&writeMutex);

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

void drawScales(char *filename) {
	FILE* scaleEntropyFile;
	int *scaleEntropies;
	int pixelsChosen = 0;
	int count = 20;

	if ((scaleEntropyFile = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing scale space values\n", filename);
		exit(-1);
	}

	scaleEntropies = (int *) malloc(width * height * sizeof(int));

	for (int k = 0; k < (width * height); k++) {
		*(scaleEntropies + k) = 0;
	}

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int ij = (i * width) + j;
			int value;
			fscanf(scaleEntropyFile, "%d ", &value);
			*(scaleEntropies + ij) = value;
		}
	}
	fclose(scaleEntropyFile);

	while (pixelsChosen < count) {
		int row = rand() % height;
		int col = rand() % width;
		int ij = (row * width) + col;
		int pixelScale = scaleEntropies[ij];

			/* Draw boundary of square neighborhood scale */
		if(tgt.inboundp(row + pixelScale, col + pixelScale) && tgt.inboundp(row + pixelScale, col - pixelScale) &&
		   tgt.inboundp(row - pixelScale, col + pixelScale) && tgt.inboundp(row - pixelScale, col + pixelScale)) {

			/* Highlight center pixel */
			tgt(row, col, RED) = 0;
			tgt(row, col, GREEN) = 198;
			tgt(row, col, BLUE) = 255;

			for (int i = -pixelScale; i <= pixelScale; i++) {
				tgt(row + i, col + pixelScale, RED) = 138;
				tgt(row + i, col + pixelScale, GREEN) = 255;
				tgt(row + i, col + pixelScale, BLUE) = 0;

				tgt(row + i, col - pixelScale, RED) = 138;
				tgt(row + i, col - pixelScale, GREEN) = 255;
				tgt(row + i, col - pixelScale, BLUE) = 0;

				tgt(row + pixelScale, col + i, RED) = 138;
				tgt(row + pixelScale, col + i, GREEN) = 255;
				tgt(row + pixelScale, col + i, BLUE) = 0;

				tgt(row - pixelScale, col + i, RED) = 138;
				tgt(row - pixelScale, col + i, GREEN) = 255;
				tgt(row - pixelScale, col + i, BLUE) = 0;

			}
			pixelsChosen++;
		}
	}
}

/**
 * Get the entropies for a pixel at different scales
 */
void processPixelForEntropy(Image& tgt, int x, int y, char *textFilename) {
	FILE *textFile;
	double previousEntropy = 0;
	int iteration = 0;
	int lastSignificantChange = 0;
	int stepsize = 16;
	bool converged = false;

	if ((textFile = fopen(textFilename, "w")) == NULL) {
		fprintf(stderr, "Cannot open text file %s for printing scale space values\n", textFilename);
		exit(-1);
	}

	int window;
	for (window = 7; window < (tgt.NC / 2)  && !converged; window += stepsize) {
		int neighbor = (int) window / 2;
		float entropy = 0;

		iteration++;
		printf("Iteration %d - window size: (%d x %d) \n", iteration, window, window);

		Group group;
		group.init(1, window, window, 10, 100);

		for (int m = -neighbor; m <= neighbor; m++) {
			for (int n = -neighbor; n <= neighbor; n++) {
				int row = x + m;
				int col = y + n;

				/* Check if neighborhood pixels are within image boundary */
				if (row > 0 && col > 0 && row < tgt.NR && col < tgt.NC) {
					int mn = (row * width) + col;
					group.addElement(row, col, gradDirections[mn], magnitudes[mn]);
				} else {
					group.addElement(x, y, M_PI, 0);
				}
			}
		}
		group.buildHistogram();
		entropy = group.calculateRenyiEntropy();

		if (fabs(entropy - previousEntropy) > 0.001) {
			lastSignificantChange = iteration;
		}
		if ((iteration - lastSignificantChange) > 5) {
			converged = true;
		}
		previousEntropy = entropy;

		/* Sanity check for unresolved bug */
		(entropy < 0) ?	entropy = 0.999 : entropy = entropy;
		fprintf(textFile, "%2.6f\n", entropy);
	}
	fprintf(textFile, "Last window size - %d x %d\n", window - stepsize, window - stepsize);
	fclose(textFile);
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


