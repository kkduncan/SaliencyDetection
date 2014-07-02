/*
 * GetWindowStructure.cpp
 *
 *  Created on: Jan 23, 2010
 *      Author: kesterduncan
 */

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <list>
#include "Group.h"
#include "image.h"


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


/** Items related to group windows */
const int NUM_COLS = 2;
const int NUM_ROWS = 3;
const int WINDOWS = NUM_COLS * NUM_ROWS;
int columnWidth = 0;
int rowHeight = 0;
Group groupWindows[NUM_ROWS][NUM_COLS];
list<Group> groupList;
Image src, tgt;
const int COLORS[5][3] = { {7, 234, 0}, {0, 187, 176}, {115, 112, 0}, {156, 66, 0}, {202, 14, 0}};
FILE *resultsFile;


/* Prototypes */
void getWindowDimensions();
void findEdges(char imageFile[100], char outFileStem[100]);
int readAttributesFile(char filename[100]);
void initializeGroups();
bool compareGroups(Group a, Group b);
void convertGroupArrayToList(Group groups[NUM_ROWS][NUM_COLS]);
void drawTopFiveGroups (char *filestem);
int getint(FILE *fp);

using namespace std;

int main(int argc, char *argv[]) {
	char imageFileName[100];
	char outFileStem[100];
	char resultsFileName[100];
	char str[50];
	FILE *imgFile;

	if (argc != 3) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		exit(-1);
	}

	strcpy(imageFileName, argv[1]);
	strcpy(outFileStem, argv[2]);

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

	getWindowDimensions();
	findEdges(imageFileName, outFileStem);

	sprintf(str, "%s.atts", outFileStem);
	readAttributesFile(str);
	initializeGroups();

	sprintf(resultsFileName, "Results/%s_results/%s_results_%d.txt", outFileStem, outFileStem, WINDOWS);
	if ((resultsFile = fopen(resultsFileName, "w")) == NULL) {
		fprintf(stderr, "Could not initialize results file\n");
	}

	/* Print entropy values to file */
	fprintf(resultsFile, "---- RESULTS FOR IMAGE FILE %s USING %d (rows) x %d (cols) WINDOWS ----\n", imageFileName, NUM_ROWS, NUM_COLS);
	for (int i = 0; i < NUM_ROWS; i++) {
		for (int j = 0; j < NUM_COLS; j++) {
			fprintf(resultsFile, "WINDOW (%d, %d) - Members = %d, Entropy= %2.4f\n", i, j, groupWindows[i][j].elemCount, groupWindows[i][j].calculateRenyiEntropy());
		}
	}
	fclose(resultsFile);

	/* Drawing */
	src.read(imageFileName);
	tgt.initialize(src.NR, src.NC);
	convertGroupArrayToList(groupWindows);
	drawTopFiveGroups(outFileStem);

	return 0;
}

void getWindowDimensions() {
	/* Revise, can be very buggy */
	int extraWidth = 0;
	int extraHeight = 0;

	if ((width % NUM_COLS) != 0) {
		extraWidth = NUM_COLS - (width % NUM_COLS);
	}
	if ((height % NUM_ROWS) != 0) {
		extraHeight = NUM_ROWS - (height % NUM_ROWS);
	}

	columnWidth = (width + extraWidth) / NUM_COLS; // Add more pels to eliminate fractions and overcompensate
	rowHeight = (height + extraHeight) / NUM_ROWS;

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

//	fprintf(stderr, "%s\n", command);
	system(command);
	sprintf(command, "mv %s.edge %s_edge.pgm", outFileStem, outFileStem);
//	fprintf(stderr, "%s\n", command);
	system(command);

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

void initializeGroups() {
	int i, j, index = 0;

	for (i = 0; i < NUM_ROWS; i++) {
		for (j = 0; j < NUM_COLS; j++) {
			groupWindows[i][j].init(index, columnWidth, rowHeight);
			index++;

		}
	}

	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			int ij = (i * width) + j;
			int rowIndex = int ((i + 1) / rowHeight);
			int colIndex = int ((j + 1) / columnWidth);

			if (PO == NULL) {
				printf("PO is NULL\n");
				exit(-1);
			}
			if (PO[ij] != 0) {
				groupWindows[rowIndex][colIndex].addElement(i, j, *(gradDirections + ((i * width) + j)));
			}
		}
	}

	for (i = 0; i < NUM_ROWS; i++) {
		for (j = 0; j < NUM_COLS; j++) {
			groupWindows[i][j].buildHistogram();
		}
	}
}

bool compareGroups(Group a, Group b) {
	return (a.entropy < b.entropy) ? true : false;
}

/**
 * Converts groups to a list and sorts it in ascending order
 * according to the entropy of each group
 * @param groups
 */
void convertGroupArrayToList(Group groups[NUM_ROWS][NUM_COLS]) {
	for (int i = 0; i < NUM_ROWS; i++) {
		for (int j = 0; j < NUM_COLS; j++) {
			if(groups[i][j].elemCount > 20) {
				groupList.push_back(groups[i][j]);
			}
		}
	}

	groupList.sort(compareGroups);
}

void drawTopFiveGroups (char *filestem) {

	list<Group>::iterator groupIt;
	list<GroupElement>::iterator elemIt;
	char outputName[100];

	for (int i = 0; i < src.NR; i++) {
		for (int j = 0; j < src.NC; j++) {
			tgt(i, j, RED) = src(i, j, RED);
			tgt(i, j, GREEN) = src(i, j, RED);
			tgt(i, j, BLUE) = src(i, j, RED);
		}
	}

	int k = 0;
	for (groupIt = groupList.begin(); groupIt != groupList.end() && k < 5; groupIt++, k++) {
		for(elemIt = groupIt->members.begin(); elemIt != groupIt->members.end(); elemIt++) {
			tgt(elemIt->x, elemIt->y, RED) = src(elemIt->x, elemIt->y, RED) + COLORS[k][0];
			tgt(elemIt->x, elemIt->y, GREEN) = src(elemIt->x, elemIt->y, RED) + COLORS[k][1];
			tgt(elemIt->x, elemIt->y, BLUE) = src(elemIt->x, elemIt->y, RED) + COLORS[k][2];

		}
	}

	sprintf(outputName, "Results/%s_results/%s_top_groups_%dx%d.ppm", filestem, filestem, NUM_ROWS, NUM_COLS);
	tgt.save(outputName);
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

