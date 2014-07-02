/*
 * K. Duncan
 * University of South Florida
 * Copyright (c) 2010
 *
 * -- GaussianPyramid.cpp --
 *
 * Implementation file for the Gaussian Pyramid technique outlined in the
 * work: "The Laplacian Pyramid as a Compact Image Code", Burt P.,
 * Adelson E., IEEE Transactions on Communications, April 1983
 *
 */

#include "GaussianPyramid.h"

/*
 * Constructors / Destructor
 */
GaussianPyramid::GaussianPyramid() {

}


GaussianPyramid::GaussianPyramid(const RImage &src) {
	setBaseLevel(src);

}


GaussianPyramid::~GaussianPyramid() {
	for (int i = 0; i < NUM_LEVELS; i++) {
		levels[i].flushMemory();
	}
}


/**
 * Sets level 0 of the pyramid, determines the dimensions of
 * the subsequent levels, and initializes them.
 */
void GaussianPyramid::setBaseLevel(const RImage& img) {
	/*
	 * Determine if the image can be processed
	 */
	if ((img.numRows / pow(2, NUM_LEVELS)) > DIM_LIMIT &&
			(img.numCols / pow(2, NUM_LEVELS)) > DIM_LIMIT) {

		for (int i = 0; i < NUM_LEVELS; i++) {
			int rows = (int) ceil (img.numRows / pow(2, i));
			int cols = (int) ceil (img.numCols / pow(2, i));

			dimensions[i][0] = rows;
			dimensions[i][1] = cols;
			levels[i].initialize(rows, cols);
			if (img.isColorImage) {
				levels[i].convertToColor();
			}
		}

		/*
		 * Set the base level
		 */
		levels[0] = img;

	} else {
		fprintf(stderr, "PYRAMID ERROR: The image cannot be processed. Too small.\n");
		exit(1);
	}

}


/**
 * Reduces the image given by srcImage to the pyramid level specified by
 * tgtLevel
 */
void GaussianPyramid::reduceImageToLevel(const int srcLevel,
										 const int tgtLevel) {
	/*
	 * Assumes that srcLevel and tgtLevel are valid
	 */
	if (tgtLevel > srcLevel && levels[srcLevel].numCols > DIM_LIMIT
			&& levels[srcLevel].numCols > DIM_LIMIT) {
		int weights[][5] = {{1, 4, 6, 4, 1},
				{4, 16, 24, 16, 4},
				{6, 24, 36, 24, 6},
				{4, 16, 24, 16, 4},
				{1, 4, 6, 4, 1}
		};

		int levelDiff = pow(2, (tgtLevel - srcLevel) - 1);

		/*
		 * These dimensions have been pre-calculated in the setBaseLevel function
		 */
		int tgtRows = dimensions[tgtLevel][0];
		int tgtCols = dimensions[tgtLevel][1];

		levels[tgtLevel].type = levels[srcLevel].type;
		levels[tgtLevel].isColorImage = levels[srcLevel].isColorImage;

		for (int i = 0; i < tgtRows; i++) {
			for (int j = 0; j < tgtCols; j++) {
				double rValue = 0;
				double gValue = 0;
				double bValue = 0;

				int rowDiff = (pow(2, levelDiff) * i);
				int colDiff = (pow(2, levelDiff) * j);

				for (int m = -2; m <= 2; m++) {
					for (int n = -2; n <= 2; n++) {
						int row = rowDiff + m;
						int col = colDiff + n;

						if (levels[srcLevel].inBounds(row, col)) {
							int temp1 = levels[srcLevel](row, col, RED);
							rValue += weights[m + 2][n + 2] * temp1;

							if (levels[tgtLevel].isColorImage == 1) {
								int temp2 = levels[srcLevel](row, col, GREEN);
								int temp3 = levels[srcLevel](row, col, BLUE);

								gValue += weights[m + 2][n + 2] * temp2;
								bValue += weights[m + 2][n + 2] * temp3;
							}
						}
					}
				}

				rValue = rValue / 256.0;
				levels[tgtLevel](i, j, RED) = (int) (rValue);

				if (levels[tgtLevel].isColorImage == 1) {
					gValue = gValue / 256.0;
					levels[tgtLevel](i, j, GREEN) = (int) (gValue);
					bValue = bValue / 256.0;
					levels[tgtLevel](i, j, BLUE) = (int) (bValue);
				}
			}
		}

	} else {
		fprintf(stderr, "PYRAMID ERROR: Image cannot be reduced.\n");
		exit(1);
	}
}


/**
 * Expands the image given by srcImage to the pyramid level specified by
 * tgtLevel
 */
void GaussianPyramid::expandImageToLevel(const int srcLevel,
										 const int tgtLevel) {

	/*
	 * Assusmes that srcLevel and tgtLevel are valid
	 */
	if (tgtLevel < srcLevel && levels[srcLevel].numCols > DIM_LIMIT
			&& levels[srcLevel].numRows > DIM_LIMIT) {
		int weights[][5] = {{1, 4, 6, 4, 1},
				{4, 16, 24, 16, 4},
				{6, 24, 36, 24, 6},
				{4, 16, 24, 16, 4},
				{1, 4, 6, 4, 1}
		};

		int levelDiff = pow(2, (srcLevel - tgtLevel) - 1);
		int tgtRows = dimensions[tgtLevel][0];
		int tgtCols = dimensions[tgtLevel][1];

		RImage tgtImg;
		tgtImg.type = levels[srcLevel].type;
		tgtImg.isColorImage = levels[srcLevel].isColorImage;
		tgtImg.numRows = tgtRows;
		tgtImg.numCols = tgtCols;

		tgtImg.initialize(tgtRows, tgtCols);

		for (int i = 0; i < tgtRows; i++) {
			for (int j = 0; j < tgtCols; j++) {
				double rValue = 0;
				double gValue = 0;
				double bValue = 0;
				int diff = 2 * levelDiff;

				for (int m = -2; m <= 2; m++) {
					for (int n = -2; n <= 2; n++) {
						int row = (i - m) / diff;
						int col = (j - n) / diff;

						if (row < 0) row = 0;
						if (col < 0) col = 0;

						/*
						 * Only terms for which the row and col values are integers
						 * are considered in this sum
						 */
						if (ceil(row) == row && ceil(col) == col) {
							int temp1 = levels[srcLevel]((int)row, (int)col, RED);
							rValue += weights[m + 2][n + 2] * temp1;

							if (tgtImg.isColorImage == 1) {
								int temp2 = levels[srcLevel]((int)row, (int)col, GREEN);
								gValue += weights[m + 2][n + 2] * temp2;

								int temp3 = levels[srcLevel]((int)row, (int)col, BLUE);
								bValue += weights[m + 2][n + 2] * temp3;
							}
						}
					}
				}

				rValue = rValue / 256.0;
				tgtImg(i, j, RED) = (int) rValue;

				if (tgtImg.isColorImage) {
					gValue = gValue / 256.0;
					tgtImg(i, j, GREEN) = (int) gValue;

					bValue = bValue / 256.0;
					tgtImg(i, j, BLUE) = (int) bValue;
				}

			}
		}

		/*
		 * Change the src image to point to the newly
		 * created expanded image.
		 */
		levels[srcLevel].flushMemory();
		levels[srcLevel] = tgtImg;

	} else {
		fprintf(stderr, "PYRAMID ERROR: Image cannot be reduced.\n");
		exit(1);
	}
}


/**
 * Reduce the base image to created the pyramid
 * structure
 */
void GaussianPyramid::reduceAll() {
	/*
	 * This assumes that the pyramid is not null and that the
	 * base level of the pyramid is set.
	 */
	for (int i = 1; i < NUM_LEVELS; i++) {
		reduceImageToLevel(i - 1, i);
	}
}


/**
 * Expand all the images in the pyramid to the original
 * source dimensions
 */
void GaussianPyramid::expandAll() {
	for (int i = 1; i < NUM_LEVELS; i++) {
		expandImageToLevel(i, 0);
	}
}


/**
 * Performs pyramid processing. The base level (level 0) must be set
 * before using this function.
 */
void GaussianPyramid::processPyramid() {

	/*
	 * This assumes that the pyramid is not null and that the
	 * base level of the pyramid is set.
	 */
	for (int i = 1; i < NUM_LEVELS; i++) {
		reduceImageToLevel(i - 1, i);
	}

	for (int i = 1; i < NUM_LEVELS; i++) {
		expandImageToLevel(i, 0);
	}

}


