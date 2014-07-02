/*
 * IntegralImage.cpp
 *
 *  Created on: Apr 9, 2011
 *      Author: kesterduncan
 */

#include "IntegralImage.h"

IntegralImage::IntegralImage() {

}

IntegralImage::~IntegralImage() {

}

void IntegralImage::getIntegralImage(const RImage& orig, RImage &intImg) {
	intImg = orig;

	int *rowSum;
	rowSum = new int [orig.numRows * orig.numCols];
	if (rowSum == NULL) {
		fprintf(stderr, "INT_IMG ERROR: Insufficient space for integral image.\n");
		exit(1);
	}

	for (int i = 0; i < orig.numRows; i++) {
		for (int j = 0; j < orig.numCols; j++) {
			int prevRow = i - 1;
			int prevCol = j - 1;
			double value;

			if (orig.isColorImage == 1) {
				value = (orig(i, j, RED) + orig(i, j, GREEN) + orig(i, j, BLUE)) / 3;
			} else {
				value = orig(i, j, RED);
			}

			if (prevRow == -1) {
				rowSum[i * orig.numCols + j] = 0;
			} else {
				rowSum[i * orig.numCols + j] = rowSum[(i - 1) * orig.numCols + j] + value;
			}

			if (prevCol == -1) {
				intImg(i, j, RED) = 0;
			} else {
				intImg(i, j, RED) = intImg(i, prevCol, RED) + rowSum[i * orig.numCols + j];
			}
		}
	}

	delete[](rowSum);
}
