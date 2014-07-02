/*
 * TestRImage.cpp
 *
 *  Created on: Feb 16, 2011
 *      Author: kesterduncan
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RImage.h"
#include "IntegralImage.h"

using namespace std;

int main (int argc, char *argv[]) {
	RImage src;
	RImage newImg;
	IntegralImage intImg;
	float *means;
	float *stddev;
	int window = 5;
	int neighbor = (int) (window / 2);
	float w2 = window * window;

	src.read("1.pgm");
	newImg.initialize(src.numRows, src.numCols);

	RImage iImg = intImg.getIntegralImage(src);

	means = (float*) malloc (sizeof(float) * src.numCols * src.numRows);
	stddev = (float*) malloc (sizeof(float) * src.numCols * src.numRows);

	for (int i = 0; i < src.numCols * src.numRows; i++) {
		means[i] = 0;
		stddev[i] = 0;
	}

	/* Calculating the means */
	for (int x = 0; x < iImg.numCols; x++) {
		for (int y = 0; y < iImg.numRows; y++) {
			int ij = y * iImg.numCols + x;

			int prevY = y - neighbor;
			int prevX = x - neighbor;
			int nextY = y + neighbor;
			int nextX = x + neighbor;

			float firstValue = 0;
			float secondValue = 0;
			float thirdValue = 0;
			float fourthValue = 0;

			if (nextY < iImg.numRows && nextX < iImg.numCols) {
				firstValue = iImg(nextY, nextX, RED);
			}

			if (prevY > 0 && prevX > 0) {
				secondValue = iImg(prevY, prevX, RED);
			}

			if (nextX < iImg.numCols && prevY > 0) {
				thirdValue = iImg(prevY, nextX, RED);
			}

			if (prevX > 0 && nextY < iImg.numRows) {
				fourthValue = iImg(nextY, prevX, RED);
			}

			if (prevX >= 0 && prevY >= 0 && nextY < iImg.numRows && nextX < iImg.numCols) {
				means[ij] = ((firstValue + secondValue) - (thirdValue + fourthValue)) / w2;
			}
		}
	}

	/* Calculating SD */
	for (int j = 0; j < src.numCols; j++) {
		for (int i = 0; i < src.numRows; i++) {
			int ij = i * src.numCols + j;
			float sum = 0;

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int y = i + m;
					int x = j + n;

					if (y > 0 && x > 0 && y < src.numRows && x < src.numCols) {
						sum += (pow(src(y, x, RED) - means[y * src.numCols + x], 2));
					}
				}
			}

			stddev[ij] = sqrt(sum / w2);
		}
	}

	/* Threshold */
	for (int j = 0; j < src.numCols; j++) {
		for (int i = 0; i < src.numRows; i++) {
			int ij = i * src.numCols + j;
			double threshold = means[ij] * (1 + 0.3 * ((stddev[ij] / 127)));

			if (src(i, j, RED) <= threshold) {
				newImg(i, j, RED) = 255;
			} else {
				newImg(i, j, RED) = 0;
			}
		}
	}


	newImg.save("thresholded_image.pgm");

	newImg.flushMemory();
	src.flushMemory();
	iImg.flushMemory();

	free(means);
	free(stddev);

	return 0;
}

