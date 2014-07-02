/*
 * ScaleSpace.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: kesterduncan
 */

#include "ScaleSpace.h"

ScaleSpace::ScaleSpace() {
	// TODO Auto-generated constructor stub

}

ScaleSpace::~ScaleSpace() {
	// TODO Auto-generated destructor stub
}

void ScaleSpace::getFirstDerivative() {
	int prevScale = 0;
	int currScale = 0;

	for (int i=0; i < 7; i++) {
		float x_diff, y_diff;
		int n = (int) pow(2, i);

		prevScale = currScale;
		currScale = (2 * n) + 1;

		if (i == 0) {
			y_diff = scaleEntropies[i];
		} else {
			y_diff = scaleEntropies[i] - scaleEntropies[i - 1];
		}

		x_diff = currScale - prevScale;
		firstDerivative[i] = (float) y_diff / x_diff;
	}
}

void ScaleSpace::getSecondDerivative() {
	int prevScale = 0;
	int currScale = 0;

	for (int i = 0; i < 7; i++) {
		float x_diff, y_diff;
		int n = (int) pow(2, i);

		prevScale = currScale;
		currScale = (2 * n) + 1;

		if (i == 0) {
			y_diff = firstDerivative[i];
		} else {
			y_diff = firstDerivative[i] - firstDerivative[i - 1];
		}

		x_diff = currScale - prevScale;
		secondDerivative[i] = (float) y_diff / x_diff;

	}
}

int ScaleSpace::getMaxScale() {
	float maxDiff = -999;
	int prevScale = 0;
	int currScale = 0;

	maxScale = 3;

	getFirstDerivative();
	getSecondDerivative();

	for (int i = 0; i < 6; i++) {
		float y_diff;
		int n = (int) pow(2, i);

		prevScale = currScale;
		currScale = (2 * n) + 1;

		if ((secondDerivative[i] > 0 && secondDerivative[i + 1] < 0) ||
			(secondDerivative[i] < 0 && secondDerivative[i + 1] > 0)) {
			y_diff = fabs(secondDerivative[i] - secondDerivative[i + 1]);
		}

		if (y_diff > maxDiff) {
			maxDiff = y_diff;
			maxScale = (int)(2 * (pow(2, i + 1))) + 1;
		}
	}

	return maxScale;
}
