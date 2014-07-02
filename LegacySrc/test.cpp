/*
 * test.cpp
 *
 *  Created on: Mar 8, 2010
 *      Author: kesterduncan
 */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>

using namespace std;

int main(int argc, char *argv[]) {
	float values[] = {0.27, 0.51, 0, 0, 0.81, 0.82, 0.71};
	int prevScale = 0;
	int currScale = 0;
	float first[7], second[7];

	for (int i=0; i < 7; i++) {
		float x_diff, y_diff;
		float n = pow(2, i);

		prevScale = currScale;
		currScale = (2 * n) + 1;

		if (i == 0) {
			y_diff = values[i];
		} else {
			y_diff = values[i] - values[i - 1];
		}

		x_diff = currScale - prevScale;
		first[i] = (float) (y_diff / x_diff);

		printf("First derivative: %2.5f \n", first[i]);


	}

	prevScale = 0;
	currScale = 0;

	for (int i=0; i < 7; i++) {
		float x_diff, y_diff;
		float n = pow(2, i);

		prevScale = currScale;
		currScale = (2 * n) + 1;

		if (i == 0) {
			y_diff = first[i];
		} else {
			y_diff = first[i] - first[i - 1];
		}

		x_diff = currScale - prevScale;
		second[i] = (float) y_diff / x_diff;

		printf("Second derivative: %2.5f \n", second[i]);
	}

	printf("Detecting strongest zero crossing\n");

	prevScale = 0;
	currScale = 0;
	float maxDiff = 0;
	int maxScale = 3;

	for (int i=0; i < 6; i++) {
		float y_diff = -999;
		float n = pow(2, i);

		prevScale = currScale;
		currScale = (2 * n) + 1;


		if ((second[i] > 0 && second[i + 1] < 0) || (second[i] < 0 && second[i + 1] > 0)) {
			y_diff = fabs(second[i] - second[i + 1]);
			printf("Diff = %2.5f\n", y_diff);
		}

		if (y_diff > maxDiff) {
			maxDiff = y_diff;
			maxScale = (2 * (pow(2, i + 1))) + 1;
		}

	}

	printf("Max scale: %d \n", maxScale);


	return 0;
}
