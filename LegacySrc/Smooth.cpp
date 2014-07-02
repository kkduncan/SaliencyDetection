/*
 * Smooth.cpp
 *
 *  Created on: Apr 1, 2010
 *      Author: kesterduncan
 */

#include <iostream>
#include <cstdio>
#include "image.h"

using namespace std;

int main(int argc, char *argv[]) {
	Image src, tgt;
	char srcName[100], tgtName[100];
	int gaussian_filter[5][5] =
			{{1, 4,  7,  4,  1},
			{ 4, 16, 26, 16, 4},
			{ 7, 26, 41, 26, 7},
			{ 4, 16, 26, 16, 4},
			{ 1, 4,  7,  4,  1}
	};

	float maxi = 0, mini = 255;
	float threshold;

	strcpy(srcName, argv[1]);
	src.read(srcName);
	tgt.initialize(src.NR, src.NC);

	for (int i=0; i < src.NR; i++) {
		for (int j = 0; j < src.NC; j++) {
			int neighbor = 2;
			int sum = 0;
			int count = 273;
			float avg = 0;

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int row = i + m;
					int col = j + n;

					/* Check if neighborhood pixels are within image boundary */
					if (row > 0 && col > 0 && row < src.NR && col < src.NC) {
						sum += (src(row, col) * gaussian_filter[m + neighbor][n + neighbor]);

					}
				}
			}

			avg = float(sum / count);
			if (avg > maxi) maxi = avg;
			if (avg < mini) mini = avg;

			tgt(i, j) = (int) avg;
		}
	}

	threshold = (maxi + mini) / 2;

	for (int i=0; i < src.NR; i++) {
		for (int j = 0; j < src.NC; j++) {
			int neighbor = 2;
			int sum = 0;
			int count = 273;
			float avg = 0;

			for (int m = -neighbor; m <= neighbor; m++) {
				for (int n = -neighbor; n <= neighbor; n++) {
					int row = i + m;
					int col = j + n;

					/* Check if neighborhood pixels are within image boundary */
					if (row > 0 && col > 0 && row < src.NR && col < src.NC) {
						sum += (src(row, col) * gaussian_filter[m + neighbor][n + neighbor]);

					}
				}
			}

			avg = float(sum / count);
			if(avg < threshold) {
				tgt(i, j) = 0;
			}
		}
	}

	sprintf(tgtName, "%s_smooth.pgm", srcName);
	tgt.save(tgtName);

	return 0;
}
