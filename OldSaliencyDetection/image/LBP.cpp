/*
 * LBP.cpp
 *
 *  Created on: Jul 17, 2011
 *      Author: kesterduncan
 */

#include "LBP.h"

LBP::LBP() {
	/*
	 * Default values for LBP
	 */
	pixels = 8;
	radius = 1;

}


LBP::LBP(int p, int r) {

	if (p >= 4) {
		this->pixels = p;
	} else {
		this->pixels = 8;
	}

	if (r > 0) {
		this->radius = r;
	} else {
		this->radius = 1;
	}
}


LBP::~LBP() {
	pixels = 0;
	radius = 0;
}


void LBP::getLBP (const RImage& src, int* results) {
	if (results != NULL) {
		for (int i=0; i < src.numRows; i++) {
			for (int j=0; j < src.numCols; j++) {
				int LBP = 0;

				for (int p=0; p < pixels; p++) {
					int tgtRow = (int) (-radius * sin ((6.283185307 * p) / pixels) + i);
					int tgtCol = (int) (radius * cos ((6.283185307 * p) / pixels) + j);

					if (tgtRow >= 0 && tgtRow < src.numRows && tgtCol >= 0 && tgtCol < src.numCols) {
						if (tgtRow != i && tgtCol != j) {
							int diff;
							if (!src.isColorImage) {
								diff = src(tgtRow, tgtCol) - src(i, j);
							} else {
								diff = (src(tgtRow, tgtCol, RED) + src(tgtRow, tgtCol, GREEN) + src(tgtRow, tgtCol, BLUE) / 3) -
										(src(i, j, RED) + src(i, j, GREEN) + src(i, j, BLUE) / 3);
							}

							if (diff >= 0) {
								LBP += pow(2.0, p);
							}
						}
					}
				}

				int ij = (i * src.numCols) + j;
				results[ij] = LBP;
			}
		}
	} else {
		fprintf(stderr, "ERROR: Cannot perform LBP processing. NULL array provided.\n");
		exit(1);
	}
}


void LBP::getOriginalLBP (const RImage& src, int* results) {
	if (results != NULL) {
		int vals[][3] = {{1, 2, 4}, {8, 0, 16}, {32, 64, 128}};

		for (int i=0; i < src.numRows; i++) {
			for (int j=0; j < src.numCols; j++) {
				int LBP = 0;

				int neighbor = 1;

				for (int m = -neighbor; m <= neighbor; m++) {
					for (int n = -neighbor; n <= neighbor; n++) {
						int tgtRow = i + m;
						int tgtCol = j + n;

						if (tgtRow >= 0 && tgtRow < src.numRows && tgtCol >= 0 && tgtCol < src.numCols) {
							int diff;
							if (!src.isColorImage) {
								diff = src(tgtRow, tgtCol) - src(i, j);
							} else {
								diff = (src(tgtRow, tgtCol, RED) + src(tgtRow, tgtCol, GREEN) + src(tgtRow, tgtCol, BLUE) / 3) -
										(src(i, j, RED) + src(i, j, GREEN) + src(i, j, BLUE) / 3);
							}

							if (diff >= 0) {
								LBP += vals[m + 1][n + 1];
							}
						}

					}
				}

				int ij = (i * src.numCols) + j;
				results[ij] = LBP;
			}
		}
	} else {
		fprintf(stderr, "ERROR: Cannot perform LBP processing. NULL array provided.\n");
		exit(1);
	}
}
