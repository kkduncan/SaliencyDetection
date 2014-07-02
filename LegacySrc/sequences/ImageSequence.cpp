/*
 * ImageSequence.cpp
 *
 *  Created on: Jun 7, 2011
 *      Author: kesterduncan
 */

#include "ImageSequence.h"

ImageSequence::ImageSequence() {
	numRows = numCols = numFrames = 0;
	cacheSize = 0;
	frameCache = NULL;
	frameGradients = NULL;
	spatialMask = NULL;
	temporalMask = NULL;
}


ImageSequence::~ImageSequence() {
	flushMemory();
}


void ImageSequence::flushMemory() {
	for (int k = 0; k < this->numFrames; k++) {
		if (frameGradients[k] != NULL){
			free(frameGradients[k]);
		}
	}

	if (frameGradients != NULL) {
		free(frameGradients);
	}

	if (frameCache != NULL) {
		free(frameCache);
	}

	if (spatialMask != NULL) {
		free(spatialMask);
	}

	if (temporalMask != NULL) {
		free(temporalMask);
	}

	numRows = numCols = numFrames = 0;
}


void ImageSequence::initialize(double sXY, double sT, int startFrame, int endFrame,
		char sequenceDir[255], char fileStem[255]) {
	this->spatialSigma = sXY;
	this->temporalSigma = sT;
	this->startFrameID = startFrame;
	this->endFrameID = endFrame;
	this->numFrames = (endFrameID - startFrameID) + 1;

	strcpy(this->sequenceDir, sequenceDir);
	strcpy(this->sequenceFileStem, fileStem);

	computeGaussianMask(sXY, &spatialMask);
	computeGaussianMask(sT, &temporalMask);
	this->spatialMaskWidth = (2 * (3 * (int) sXY)) + 1;
	this->temporalMaskWidth = (2 * (3 * (int) sT)) + 1;
	cacheSize = temporalMaskWidth + 3;
	window = temporalMaskWidth / 2;

	if ((frameCache = (RImage*) calloc(cacheSize, sizeof(RImage))) == NULL) {
		fprintf(stderr, "ERROR: Cannot allocate space for frame cache.\n");
		exit(1);
	}

	/* Retrieve Frame Dimensions */
	formFrameName(startFrameID);
	frameCache[0].read(sequenceFrameName);
	this->numRows = frameCache[0].numRows;
	this->numCols = frameCache[0].numCols;

	/* Set up frame gradients */
	if (numFrames > 0) {
		frameGradients = (ThreeDGradient**) calloc(numFrames, sizeof(ThreeDGradient*));
		if (frameGradients == NULL) {
			fprintf(stderr, "ERROR: Cannot allocate space for frame gradients array.\n");
			exit(1);
		}

		for (int i = 0; i < numFrames; i++) {
			frameGradients[i] = (ThreeDGradient*) calloc(this->numRows * this->numCols, sizeof(ThreeDGradient));

			if (frameGradients[i] == NULL) {
				printf("ERROR: Cannot allocate memory for the specified frame gradient array.\n");
				exit(1);
			}
		}
	} else {
		printf("ERROR: The number of frames must exceed 3.\n");
		exit(1);
	}
}


void ImageSequence::detectEdges() {
	for (int t = startFrameID, frameNum = 0; t <= endFrameID; t++, frameNum++) {
		loadBlock(t);
		smoothXYT(t);
		computeGradient(t);
		getGradients(frameNum);
	}
}


void ImageSequence::loadBlock(int frameNum) {
	int t, k;

	timeIndex = frameNum;
	if (frameNum == startFrameID) {
		/* This goes for the frame cache size, in this case 16 */
		for (t = frameNum - window - 2, k = 0; t <= frameNum + window + 1; t++, k++) {
			if (t <= startFrameID) {
				formFrameName(startFrameID);

			} else if (t >= endFrameID) {
				formFrameName(endFrameID);

			} else {
				formFrameName(t);
			}

			frameCache[k].read(sequenceFrameName);
		}

		this->numRows = frameCache[0].numRows;
		this->numCols = frameCache[0].numCols;

	} else {
		for (k = 0; k < cacheSize - 1; k++) {
			frameCache[k] = frameCache[k + 1];
		}

		t = frameNum + window + 1;
		if (t >= endFrameID) {
			formFrameName(endFrameID);

		} else {
			formFrameName(t);

		}

		frameCache[k].read(sequenceFrameName);
	}

}


void ImageSequence::smoothXY(RImage &img) {
	RImage output;
	int offset, i, j, k;

	output.initialize(1, max(img.numCols, img.numRows));
	offset = spatialMaskWidth / 2;

	/* First Smooth along the Rows */
	for (i = 0; i < this->numRows; i++) {
		for (j = 0; j < this->numCols; j++) {
			output(0, j) = 0;
			for (k = 0; k < spatialMaskWidth; k++) {
				output(0, j) += spatialMask[k] * img(i, max(0, min(img.numCols - 1, j - offset + k)));
			}
		}
		for (j = 0; j < img.numCols; j++)
			img(i, j) = output(0, j);
	}

	/* Smooth along the Columns */
	for (j = 0; j < this->numCols; j++) {
		for (i = 0; i < this->numRows; i++) {
			output(0, i) = 0;

			for (k = 0; k < spatialMaskWidth; k++) {
				output(0, i) += spatialMask[k] * img(max(0, min(img.numRows - 1, i - offset + k)), j);
			}
		}

		for (i = 0; i < this->numRows; i++) {
			img(i, j) = output(0, i);
		}
	}

	output.flushMemory();
}


void ImageSequence::smoothXYT(int frameNum) {
	int k;
	if (timeIndex != frameNum) {
		fprintf(stderr, "Wrong time frame: Expected %d, Got Request for %d\n",
				timeIndex, frameNum);
		exit (1);
	}

	if (frameNum == startFrameID) {
		/* For the first frame, we have to smooth all the frames in x, y, and t directions*/

		for (k = 0; k < cacheSize; k++) {
			smoothXY(frameCache[k]);
		}

		for (k = 0; k < 4; k++) {
			smoothed[k].initialize(frameCache[k].numRows, frameCache[k].numCols);

			for (int i = 0; i < frameCache[k].numRows; i++) {
				for (int j = 0; j < frameCache[k].numCols; j++) {
					smoothed[k](i, j) = 0;
					for (int t = 0; t < temporalMaskWidth; t++) {
						smoothed[k](i, j) += frameCache[k + t](i, j) * temporalMask[t];
					}
				}
			}
		}
	} else {
		/* After the first frame, processing is incremental*/
		/* Last Frame */
		smoothXY(frameCache[cacheSize - 1]);
		for (k = 0; k < 3; k++)
			smoothed[k] = smoothed[k + 1];

		for (int i = 0; i < frameCache[0].numRows; i++) {
			for (int j = 0; j < frameCache[0].numCols; j++) {
				smoothed[k](i, j) = 0;
				for (int t = 0; t < temporalMaskWidth; t++) {
					smoothed[k](i, j) += frameCache[k + t](i, j) * temporalMask[t];
				}
			}
		}
	}
}


void ImageSequence::computeGradient(int frameNum) {
	if (timeIndex != frameNum) {
		fprintf(stderr, "Wrong time frame: Expected %d, Got Request for %d\n",
				timeIndex, frameNum);
		exit (1);
	}
	for (int t = 0; t < 3; t++) {
		if (frameNum == startFrameID) {
			this->delX[t].initialize(smoothed[0].numRows, smoothed[0].numCols);
			this->delY[t].initialize(smoothed[0].numRows, smoothed[0].numCols);
			this->delT[t].initialize(smoothed[0].numRows, smoothed[0].numCols);
		}

		for (int i = 1; i < smoothed[0].numRows - 1; i++) {
			for (int j = 1; j < smoothed[0].numCols - 1; j++) {

				/* Simply calculates the intensity changes */
				this->delT[t](i, j) = smoothed[t + 1](i, j) - smoothed[t](i, j);
				this->delY[t](i, j) = smoothed[t + 1](i, j) - smoothed[t + 1](i, j - 1);
				this->delX[t](i, j) = smoothed[t + 1](i, j) - smoothed[t + 1](i - 1, j);
			}
		}
	}
}


void ImageSequence::getGradients(int frameNum) {
	for (int i = 0; i < this->numRows; i++) {
		for (int j = 0; j < this->numCols; j++) {
			int ij = (i * this->numCols) + j;
			frameGradients[frameNum][ij].delX = this->delX[1](i, j);
			frameGradients[frameNum][ij].delY = this->delY[1](i, j);
			frameGradients[frameNum][ij].delT = this->delT[1](i, j);
		}
	}
}


void ImageSequence::formFrameName(int t) {
	/* Forms the name of the file using the FileStem and the Frame time id */
	sprintf(sequenceFrameName, "%s/%s%04d.pgm", this->sequenceDir, this->sequenceFileStem, t);

}


void ImageSequence::computeGaussianMask(double sigma, double** mask) {
	double sum = 0.0;
	int i, W;

	W = (int) (3 * sigma);
	if (W > 500) {
		fprintf(stderr, "ERROR: The Gaussian operator mask is too large (%s, %d)\n", __FILE__, __LINE__);
		exit(1);
	}

	*mask = (double*) calloc(2 * W + 1, sizeof(double));
	if (*mask == NULL) {
		fprintf(stderr, "ERROR: Cannot allocate space for Gaussian mask.\n");
		exit(1);
	}

	//mask.initialize(1, 2 * W + 1);
	for (i = W; i >= -W; i--) {
		(*mask)[i + W] = exp(-i * i / (2 * sigma * sigma));
		sum += (*mask)[i + W];
	}

	for (i = W; i >= -W; i--) {
		(*mask)[i + W] /= sum;
	}
}
