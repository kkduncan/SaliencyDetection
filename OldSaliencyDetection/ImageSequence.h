/*
 * ImageSequence.h
 *
 *  Created on: Jun 7, 2011
 *      Author: kesterduncan
 *
 *	DISCLAIMER:
 *	==============================================================
 *  This code was not designed well. It is not modular and doesn't
 *  obey the Object Oriented Principles. Therefore it may contain
 *  many bugs. Use at your own risk!
 *  ==============================================================
 */

#ifndef IMAGESEQUENCE_H_
#define IMAGESEQUENCE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../image/RImage.h"

typedef struct gradient {
	double delX;
	double delY;
	double delT;

} ThreeDGradient;

class ImageSequence {

private:
	int max(int i, int j) {
		if (i > j)
			return (i);
		else
			return (j);
	}

	int min(int i, int j) {
		if (i <= j)
			return (i);
		else
			return (j);
	}

public:
	int numRows;
	int numCols;
	int numFrames;
	int startFrameID;
	int endFrameID;
	int timeIndex;
	int cacheSize;
	int window;
	int spatialMaskWidth;
	int temporalMaskWidth;
	double spatialSigma;
	double temporalSigma;
	char sequenceFrameName[255];
	char sequenceDir[255];
	char sequenceFileStem[255];

	double* spatialMask;
	double* temporalMask;

	RImage smoothed[4];
	RImage delX[3];
	RImage delY[3];
	RImage delT[3];
	RImage *frameCache;
	ThreeDGradient **frameGradients;

	void initialize(double sXY, double sT, int startFrame, int endFrame, char sequenceDir[255], char fileStem[255]);
	void detectEdges();
	void flushMemory();
	void formFrameName(int frameNum);
	void computeGaussianMask(double sigma, double **mask);
	void smoothXY(RImage &img);
	void smoothXYT(int frameNum);
	void loadBlock(int frameNum);
	void computeGradient(int frameNum);
	void getGradients(int frameNum);

	ImageSequence();
	virtual ~ImageSequence();
};

#endif /* IMAGESEQUENCE_H_ */
