/*
 * Copyright (c) 2010 Kester K. Duncan.
 * All rights reserved.
 *
 *  Created on: Nov 27, 2009
 *      Author: kesterduncan
 *
 *  REM class header file:
 *  	This class forms a relational distribution based on the
 *  	pairwise comparisons of its members, then calculates the
 *  	entropy of that distribution.
 *
 *  $Id$
 */

#ifndef REM_H_
#define REM_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <list>
#include "histogram/histogram.h"
#include "image/RImage.h"

using namespace std;

/**
 * Element encapsulates the low-level properties for an edge
 * pixel that we are interested in, particularly its location and
 * its gradient angle.
 */
typedef struct GroupElement {
	/*
	 * For Euclidean Distance
	 */
	int x;
	int y;
	/*
	 * Intensity Gradient Angle
	 */
	double angle;
	/*
	 * Intensity Gradient Magnitude
	 */
	double magnitude;

} GroupElement;


class REM {

private:
	double minDist;
	double maxDist;
	double minAngle;
	double maxAngle;
	double distanceThreshold;
	double distanceBinWidth;
	double angleBinWidth;
	double samplePercentage;
	double log2(double value);

public:
	int binDimension;
	double entropy;
	histogram *hist;
	list<GroupElement> members;
	GroupElement *membersArray;

	void init(int width, int height);
	void init(int imageWidth, int imageHeight, int dimension);
	void init(int width, int height, int dimension, int percentage);
	void addElement(int x, int y, double angle);
	void addElement(int x, int y, double angle, double magnitude);
	void initializeHistogram();
	void buildHistogram();
	void buildSampledHistogram();
	void clearHistogram();
	void processHistogram();
	void updateHistogram(double distance, double angle);
	void updateHistogram(double distance, double angle, double magnitude);
	void printHistogramMatrix(char * filename);
	void printHistogramFile(char * filename);
	void readHistogramFile(char * filename);
	void drawHistogram(char * filename);

	/*
	 * Kernel Density Estimation
	 */
	double estimateEntropy();
	double estimateEntropyFromSample();

	/*
	 * Entropy
	 */
	double calculateRenyiEntropy();
	double calculateShannonEntropy();
	void copyElementsToArray();

	REM();
	REM(int width, int height, int dimension);
	REM(int width, int height, int dimension, int percentage);
	virtual ~REM();

};

#endif /* REM_H_ */
