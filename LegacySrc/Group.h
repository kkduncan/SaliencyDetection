/*
 * Copyright (c) 2010 Kester K. Duncan.
 * All rights reserved.
 *
 *  Created on: Nov 27, 2009
 *      Author: kesterduncan
 *
 *  Group class header file:
 *  	This class forms a relational distribution based on the
 *  	pairwise comparisons of its members, then calculates the
 *  	entropy of that distribution.
 *
 *  $Id$
 */

#ifndef GROUP_H_
#define GROUP_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <time.h>
#include <list>
#include <gsl/gsl_histogram2d.h>
#include "GroupElement.h"
#include "image.h"

using namespace std;

class Group {

private:
	double minDist;
	double maxDist;
	double minAngle;
	double maxAngle;
	double distanceThreshold;
	double distanceBinWidth;
	double angleBinWidth;
	double samplePercentage;
	double groupLog2(double value);

public:
	int id;
	int rank;
	int binDimension;
	int elemCount;
	double entropy;
	gsl_histogram2d *histogram;
	list<GroupElement> members;
	GroupElement *membersArray;

	void init(int index, int width, int height);
	void init(int index, int width, int height, int dimension, int percentage);
	void init(int index, int imageWidth, int imageHeight, int dimension);
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
	double estimateEntropy();
	double calculateRenyiEntropy();
	double calculateShannonEntropy();
	void copyElementsToArray();

	Group();
	virtual ~Group();

};

#endif /* GROUP_H_ */
