/*
 * ScaleSpace.h
 *
 *  Created on: Mar 8, 2010
 *      Author: kesterduncan
 */

#ifndef SCALESPACE_H_
#define SCALESPACE_H_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "image.h"
#include "Group.h"

class ScaleSpace {

private:
	void getFirstDerivative();
	void getSecondDerivative();

public:
	static const int MIN_SCALE = 3;
	static const int MAX_SCALE = 129;
	static const int NUM_SCALES = 7;
	int maxScale;
	float scaleEntropies[NUM_SCALES];
	Group scaleGroups[NUM_SCALES];
	float firstDerivative[NUM_SCALES];
	float secondDerivative[NUM_SCALES];

	ScaleSpace();
	virtual ~ScaleSpace();
	int getMaxScale();

};

#endif /* SCALESPACE_H_ */
