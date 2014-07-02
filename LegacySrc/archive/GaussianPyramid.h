/*
 * K. Duncan
 * University of South Florida
 * Copyright (c) 2010
 *
 * -- GaussianPyramid.h --
 *
 * Header file for the Gaussian Pyramid technique outlined in the
 * work: "The Laplacian Pyramid as a Compact Image Code", Burt P.,
 * Adelson E., IEEE Transactions on Communications, April 1983
 *
 */

#ifndef GAUSSIANPYRAMID_H_
#define GAUSSIANPYRAMID_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image.h"

/*
 * This constant specifies the smallest image dimension that can be
 * used to create Gaussian Pyramids. This pyramid definition only
 * uses 5 levels.
 */
#define DIM_LIMIT 2

/*
 * This constant specifies the number of levels in this Gaussian pyramid
 */
#define NUM_LEVELS 5

class GaussianPyramid {

public:
	/*
	 * Keeps a record of the relevant dimensions of the images
	 * at each pyramid level in row, column format
	 */
	int dimensions[NUM_LEVELS][2];

	/*
	 * Stores the images at each pyramid level
	 */
	Image levels[NUM_LEVELS];

	void setBaseLevel(const Image &img);
	void reduceImageToLevel(const int srcLevel, const int tgtLevel);
	void expandImageToLevel(const int srcLevel, const int tgtLevel);
	void reduceAll();
	void expandAll();
	void processPyramid();


	GaussianPyramid();
	GaussianPyramid(Image src);
	virtual ~GaussianPyramid();
};

#endif /* GAUSSIANPYRAMID_H_ */
