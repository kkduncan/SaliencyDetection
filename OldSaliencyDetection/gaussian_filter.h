/*
 * Copyright (c) 2010
 * K. Duncan
 * University of South Florida
 *
 * This header file contains functions for performing Gaussian smoothing of an image
 *
 *	DISCLAIMER:
 *	==============================================================
 *  This code was not designed well. It is not modular and doesn't
 *  obey the Object Oriented Principles. Therefore it may contain
 *  many bugs. Use at your own risk!
 *  ==============================================================
 *
 */

#ifndef GAUSSIAN_FILTER_H_
#define GAUSSIAN_FILTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void gaussianSmooth(unsigned int *image, int rows, int cols, float sigma, int **smoothedim);
void makeGaussianKernel(float sigma, float **kernel, int *windowsize);


#endif /* GAUSSIAN_FILTER_H_ */
