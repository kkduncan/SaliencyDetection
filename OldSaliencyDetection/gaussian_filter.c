/*
 * Copyright (c) 2010
 * K. Duncan
 * University of South Florida
 *
 * This implementation file contains functions for performing Gaussian smoothing of an image
 *
 */

#include "gaussian_filter.h"

/**
 * Blur an image with a Gaussian filter.
 * @param image
 * @param rows
 * @param cols
 * @param sigma
 * @param smoothedim
 */
void gaussianSmooth(unsigned int *image, int rows, int cols, float sigma, int **smoothedim) {
	int r, c, rr, cc, 	/* Counter variables. */
	windowsize, 		/* Dimension of the gaussian kernel. */
	center; 			/* Half of the windowsize. */

	float *tempim, 		/* Buffer for separable filter gaussian smoothing. */
	*kernel, 			/* A one dimensional gaussian kernel. */
	dot, 				/* Dot product summing variable. */
	sum; 				/* Sum of the kernel weights variable. */

	/****************************************************************************
	 * Create a 1-dimensional gaussian smoothing kernel.
	 ****************************************************************************/
	makeGaussianKernel(sigma, &kernel, &windowsize);
	center = windowsize / 2;

	/****************************************************************************
	 * Allocate a temporary buffer image and the smoothed image.
	 ****************************************************************************/
	if ((tempim = (float *) malloc(rows * cols * sizeof(float))) == NULL) {
		fprintf(stderr, "Error allocating the buffer image.\n");
		exit(1);
	}
	if (((*smoothedim) = (int *) malloc(rows * cols * sizeof(int))) == NULL) {
		fprintf(stderr, "Error allocating the smoothed image.\n");
		exit(1);
	}

	/****************************************************************************
	 * Blur in the x - direction.
	 ****************************************************************************/
	for (r = 0; r < rows; r++) {
		for (c = 0; c < cols; c++) {
			dot = 0.0;
			sum = 0.0;
			for (cc = (-center); cc <= center; cc++) {
				if (((c + cc) >= 0) && ((c + cc) < cols)) {
					dot += (float) image[r * cols + (c + cc)] * kernel[center + cc];
					sum += kernel[center + cc];
				}
			}
			tempim[r * cols + c] = dot / sum;
		}
	}

	/****************************************************************************
	 * Blur in the y - direction.
	 ****************************************************************************/
	/*
	 * Get the maximum value for normalization
	 */
	int maximum = -999;

	for (c = 0; c < cols; c++) {
		for (r = 0; r < rows; r++) {
			sum = 0.0;
			dot = 0.0;
			for (rr = (-center); rr <= center; rr++) {
				if (((r + rr) >= 0) && ((r + rr) < rows)) {
					dot += tempim[(r + rr) * cols + c] * kernel[center + rr];
					sum += kernel[center + rr];
				}
			}
			int value = (int) (dot * 1.0 / sum + 0.5);
			(*smoothedim)[r * cols + c] = value;

			if (value > maximum) {
				maximum = value;
			}
		}
	}

	// Normalization of the smoothed image
	for (c = 0; c < cols; c++) {
		for (r = 0; r < rows; r++) {
			float value = (float) (*smoothedim)[r * cols + c];
			(*smoothedim)[r * cols + c] = (int) ((value / maximum) * 255);
		}
	}

	free(tempim);
	free(kernel);
}


/**
 * Create a one-dimensional Gaussian kernel
 * @param sigma
 * @param kernel
 * @param windowsize
 */
void makeGaussianKernel(float sigma, float **kernel, int *windowsize) {
	int i, center;
	float x, fx, sum = 0.0;

	/* changed from 2.5 to 4.0, (1 + 2 * ) removed */
	*windowsize = (int) (1 + ceil(4.0 * sigma));
	center = (*windowsize) / 2;

	if ((*kernel = (float *) malloc((*windowsize) * sizeof(float))) == NULL) {
		fprintf(stderr, "Error callocing the gaussian kernel array.\n");
		exit(1);
	}

	for (i = 0; i < (*windowsize); i++) {
		x = (float) (i - center);
		fx = pow(2.71828, -0.5 * x * x / (sigma * sigma)) / (sigma * sigma * sqrt(6.2831853));
		(*kernel)[i] = fx;
		sum += fx;
	}

	for (i = 0; i < (*windowsize); i++) {
		(*kernel)[i] /= sum;
	}

}

