/*
 * Copyright (c) 2010
 * University of South Florida
 *
 * CannyEdge.h
 *
 * Header file for Mike Heath's Canny edge detector implementation
 *
 */

#ifndef CANNYEDGE_H_
#define CANNYEDGE_H_

#define VERBOSE 0
#define BOOSTBLURFACTOR 1.0 /* was 90.0 */

/*
 * This assumes that the grad_dirs, mags arrays are already allocated.
 */
void
canny(unsigned char *image, int rows, int cols, float sigma, float *grad_dirs, float *mags);

void
gaussian_smooth(unsigned char *image, int rows, int cols, float sigma, short int **smoothedim);

void
make_gaussian_kernel(float sigma, float **kernel, int *windowsize);

void
derrivative_x_y(short int *smoothedim, int rows, int cols, short int **delta_x, short int **delta_y);

void
radian_direction(short int *delta_x, short int *delta_y, int rows, int cols, float **dir_radians, int xdirtag, int ydirtag);

double
angle_radians(double x, double y);


#endif /* CANNYEDGE_H_ */
