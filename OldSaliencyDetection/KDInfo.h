/*
 * KDInfo.h
 *
 * Kernel Density Information
 *
 *  Created on: Nov 14, 2010
 *      Author: kesterduncan
 *
 *  Copyright (c) 2010
 *
 *	DISCLAIMER:
 *	==============================================================
 *  This code was not designed well. It is not modular and doesn't
 *  obey the Object Oriented Principles. Therefore it may contain
 *  many bugs. Use at your own risk!
 *  ==============================================================
 */

#ifndef KDINFO_H_
#define KDINFO_H_

#include <cstdlib>

/*
 * Structure to store per pixel information required
 * for performing iterative kernel density estimation
 * for local neighborhoods
 */
typedef struct KDInfo {
	float kernelSum;
	float entropy;
	float magSum1;
	float magSum2;
	int sampleCount;

} KDInfo;


/*
 * Structure for storing the information that is returned
 * after calculating the kernel sum
 */
typedef struct KSumInfo {
	float kernelSum;
	float magSum1;
	float magSum2;

} KSumInfo;


/*
 * Structure for storing kernel descriptor information
 */
typedef struct KernelDescriptors {
	double gradient;
	double lbp;
	double color;

} KernelDescriptors;


/*
 * 2D Pixel Location
 */
typedef struct Location {
	int col;
	int row;

} Location;


/*
 * 3D Pixel Location
 */
typedef struct Location3D {
	int x;
	int y;
	int t;
} Location3D;


/*
 * Stores the applicable bounds region of the 4 pixels forming
 * the two samples
 */
typedef struct Bounds {
	Location topLeft;
	Location topRight;
	Location botLeft;
	Location botRight;

} Bounds;



/*
 * Initialize an array of KDInfos
 */
void initializeKDInfo (KDInfo *pi, const int size) {
	if (pi != NULL) {
		for (int i = 0; i < size; i++) {
			pi[i].entropy = 0;
			pi[i].kernelSum = 0;
			pi[i].sampleCount = 0;
			pi[i].magSum1 = 0;
			pi[i].magSum2 = 0;

		}
	} else {
		fprintf(stderr, "Error: KDInfo structure cannot be init. It is NULL.\n");
		exit(1);
	}
}





#endif /* KDINFO_H */
