/*
 * LBP.h
 *
 *  Created on: Jul 17, 2011
 *      Author: kesterduncan
 *
 *  Calculates the Local Binary Pattern of an image (RImage)
 */

#ifndef LBP_H_
#define LBP_H_

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "RImage.h"


class LBP {

	int pixels, radius;

public:
	LBP ();
	LBP (int p, int r);
	virtual ~LBP();
	void getLBP (const RImage& src, int* results);
	void getOriginalLBP (const RImage& src, int* results);
};

#endif /* LBP_H_ */
