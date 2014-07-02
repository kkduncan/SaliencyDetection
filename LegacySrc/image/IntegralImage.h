/*
 * IntegralImage.h
 *
 *  Created on: Apr 9, 2011
 *      Author: kesterduncan
 *
 *  This class creates the Integral image of an image represented by
 *  an RImage.
 *
 *  Copyright (c) 2011 Kester Duncan
 */

#ifndef INTEGRALIMAGE_H_
#define INTEGRALIMAGE_H_

#include <stdlib.h>
#include <stdio.h>

#include "RImage.h"

class IntegralImage {

public:
	IntegralImage();
	virtual ~IntegralImage();
	void getIntegralImage(const RImage& orig, RImage &intImg);
};

#endif /* INTEGRALIMAGE_H_ */
