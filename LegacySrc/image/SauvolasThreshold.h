/*
 * SauvolasThreshold.h
 *
 *  Created on: May 16, 2011
 *      Author: kesterduncan
 */

#ifndef SAUVOLASTHRESHOLD_H_
#define SAUVOLASTHRESHOLD_H_

#include <stdlib.h>
#include <stdio.h>
#include "RImage.h"
#include "IntegralImage.h"

class SauvolasThreshold {
public:
	SauvolasThreshold();
	virtual ~SauvolasThreshold();
	void threshold(RImage& src);

};

#endif /* SAUVOLASTHRESHOLD_H_ */
