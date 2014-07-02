/*
 * GroupElement.h
 *
 *  Created on: Dec 17, 2009
 *      Author: kesterduncan
 */

#ifndef GROUPELEMENT_H_
#define GROUPELEMENT_H_

/**
 * Element encapsulates the low-level properties for an edge
 * pixel that we are interested in, particularly its location and
 * its gradient angle.
 */
typedef struct GroupElement {
	/*
	 * For Euclidean Distance
	 */
	int x;
	int y;
	/*
	 * Gradient Angle
	 */
	double angle;
	/*
	 * Gradient Magnitude
	 */
	double magnitude;

} GroupElement;


#endif /* GROUPELEMENT_H_ */
