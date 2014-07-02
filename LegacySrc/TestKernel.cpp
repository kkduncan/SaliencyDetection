/*
 * TestKernel.cpp
 *
 *  Created on: Mar 26, 2010
 *      Author: kesterduncan
 */

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include "Group.h"

using namespace std;

int main(int argc, char *argv[]) {

	Group g;
	g.init(1, 3, 3);

	srand(time(NULL));

	for (int i = 1; i <= 3; i++) {
		for (int j = 1; j <= 3; j++) {
			g.addElement(i, j, (M_PI / i), 1);
		}
	}

	g.buildHistogram();
	fprintf(stderr, "Real entropy: %2.5f\n", g.calculateRenyiEntropy());

	fprintf(stderr, "Estimation: %2.5f\n", g.estimateEntropy());

	return 0;
}
