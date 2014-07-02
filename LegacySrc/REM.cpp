/*
 * Copyright (c) 2010 Kester K. Duncan.
 * All rights reserved.
 *
 *  Created on: Nov 27, 2009
 *
 *  REM class implementation:
 *  	This class forms a relational distribution based on the
 *  	pairwise comparisons of its members, then calculates the
 *  	entropy of that distribution.
 */

#include "REM.h"

namespace {
	const char rcs_id[] = "$Id$";
}


REM::REM() {
	entropy = 0;
	membersArray = NULL;
	srand(time(NULL));
}


REM::REM(int columnWidth, int rowHeight, int dimension) {
	REM(columnWidth, rowHeight, dimension, 100);
}


REM::REM(int width, int height, int dimension, int percentage) {
	entropy = 0;
	membersArray = NULL;

	minDist = 0;
	maxDist = sqrt(pow(width, 2) + pow(height, 2));
	minAngle = 0;
	maxAngle = M_PI;
	binDimension = dimension;
	distanceBinWidth = (double)(maxDist / dimension);
	angleBinWidth = (double)(maxAngle / dimension);
	samplePercentage = (double) percentage / 100;

	if (binDimension > 1) {
		initializeHistogram();
	} else {
		fprintf(stderr, "ERROR: Invalid histogram dimension provided\n");
		exit(-1);
	}

	srand(time(NULL));

}


REM::~REM() {
	histogram_free(hist);
	if (membersArray != NULL) {
		free(membersArray);
	}
	membersArray = NULL;
}


/**
 * Initialize a group that corresponds to the whole image
 */
void REM::init(int columnWidth, int rowHeight, int dimension) {
	init(columnWidth, rowHeight, dimension, 100);
}



/**
 * Initialize a group that corresponds to an image
 */
void REM::init(int width, int height, int dimension, int percentage) {
	entropy = 0;
	minDist = 0;
	maxDist = sqrt(pow(width, 2) + pow(height, 2));
	minAngle = 0;
	maxAngle = M_PI;
	binDimension = dimension;
	distanceBinWidth = (double)(maxDist / dimension);
	angleBinWidth = (double)(maxAngle / dimension);
	samplePercentage = (double) percentage / 100;

	if (binDimension > 1) {
		initializeHistogram();
	} else {
		fprintf(stderr, "ERROR: Invalid histogram dimension provided\n");
		exit(-1);
	}

}


/**
 * Adds an element to the members list. These elements have equal weight
 * in building the histogram.
 */
void REM::addElement(int x, int y, double angle) {
	GroupElement newMember;
	newMember.x = x;
	newMember.y = y;
	newMember.angle = angle;
	newMember.magnitude = 1;
	members.push_back(newMember);
}


/**
 * Adds an element to the members list. The magnitude can vary.
 */
void REM::addElement(int x, int y, double angle, double magnitude) {
	GroupElement newMember;
	newMember.x = x;
	newMember.y = y;
	newMember.angle = angle;
	newMember.magnitude = magnitude;
	members.push_back(newMember);
}


/**
 * Initialize the relational histogram for the distance and angle parameters
 */
void REM::initializeHistogram() {
	hist = histogram_alloc(binDimension, binDimension);
	histogram_set_ranges_uniform(hist, minDist, maxDist, minAngle, maxAngle);
}


void REM::buildHistogram() {
	list<GroupElement>::iterator xElem, yElem;
	double distance, angle, magnitude;

	if (members.size() > 0) {
		for (xElem = members.begin(); xElem != members.end(); xElem++) {
			for (yElem = members.begin(); yElem != members.end(); yElem++) {

				distance = sqrt(pow((xElem->x - yElem->x), 2) + pow((xElem->y - yElem->y), 2));
				angle = fmod(((2 * M_PI) + (xElem->angle - yElem->angle)), (2 * M_PI));

				if (angle > M_PI) {
					angle = (2 * M_PI) - angle;
				}

				magnitude = fabs(xElem->magnitude - yElem->magnitude); //changed from fabs
				updateHistogram(distance, angle, magnitude);
			}
		}
	} else {
		fprintf(stderr, "REM_ERROR: No members found to build a distribution.\n");
		exit(1);
	}
}


void REM::buildSampledHistogram() {
	double distance, angle, magnitude;
	int randElem1, randElem2;
	int count, samples;

	samples = (int)((samplePercentage * members.size()) * 2);
	count = 0;

	while (count < samples) {
		randElem1 = rand() % members.size();
		randElem2 = rand() % members.size();

		distance = sqrt(pow(membersArray[randElem1].x - membersArray[randElem2].x, 2))
					+ sqrt(pow(membersArray[randElem1].y - membersArray[randElem2].y, 2));
		angle = fmod(((2 * M_PI) + (membersArray[randElem1].angle - membersArray[randElem2].angle)), (2 * M_PI));

		if (angle > M_PI) {
			angle = (2 * M_PI) - angle;
		}

		magnitude = fabs(membersArray[randElem1].magnitude - membersArray[randElem2].magnitude);
		updateHistogram(distance, angle, magnitude);
		count++;
	}
}


void REM::clearHistogram() {
	histogram_reset(hist);
}


/**
 * Update a histogram bin based on distance and angle with equal weights
 * for each group element
 */
void REM::updateHistogram(double distance, double angle) {
	histogram_increment(hist, distance, angle);
}

/**
 * Update a histogram bin based on distance and angle with the weights
 * given by the root of the product of two element magnitudes
 */
void REM::updateHistogram(double distance, double angle, double magnitude) {
	/*
	 * Since the magnitude differences are used to weight the histogram bins,
	 * we cannot allow differences of 0, because the histograms would
	 * have no values in the bins. Therefore, when there is a zero magnitude
	 * difference, we use the value of 1 instead in order to allow
	 * the distribution to be constructed properly.
	 */
	//(magnitude == 0) ? magnitude = 0.1 : NULL;
	histogram_accumulate(hist, distance, angle, magnitude);
}


/**
 * Calculate the Renyi entropy of this group's relational histogram
 */
double REM::calculateRenyiEntropy() {
	int i, j;
	double totalElems = 0;
	double totalBins = binDimension * binDimension;

	totalElems = histogram_sum(hist);

	if (entropy != 0) entropy = 0;

	if (totalElems != 0 && histogram_max_val(hist) != histogram_min_val(hist)) {
		for(i = 0; i < binDimension; i++) {
			for(j = 0; j < binDimension; j++) {
				double e = 0;
				double p = (double) (histogram_get(hist, i, j) / totalElems);

				if(p > 0) {
					e = p * p;
				} else {
					e = 0;
				}
				entropy += e;
			}
		}
	} else {
		entropy = -1;
		return entropy;
	}

	/* Extra term to penalize small groups - Abe's term */
	entropy = -1 * (log2(entropy));
	entropy += (( (totalBins) - 1) * ((log2(exp(1))) / (2 * totalElems))) / totalBins;

	if(entropy >= 0) {
		entropy /= (double) (log2(totalBins));
	}

	return entropy;
}


/**
 * Calculate the Shannon entropy of this group's relational histogram
 */
double REM::calculateShannonEntropy() {
	int i, j;
	double totalElems = 0;
	double totalBins = binDimension * binDimension;

	totalElems = histogram_sum(hist);

	if (totalElems != 0) {
		for(i = 0; i < binDimension; i++) {
			for(j = 0; j < binDimension; j++) {
				double e = 0;
				double p = (double) (histogram_get(hist, i, j) / totalElems);

				if(p > 0) {
					e = -p * log2(p);
				} else {
					e = 0;
				}
				entropy += e;
			}
		}
	}
	else {
		/*
		 * This is a hack for constant intensity regions
		 *
		 * Using this, the entropy value computed after normalization
		 * below is high (close to 1)
		 */
		entropy = 6.5;
	}

	/* Abe's term to penalize small groups */
//	entropy += (((totalBins) - 1) * ((log2(exp(1))) / (2 * totalElems))) / totalBins;

	if(entropy >= 0) {
		entropy /= (double) (log2(totalBins));
	}

	return entropy;
}

/**
 * This function works by estimating the probability density function
 * for the local pixel neighborhood of a central pixel. The neighborhood
 * is represented by the list "members" which is part of the REM class.
 *
 */
double REM::estimateEntropy() {
	list<GroupElement>::iterator xElem1, yElem1, xElem2, yElem2;
	double sampleDistance1, sampleAngle1;
	double sum = 0;
	double distanceKernel = 0, angleKernel = 0;
	double estimation = 0, entropyValue = 0;
	double numSamples = pow(members.size(), 4); // Number of sample pairs
	double twoPI = 2 * M_PI;
	double dNorm = (sqrt(twoPI) * distanceBinWidth);
	double aNorm = (sqrt(twoPI) * angleBinWidth);

	for (xElem1 = members.begin(); xElem1 != members.end(); xElem1++) {
		for (yElem1 = members.begin(); yElem1 != members.end(); yElem1++) {

			sampleDistance1 = sqrt(pow((xElem1->x - yElem1->x), 2) + pow((xElem1->y - yElem1->y), 2));
			sampleAngle1 = fmod((twoPI + (xElem1->angle - yElem1->angle)), twoPI);

			if (sampleAngle1 > M_PI) {
				sampleAngle1 = twoPI - sampleAngle1;
			}

			double sampleDistance2, sampleAngle2;

			for (xElem2 = members.begin(); xElem2 != members.end(); xElem2++) {
				for (yElem2 = members.begin(); yElem2 != members.end(); yElem2++) {

					sampleDistance2 = sqrt(pow((xElem2->x - yElem2->x), 2) + pow((xElem2->y - yElem2->y), 2));
					sampleAngle2 = fmod((twoPI + (xElem2->angle - yElem2->angle)), twoPI);

					if (sampleAngle2 > M_PI) {
						sampleAngle2 = twoPI - sampleAngle2;
					}

					distanceKernel = (1/dNorm) * exp((pow(sampleDistance1 - sampleDistance2, 2) / (-2 * pow(distanceBinWidth, 2))));
					angleKernel = (1/aNorm) * exp((pow(sampleAngle1 - sampleAngle2, 2) / (-2 * pow(angleBinWidth, 2))));

					sum = sum + distanceKernel * angleKernel;
				}
			}
		}
	}

	estimation = sum / numSamples;

	/* Sanity check	 */
	if (estimation > 1 || estimation < 0) {
		printf("Estimation greater than 1 or less than 0 at pixel: %2.4f\n", estimation);
	}

	entropyValue = -1 * log2(estimation * estimation);

	/* Sanity check	 */
	if (entropyValue < 0) {
		entropyValue = -1;
	}

	return entropyValue;
}


/**
 * Perform Kernel Density Estimation on samples from a local
 * pixel neighborhood
 */
double REM::estimateEntropyFromSample() {
	double sum = 0, magSum1 = 0, magSum2 = 0;
	double distanceKernel = 0, angleKernel = 0;
	double estimation = 0, entropyValue = 0;
	const double twoPI = 6.283;
	double dNorm = (sqrt(twoPI) * distanceBinWidth);
	double aNorm = (sqrt(twoPI) * angleBinWidth);
	int numSamples = (int)((samplePercentage * members.size()));
	int i, j, count = 0;

	i = 0;
	while (i < numSamples) {
		int randElem1 = rand() % members.size();
		int randElem2 = rand() % members.size();
		double sampleDistance1, sampleAngle1, sampleMag1;

		sampleDistance1 = sqrt(pow((membersArray[randElem1].x - membersArray[randElem2].x), 2) +
							   pow((membersArray[randElem1].y - membersArray[randElem2].y), 2));
		sampleAngle1 = fmod((fabs(membersArray[randElem1].angle - membersArray[randElem2].angle)), 6.283);

		//if (sampleAngle1 > M_PI) {
		//	sampleAngle1 = twoPI - sampleAngle1;
		//}

		sampleMag1 = fabs(membersArray[randElem1].magnitude - membersArray[randElem2].magnitude);
		magSum1 += sampleMag1;

		j = 0;
		while (j < numSamples) {
			int randElem3 = rand() % members.size();
			int randElem4 = rand() % members.size();
			double sampleDistance2, sampleAngle2, sampleMag2;

			sampleDistance2 = sqrt(pow((membersArray[randElem3].x - membersArray[randElem4].x), 2) +
							       pow((membersArray[randElem3].y - membersArray[randElem4].y), 2));
			sampleAngle2 = fmod((fabs(membersArray[randElem3].angle - membersArray[randElem4].angle)), 6.283);

			//if (sampleAngle2 > M_PI) {
			//	sampleAngle2 = twoPI - sampleAngle2;
			//}

			sampleMag2 = fabs(membersArray[randElem3].magnitude - membersArray[randElem4].magnitude);
			magSum2 += sampleMag2;

			/*
			 * Kernels
			 */
			distanceKernel = (1/dNorm) * exp((pow(sampleDistance1 - sampleDistance2, 2) / (-2 * pow(distanceBinWidth, 2))));
			angleKernel = (1/aNorm) * exp((pow(sampleAngle1 - sampleAngle2, 2) / (-2 * pow(angleBinWidth, 2))));

			if (sampleMag1 > 0 && sampleMag2 > 0 && distanceKernel > 0 && angleKernel > 0) {
				sum += (sampleMag1 * sampleMag2 * distanceKernel * angleKernel);
			}

			count++;
			j++;
		}
		i++;
	}

	double totalMagSum = magSum1 * magSum2;

	/* Special Case Handling */
	if (totalMagSum <= 0) {
		totalMagSum = count;
	}
	if (sum < 0 || isnan(sum)/* || !isfinite(sum) */) {
		printf("Sum: %2.5f\n", sum);
		sum = 0;

	}

	estimation = (double) sum / totalMagSum;

	if (estimation <= 0 || isnan(estimation)/* || !isfinite(estimation) */) {
		/*
		 * This is a hack in case the estimation value is invalid.
		 * In ComputeSaliency.cpp we give pixels with an entropy value
		 * of -1 the maximum entropy found in the image.
		 */
		entropyValue = -1;
	} else {
		entropyValue = -1 * this->log2(estimation * estimation);
	}

	return entropyValue;
}


/**
 * Draw a PGM image of the histogram
 */
void REM::drawHistogram(char * filename) {
	RImage img;
	int i, j;
	double maxValue = histogram_max_val(hist);

	img.initialize(binDimension, binDimension);

	for(i = 0; i < binDimension; i++) {
		for (j = 0; j < binDimension; j++) {
			img(i, j, RED) = (int) ceil(((double) histogram_get(hist, i, j) / maxValue) * 255);
		}
	}
	img.save(filename);

}


/**
 * Copy the elements of 'members' to 'membersArray' to facilitate
 * direct access
 */
void REM::copyElementsToArray() {
	list<GroupElement>::iterator xElem;
	int i;

	if(membersArray) free(membersArray);

	membersArray = (GroupElement*) calloc(members.size(), sizeof(GroupElement));

	if (membersArray == NULL) {
		fprintf(stderr, "ERROR: Cannot allocate memory for internal group list\n");
		exit(-1);
	}

	i = 0;
	for (xElem = members.begin(); xElem != members.end(); xElem++) {
		membersArray[i].angle = xElem->angle;
		membersArray[i].x = xElem->x;
		membersArray[i].y = xElem->y;
		membersArray[i].magnitude = xElem->magnitude;
		i++;
	}
}


/**
 * Print a matrix text file of the histogram bin values
 */
void REM::printHistogramMatrix(char * filename) {
	FILE *fp;
	int i, j;
	double totalElems = histogram_sum(hist);

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s to print histogram matrix\n", filename);
		exit(-1);
	}

	for(i = 0; i < binDimension; i++) {
		for(j = 0; j < binDimension; j++) {
			fprintf(fp, "%2.5f ", (double) (histogram_get(hist, i, j) / totalElems));
		}

		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
	fclose(fp);

}


/**
 * Print a histogram text file in histogram format
 */
void REM::printHistogramFile(char * filename) {
	FILE *fp;

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s to print histogram file\n", filename);
		exit(-1);
	}

	histogram_fprintf(fp, hist, "%2.3f", "%2.5f");
}


/**
 * Read a histogram from a file
 */
void REM::readHistogramFile(char * filename) {
	FILE *fp;

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Unable to open file %s to read histogram\n", filename);
		exit(-1);
	}

	histogram_fscanf(fp, hist);

}


/**
 * This class' implementation of log2 for portability
 */
double REM::log2(double value) {
	return (log(value) / log(2.0));
}

