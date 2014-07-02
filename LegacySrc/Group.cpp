/*
 * Copyright (c) 2010 Kester K. Duncan.
 * All rights reserved.
 *
 *  Created on: Nov 27, 2009
 *
 *  Group class implementation:
 *  	This class forms a relational distribution based on the
 *  	pairwise comparisons of its members, then calculates the
 *  	entropy of that distribution.
 */

#include "Group.h"

namespace {
	const char rcs_id[] = "$Id$";
}


Group::Group() {
	entropy = 0;
	membersArray = NULL;
	srand(time(NULL));
}

Group::~Group() {
	gsl_histogram2d_free(histogram);
	if (membersArray != NULL) free(membersArray);
	membersArray = NULL;
}

/**
 * Initialize a group that corresponds to the whole image
 */
void Group::init(int index, int columnWidth, int rowHeight, int dimension) {
	init(index, columnWidth, rowHeight, dimension, 100);
}

/**
 * Initialize a group that corresponds to a group from Perceptual Org. Grouping
 * 	o the groups members list must not be empty before this
 * 	  operation. This init function is ONLY used after the
 *    group members have been set.
 */
void Group::init(int index, int width, int height) {
	list<GroupElement>::iterator xElem, yElem;
	double distance, maximumDist = -999;

	id = index;
	entropy = 0;
	minAngle = 0;
	maxAngle = M_PI;
	minDist = 0;
	binDimension = 7;
	samplePercentage = 1;

	/* For getting the maximum distance possible in group */
	for (xElem = members.begin(); xElem != members.end(); xElem++) {
		for (yElem = members.begin(); yElem != members.end(); yElem++) {
			distance = sqrt(pow((xElem->x - yElem->x), 2) +	pow((xElem->y - yElem->y), 2));
			if (distance > maximumDist) {
				maximumDist = distance;
			}
		}
	}
	/* Sanity check */
	maxDist = (maximumDist > 0) ? maximumDist : sqrt(pow(width, 2) + pow(height, 2));
	initializeHistogram();

	srand(time(0));
}

/**
 * Initialize a group that corresponds to an image
 */
void Group::init(int index, int width, int height, int dimension, int percentage) {
	id = index;
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
void Group::addElement(int x, int y, double angle) {
	GroupElement newMember;
	newMember.x = x;
	newMember.y = y;
	newMember.angle = angle;
	newMember.magnitude = 1;
	members.push_back(newMember);
}

/**
 * Adds an element to the members list
 */
void Group::addElement(int x, int y, double angle, double magnitude) {
	GroupElement newMember;
	newMember.x = x;
	newMember.y = y;
	newMember.angle = angle;

	if (magnitude == 0) magnitude = 1;
	newMember.magnitude = magnitude;
	members.push_back(newMember);
}

/**
 * Initialize the relational histogram for the distance and angle parameters
 */
void Group::initializeHistogram() {
	histogram = gsl_histogram2d_alloc(binDimension, binDimension);
	gsl_histogram2d_set_ranges_uniform(histogram, minDist, maxDist, minAngle, maxAngle);
	elemCount = 0;
}


/**
 * Build up the histogram by comparing each element of the members
 * array with each other
 */
void Group::buildHistogram() {
	list<GroupElement>::iterator xElem, yElem;
	double distance, angle, magnitude;

	for (xElem = members.begin(); xElem != members.end(); xElem++) {
		for (yElem = members.begin(); yElem != members.end(); yElem++) {
			distance = sqrt(pow((xElem->x - yElem->x), 2) + pow((xElem->y - yElem->y), 2));
			angle = fmod(((2 * M_PI) + (xElem->angle - yElem->angle)), (2 * M_PI));

			if (angle > M_PI) {
				angle = (2 * M_PI) - angle;
			}

			magnitude = sqrt(xElem->magnitude * yElem->magnitude);

			if (magnitude != 0) {
				updateHistogram(distance, angle, magnitude);
			} else {
				updateHistogram(distance, angle, 10);
			}
		}
	}
}

/**
 * Build up the histogram by comparing each element of the members
 * array with each other - using sampling of points
 * This function assumes that copyMembersToArray() has been called
 */
void Group::buildSampledHistogram() {
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
		magnitude = fabs(membersArray[randElem1].magnitude * membersArray[randElem2].magnitude);
		if (magnitude != 0) {
			updateHistogram(distance, angle, magnitude);
		} else {
			updateHistogram(distance, angle, 10);
		}
		count++;
	}
}

void Group::clearHistogram() {
	gsl_histogram2d_reset(histogram);
}

/**
 * Update a histogram bin based on distance and angle with equal weights
 * for each group element
 */
void Group::updateHistogram(double distance, double angle) {
	gsl_histogram2d_increment(histogram, distance, angle);
	elemCount++;

}

/**
 * Update a histogram bin based on distance and angle with the weights
 * given by the root of the product of two element magnitudes
 */
void Group::updateHistogram(double distance, double angle, double magnitude) {
	gsl_histogram2d_accumulate(histogram, distance, angle, magnitude);
	elemCount++;

}

/**
 * Calculate the Renyi entropy of this group's relational histogram
 */
double Group::calculateRenyiEntropy() {
	int i, j;
	double totalElems = 0;
	double totalBins = binDimension * binDimension;

	totalElems = gsl_histogram2d_sum(histogram);

	if (entropy != 0) entropy = 0;

	if (totalElems != 0 && gsl_histogram2d_max_val(histogram) != gsl_histogram2d_min_val(histogram)) {
		for(i = 0; i < binDimension; i++) {
			for(j = 0; j < binDimension; j++) {
				double e = 0;
				double p = (double) (gsl_histogram2d_get(histogram, i, j) / totalElems);

				if(p > 0) {
					e = p * p;
				} else {
					e = 0;
				}
				entropy += e;
			}
		}
	} else {
		for(i = 0; i < binDimension; i++) {
			for(j = 0; j < binDimension; j++) {
				double e = 0;
				double p = (double) 1 / (totalBins);

				e = -p * groupLog2(p);
				entropy += e;

				return entropy;
			}
		}
	}

	/* Extra term to penalize small groups */
	entropy = -1 * (groupLog2(entropy));
	entropy += (( (totalBins) - 1) * ((groupLog2(exp(1))) / (2 * totalElems))) / totalBins;

	if(entropy >= 0) {
		entropy /= (double) (groupLog2(totalBins));
	}

	return entropy;
}

/**
 * Calculate the Renyi entropy of this group's relational histogram
 */
double Group::calculateShannonEntropy() {
	int i, j;
	double totalElems = 0;
	double totalBins = binDimension * binDimension;

	totalElems = gsl_histogram2d_sum(histogram);

	if (entropy != 0) entropy = 0;

	if (totalElems != 0 && gsl_histogram2d_max_val(histogram) != gsl_histogram2d_min_val(histogram)) {
		for(i = 0; i < binDimension; i++) {
			for(j = 0; j < binDimension; j++) {
				double e = 0;
				double p = (double) (gsl_histogram2d_get(histogram, i, j) / totalElems);

				if(p > 0) {
					e = -p * groupLog2(p);
				} else {
					e = 0;
				}
				entropy += e;
			}
		}
	} else {
		for(i = 0; i < binDimension; i++) {
			for(j = 0; j < binDimension; j++) {
				double e = 0;
				double p = (double) 1 / (totalBins);

				e = -p * groupLog2(p);
				entropy += e;

				return 0;
			}
		}
	}

	/* Extra term to penalize small groups */
	entropy += (((totalBins) - 1) * ((groupLog2(exp(1))) / (2 * totalElems))) / totalBins;

	if(entropy >= 0) {
		entropy /= (double) (groupLog2(totalBins));
	}

	return entropy;
}

double Group::estimateEntropy() {
	list<GroupElement>::iterator xElem, yElem;
	double distance, angle, sum = 0;
	double part1 = 0, part2 = 0, part3 = 0;
	double estimation = 0, entropyValue = 0;

	for (xElem = members.begin(); xElem != members.end(); xElem++) {
		for (yElem = members.begin(); yElem != members.end(); yElem++) {
			distance = abs(xElem->x - yElem->x) +	abs(xElem->y - yElem->y);
			angle = fmod(((2 * M_PI) + (xElem->angle - yElem->angle)), (2 * M_PI));
			if (angle > M_PI) {
				angle = (2 * M_PI) - angle;
			}

			/* Estimation */
			if (distance != 0 && angle != 0) {
				part1 = (1/(sqrt(2 * M_PI) * distanceBinWidth) * (1/(sqrt(2 * M_PI) * angleBinWidth)));
				part2 = exp(-1 * (pow(distance, 2) / (2 * distanceBinWidth)));
				part3 = exp(-1 * (pow(angle, 2) / (2 * angleBinWidth)));
				sum += (part1 * part2 * part3);
			}

			printf("Part 1 = %2.5f\n", part1);
			printf("Part 2 = %2.5f\n", part2);
			printf("Part 3 = %2.5f\n", part3);

		}
	}

	estimation = (double) sum / (members.size() * members.size());
	entropyValue = -groupLog2(estimation);

	return entropyValue;
}


/**
 * Draw a PGM image of the histogram
 */
void Group::drawHistogram(char * filename) {
	Image img;
	int i, j;
	double maxValue = gsl_histogram2d_max_val(histogram);

	img.initialize(binDimension, binDimension);

	for(i = 0; i < binDimension; i++) {
		for (j = 0; j < binDimension; j++) {
			img(i, j, RED) = (int) ceil(((double) gsl_histogram2d_get(histogram, i, j) / maxValue) * 255);
		}
	}
	img.save(filename);

}

/**
 * Copy the elements of 'members' to 'membersArray'
 */
void Group::copyElementsToArray() {
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
void Group::printHistogramMatrix(char * filename) {
	FILE *fp;
	int i, j;
	double totalElems = gsl_histogram2d_sum(histogram);

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s to print histogram matrix\n", filename);
		exit(-1);
	}

	for(i = 0; i < binDimension; i++) {
		for(j = 0; j < binDimension; j++) {
			fprintf(fp, "%2.5f ", (double) (gsl_histogram2d_get(histogram, i, j) / totalElems));
		}

		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
	fclose(fp);

}

/**
 * Print a histogram text file in gsl_histogram2d format
 */
void Group::printHistogramFile(char * filename) {
	FILE *fp;

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Unable to open file %s to print histogram file\n", filename);
		exit(-1);
	}

	gsl_histogram2d_fprintf(fp, histogram, "%2.3f", "%2.5f");
}

/**
 * Read a histogram from a file
 */
void Group::readHistogramFile(char * filename) {
	FILE *fp;

	if ((fp = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "Unable to open file %s to read histogram\n", filename);
		exit(-1);
	}

	gsl_histogram2d_fscanf(fp, histogram);

}

/**
 * Group class implementation of log2 for portability
 */
double Group::groupLog2(double value) {
	return (log(value) / log(2));
}

