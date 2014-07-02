/*
 * CalculatePerformance.cpp
 *
 *  Created on: Apr 4, 2010
 *      Author: kesterduncan
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "image/RImage.h"

using namespace std;

double calculateCorrelation(const RImage& hmap, const RImage& smap) {
	double humanAvg = 0, salAvg = 0;
	double humanSum = 0, salSum = 0;
	double hdenom = 0, sdenom = 0;
	double numerator = 0;
	int count = 0;

	for (int i = 0; i < hmap.numRows; i++) {
		for (int j = 0; j < hmap.numCols; j++) {
			humanSum += hmap(i, j);
			salSum += smap(i, j);
			count++;
		}
	}

	humanAvg = (float) humanSum / count;
	salAvg = (float) salSum / count;

	for (int i = 0; i < hmap.numRows; i++) {
		for (int j = 0; j < hmap.numCols; j++) {
			numerator += ((hmap(i, j) - humanAvg) * (smap(i, j) - salAvg));
			hdenom += pow(hmap(i, j) - humanAvg, 2);
			sdenom += pow(smap(i, j) - salAvg, 2);
		}
	}

	return numerator / sqrt(hdenom * sdenom);

}

void getComparisonMap(const RImage& hmap, const RImage& smap, RImage& cmap) {
	for (int i = 0; i < hmap.numRows; i++) {
		for (int j = 0; j < hmap.numCols; j++) {
			cmap(i, j, RED) = 0;
			cmap(i, j, GREEN) = hmap(i, j);
			cmap(i, j, BLUE) = smap(i, j);
		}
	}
}

void histogramEqualization(RImage& img) {
	int hist[256];
	double s_hist_eq[256] = { 0.0 }, sum_of_hist[256] = { 0.0 };
	int i, j, k, n;

	n = img.numRows * img.numCols;

	for (i = 0; i < 256; i++) {
		hist[i] = 0;
	}

	for (i = 0; i < img.numRows; i++) {
		for (j = 0; j < img.numCols; j++) {
			hist[img(i, j)]++;
		}
	}

	for (i = 0; i < 256; i++) {
		s_hist_eq[i] = (double) hist[i] / (double) n;
	}

	sum_of_hist[0] = s_hist_eq[0];

	for (i = 1; i < 256; i++) // cdf of image
	{
		sum_of_hist[i] = sum_of_hist[i - 1] + s_hist_eq[i];
	}

	for (i = 0; i < img.numRows; i++) {
		for (j = 0; j < img.numCols; j++) {
			k = img(i, j);
			img(i, j) = (int) ((sum_of_hist[k] * 255) + 0.5);
		}
	}

}

int main(int argc, char *argv[]) {

	if (argc != 4) {
		fprintf(stderr, "Error: Incorrect arguments provided\n");
		exit(0);
	}

	char humanMapName[50];
	char saliencyMapName[50];
	char compImageName[50];
	char perfFilename[50] = "performance.csv";
	RImage humanMap;
	RImage saliencyMap;
	RImage compMap;
	FILE *perfFile;

	strcpy(humanMapName, argv[1]);
	strcpy(saliencyMapName, argv[2]);
	strcpy(compImageName, argv[3]);
//	sprintf(compImageName, "%s_cmap.ppm", saliencyMapName);

	if ((perfFile = fopen(perfFilename, "a+")) == NULL) {
		fprintf(stderr, "Cannot open file, %s for reading.\n", perfFilename);
		return (0);
	}

	humanMap.read(humanMapName);
	saliencyMap.read(saliencyMapName);
	compMap.read(compImageName);
//	compMap.initialize(saliencyMap.NR, saliencyMap.NC);

	double iLabVal = calculateCorrelation(humanMap, saliencyMap);
	if (iLabVal < 0) {
		iLabVal = calculateCorrelation(saliencyMap, humanMap);
		if (iLabVal < 0) iLabVal = -1 * iLabVal;
	}

	double resultVal = calculateCorrelation(humanMap, compMap);
	if (resultVal < 0) {
		resultVal = calculateCorrelation(compMap, humanMap);
		if (resultVal < 0) resultVal = -1 * resultVal;
	}

	int winner = -1;
	if (resultVal > iLabVal) {
		winner = 1;
	} else {
		winner = 0;
	}

	fprintf(perfFile, "%s, %s, %s, %2.5f, %2.5f, %d\n", humanMapName, saliencyMapName,
			compImageName, iLabVal, resultVal, winner);

	fclose(perfFile);

//	getComparisonMap(humanMap, saliencyMap, compMap);
//	compMap.save(compImageName);
//	histogramEqualization(saliencyMap);
//	saliencyMap.save(compImageName);

	return 0;
}
