/*
 * ComputeVideoSaliency.cpp
 *
 *  Created on: Feb 19, 2011
 *      Author: kesterduncan
 *
 *  This file computes the saliency of video image pixels using
 *  an Iterative Kernel Density Estimation method utilizing the
 *  3D gradient.
 *
 *  Copyright (c) 2011 Kester Duncan
 *
 *  $Id: ComputeVideoSaliency$
 *
 *
 *	DISCLAIMER:
 *	==============================================================
 *  This code was not designed well. It is not modular and doesn't
 *  obey the Object Oriented Principles. Therefore it may contain
 *  many bugs. Use at your own risk!
 *  ==============================================================
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "image/RImage.h"
#include "gaussian_filter.h"
#include "KDInfo.h"
#include "ImageSequence.h"
#include "CommonUtils.h"

#define LOG2(X) (log(X) / 0.30103)


// FIXME: remove global variables (poor design)
/** Global variables */
int		numFrames;
int		slabFrames = 5;			// number of frames in a slab of focus
int 	width, height; 			// Dimensions of the video frame
ThreeDGradient 	**gradients; 	// Gradient directions map for each frame in the frame cache
float 	**magnitudes; 			// Edge magnitudes map for each frame in the frame cache
RImage 	map; 					// Saliency map for saliency of video slab
KDInfo	*pixels;				// Stores relevant information for pixels


/** Function Prototypes */
void 			initializeProcess(char* directory, char* fileStem, int startFrameNum, int endFrameNum);
void 			magnitudeQuantization();
void 			processVideo(const int slabStartFrame, const int slabEndFrame, const int percentage, char sequenceDir[255], char fileStem[255]);
void			computeSequenceSaliency(int percentage);
void 			postProcessMap(float sigma);
static void 	allocateMemory();
static void 	freeMemory();
static void 	readGradientsAndMagnitudes(const ImageSequence *sequence);

void 			updateApplicableRegion(Bounds bounds, KSumInfo ksum);
Bounds 			getApplicableBounds(Location3D locs[4], int neighborhood);
KSumInfo 		calculateKernelSum (Location3D locs[4]);
void 			updateEntropy(KDInfo *p);
void 			updateSaliencyMap();
int 			inValidImageBounds(Location3D loc[4]);
float 			calculateGradientAngle(ThreeDGradient a, ThreeDGradient b);


using namespace std;

int main(int argc, char *argv[]) {
	char	sequenceDir[255];
	char 	mapFileName[255];
	char	fileName[255];
	char	fileStem[255];
	char 	origName[255];
	int 	percentage; 		// sampling percentage
	int 	startFrameNum;		// start of frame cache sequence
	int		endFrameNum;		// end of frame cache sequence
	int		slabStartFrame;
	int		slabEndFrame;
	time_t	start, stop;

	if (argc != 5) {
		usageMsg("./ComputeVideoSaliency [SEQUENCE_DIRECTORY] [START_FRAME_NUM] [END_FRAME_NUM] "
				"[SAMPLING_PERCENTAGE]\n" \
				"\t e.g. ./ComputeVideoSaliency video_directory 001 012 25\n");
	}

	/* Mark the start time for all processing */
	STARTTIMER

	strcpy(sequenceDir, argv[1]);
	stripDirectory(sequenceDir, fileName);
	sprintf(fileStem, "%s_", fileName);

	infoMsg("Evaluating the sequence saliency of video in: %s", sequenceDir);

	startFrameNum = atoi(argv[2]);
	endFrameNum = atoi(argv[3]);
	percentage = atoi(argv[4]);

	slabStartFrame = 0;
	slabEndFrame = 0;

	for (int i = startFrameNum, count = 1; i < endFrameNum; i += slabFrames, count++) {
		if (i == startFrameNum) {
			slabStartFrame = startFrameNum;
		} else {
			slabStartFrame += slabFrames;
		}

		if ((slabStartFrame + slabFrames) > endFrameNum) {
			slabEndFrame = endFrameNum;
		} else {
			slabEndFrame = slabStartFrame + slabFrames;
		}

		RImage comp;
		RImage orig;

		infoMsg("[%d] - Evaluating the slab from %d to %d", count, slabStartFrame, slabEndFrame);
		// TODO: remove global numFrames
		numFrames = (slabEndFrame - slabStartFrame);

		if (i == startFrameNum) {
			initializeProcess(sequenceDir, fileStem, slabStartFrame, slabEndFrame);
		}
		processVideo(slabStartFrame, slabEndFrame, percentage, sequenceDir, fileStem);
		postProcessMap(15);

		sprintf(mapFileName, "%s/results/%ssmap_%03d.ppm", sequenceDir, fileStem, count);


		/*
		 * Make composite map
		 */
		int midFrame = ((int)(slabStartFrame + (slabFrames / 2)));
		if (midFrame > slabEndFrame || midFrame < slabStartFrame) {
			midFrame = slabEndFrame;
		}
		sprintf(origName, "%s/%s%04d.pgm", sequenceDir, fileStem, midFrame);
		orig.read(origName);
		comp.isColorImage = 1;
		comp.initialize(map.numRows, (map.numCols * 2) + 1);

		for (int m=0; m < map.numRows; m++) {
			for (int n=0; n < map.numCols; n++) {
				comp(m, n, RED) = orig(m, n, RED);
				comp(m, n, GREEN) = orig(m, n, RED);
				comp(m, n, BLUE) = orig(m, n, RED);
			}
		}

//		map.makeThermalColorImage();

		// Saliency map results
		for (int m=0; m < map.numRows; m++) {
			for (int n=0; n < map.numCols; n++) {
				int col = n + map.numCols;
				double value = map(m, n, RED);
				if (value > 128) {
					comp(m, col, RED) = orig(m, n, RED);
					comp(m, col, GREEN) = orig(m, n, RED);
					comp(m, col, BLUE) = orig(m, n, RED);
				}
			}
		}

		infoMsg("\t...Saving map");
		comp.save(mapFileName);

		comp.flushMemory();
		orig.flushMemory();
		map.clearImage();

	}

	freeMemory();

	/* Mark the end time for all processing and calculate */
	ENDTIMER

	return 0;
}


void processVideo(const int slabStartFrame, const int slabEndFrame, const int percentage,
		char sequenceDir[255], char fileStem[255]) {
	ImageSequence sequence; // Sequence of frames etc

	sequence.initialize(3.0, 2.0, slabStartFrame, slabEndFrame, sequenceDir, fileStem);
	sequence.detectEdges();

	readGradientsAndMagnitudes(&sequence);
	magnitudeQuantization();

	infoMsg("\t...Computing saliency");
	initializeKDInfo(pixels, (width * height));
	computeSequenceSaliency(percentage);
}


/**
 * Initializes global variables etc. in order to start the
 * sequence saliency process.
 */
void initializeProcess(char* directory, char* fileStem, int startFrameNum, int endFrameNum) {
	/* Get the frame width and height */
	char frameName[50];

	sprintf(frameName, "%s/%s%04d.pgm", directory, fileStem, startFrameNum);
	map.read(frameName);
	width = map.numCols;
	height = map.numRows;

	allocateMemory();
}



/**
 * Allocates memory for the global arrays. This function assumes that
 * frameCacheSize has been set to a value greater than or equal to 3.
 */
static void allocateMemory() {
	// FIXME: should read in value of slabFrames
	if (slabFrames >= 3) {
		gradients = (ThreeDGradient**) calloc(slabFrames, sizeof(ThreeDGradient*));
		magnitudes = (float**) calloc(slabFrames, sizeof(float*));

		for (int i=0; i < slabFrames; i++) {
			gradients[i] = (ThreeDGradient*) calloc(width * height, sizeof(ThreeDGradient));
			if (gradients[i] == NULL) {
				fprintf(stderr, "Cannot allocate memory for the gradDirections map.\n");
				exit(1);
			}

			magnitudes[i] = (float*) calloc(width * height, sizeof(float));
			if (magnitudes[i] == NULL) {
				fprintf(stderr, "Cannot allocate memory for the magnitudes map.\n");
				exit(1);
			}

		}

		pixels = (KDInfo *) calloc(width * height, sizeof(KDInfo));
		if (pixels == NULL) {
			fprintf(stderr, "Cannot allocate memory for the pixels map.\n");
			exit(1);
		}

	} else {
		printf("ERROR: The number of slab frames MUST be greater than or equal to 3 frames to "
				"be able to process sequence saliency!\n");
		exit(1);
	}
}

/**
 * Read in the gradient vectors and compute the magnitudes
 */
static void readGradientsAndMagnitudes(const ImageSequence *sequence) {
	/* Get the 3D Gradient vectors */
	// FIXME: Should read in value of slabFrames
	for (int k = 0; k < slabFrames; k++) {
		// TODO: remove global numFrames
		if (numFrames >= slabFrames) {
			ThreeDGradient* grads = sequence->frameGradients[k];
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					int ij = (i * width) + j;
					float delX = grads[ij].delX;
					float delY = grads[ij].delY;
					float delT = grads[ij].delT;

					// TODO - time is weighted more.
					magnitudes[k][ij] = sqrt((delX * delX * 0.30) + (delY * delY * 0.30) + (delT * delT * 0.40));

					gradients[k][ij].delX = delX;
					gradients[k][ij].delY = delY;
					gradients[k][ij].delT = delT;

				}
			}
		}
	}
}

/**
 * Clear memory
 */
static void freeMemory() {
	for (int i=0; i < slabFrames; i++) {
		if (gradients[i] != NULL) {
			free(gradients[i]);
		}
		if (magnitudes[i] != NULL) {
			free(magnitudes[i]);
		}
	}

	if (gradients != NULL) free(gradients);
	if (magnitudes != NULL) free(magnitudes);
	if (pixels != NULL) free(pixels);

	gradients = NULL;
	magnitudes = NULL;
	pixels = NULL;
}


/**
 * Suppresses low magnitudes from the 'magnitudes' array
 * and increases magnitudes closer to the maximum magnitude
 * value
 */
void magnitudeQuantization() {
	if (magnitudes != NULL) {
		double maxMag = -999;

		for (int m = 0; m < slabFrames; m++) {
			for (int k = 0; k < (width * height); k++) {
				if (magnitudes[m][k] > maxMag) {
					maxMag = magnitudes[m][k];
				}
			}
		}

		/* Thresholds */
		double t1 = (double) maxMag / 3.0;
		double t2 = (double) maxMag / 1.5; 	// was 2.0
		double t3 = (double) maxMag / 0.75;	// was 1.5

		for (int m = 0; m < slabFrames; m++) {
			for (int k = 0; k < (width * height); k++) {

				if (magnitudes[m][k] >= 0 && magnitudes[m][k] < t1) {
					magnitudes[m][k] = 0.0005;

				} else if (magnitudes[m][k] >= t1 && magnitudes[m][k] < t2) {
					magnitudes[m][k] = 0.0025; // was 0.0025

				} else if (magnitudes[m][k] >= t2 && magnitudes[m][k] < t3) {
					magnitudes[m][k] = 0.0125; // was 0.0125

				} else if (magnitudes[m][k] >= t3){
					magnitudes[m][k] = 1.0;

				}
			}
		}
	} else {
		fprintf(stderr, "ERROR: The magnitudes map is NULL.\n");
		exit(1);
	}
}


/**
 * Applies Gaussian filtering and image adjustment to a
 * saliency map
 */
void postProcessMap(float sigma) {
	int *filteredImage = NULL;
	float sig2 = 0.6 * sigma;

	/*
	 * First smoothing operation
	 */
	gaussianSmooth((unsigned int*) map.getChannel(RED), map.numRows, map.numCols, sigma, &filteredImage);
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			int ij = (i * width) + j;
			map(i, j, RED) = (int) (filteredImage[ij]);
		}
	}

	free(filteredImage);
	filteredImage = NULL;

	/*
	 * Smoothing another time effectively gets rid of
	 * unsightly rings caused by the Gaussian
	 */
	gaussianSmooth((unsigned int*) map.getChannel(RED), map.numRows, map.numCols, sig2, &filteredImage);
	for (int i=0; i < height; i++) {
		for (int j=0; j < width; j++) {
			int ij = (i * width) + j;
			map(i, j, RED) = (int) (filteredImage[ij]);
		}
	}

	free(filteredImage);

	float maxVal = -999;
	float minVal = 999;

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (map(i, j, RED) > maxVal) {
				maxVal = map(i, j);
			}
			if (map(i, j, RED) < minVal) {
				minVal = map(i, j);
			}
		}
	}

	/* Normalization */
	float sum = 0;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			map(i, j, RED) = (int) ((((map(i, j, RED) - minVal) / (maxVal - minVal)) * 255.0));
			sum += map(i, j, RED);
		}
	}
}


/**
 * Process the saliency of each pixel in an image using an Iterative Approach
 * for calculating the Kernel Density Estimate of the local neighborhood
 * distribution, then calculating the entropy of that distribution for
 * determining the saliency.
 */
void computeSequenceSaliency(int percentage) {
	int pixelSamples;		// Maximum number of pixels sampled randomly
	int mainCount = 0;		// Counter for the main loop
	int halfM = slabFrames / 2;

	pixelSamples = (int)((percentage / 100.0) * width * height * slabFrames * 3);

	/* Ensure that the dimensions are valid before any processing */
	if (width <= 0 || height <= 0) {
		fprintf(stderr, "Error: Image dimensions are invalid.\n");
		exit(1);
	}

	/* Set the seed for randomization */
	srand(time(NULL));

	while (mainCount < pixelSamples) {
		KSumInfo ksum;
		Bounds bounds;
		Location3D locs[4];			// Array to store randomly selected pixel locations

		/*
		 * Randomly select location of first pixel
		 */
		locs[0].y = rand() % height;
		locs[0].x = rand() % width;
		locs[0].t = rand() % slabFrames;

		/*
		 * The other 3 pixels must be selected in the neighborhood
		 * of pixel 1 and within the range of the slab frames
		 */
		locs[1].y = (rand() % slabFrames) + (locs[0].y - halfM);
		locs[1].x = (rand() % slabFrames) + (locs[0].x - halfM);
		locs[1].t = rand() % slabFrames;

		locs[2].y = (rand() % slabFrames) + (locs[0].y - halfM);
		locs[2].x = (rand() % slabFrames) + (locs[0].x - halfM);
		locs[2].t = rand() % slabFrames;

		locs[3].y = (rand() % slabFrames) + (locs[0].y - halfM);
		locs[3].x = (rand() % slabFrames) + (locs[0].x - halfM);
		locs[3].t = rand() % slabFrames;

		ksum = calculateKernelSum(locs);
		bounds = getApplicableBounds(locs, slabFrames);
		updateApplicableRegion(bounds, ksum);
		mainCount++;
	}

	updateSaliencyMap();

}


/**
 * Calculate the intermediate kernel sum from the contribution of the
 * pixels in 'randomPixels'
 */
KSumInfo calculateKernelSum (Location3D locs[4]) {
	KSumInfo ksum;
	double sampleDistance1, sampleAngle1, sampleMag1;
	double sampleDistance2, sampleAngle2, sampleMag2;
	double distanceKernel = 0, angleKernel = 0;
	double binDimension = 10.0; // Conceptual histogram's bin dimension
	double distanceBinWidth = sqrt(pow(width, 2) + pow(height, 2)) / binDimension;
	double angleBinWidth = 3.142 / binDimension;
	double twoPI = 6.283;
	double dNorm = (sqrt(twoPI) * distanceBinWidth);
	double aNorm = (sqrt(twoPI) * angleBinWidth);

	if (inValidImageBounds(locs) == 1) {
		int loc1[2] = { (locs[0].y * width) + locs[0].x,  locs[0].t};
		int loc2[2] = { (locs[1].y * width) + locs[1].x,  locs[1].t};
		int loc3[2] = { (locs[2].y * width) + locs[2].x,  locs[2].t};
		int loc4[2] = { (locs[3].y * width) + locs[3].x,  locs[3].t};

		sampleDistance1 = sqrt(pow((locs[0].y - locs[1].y), 2) + pow((locs[0].x - locs[1].x), 2));
		sampleDistance2 = sqrt(pow((locs[2].y - locs[3].y), 2) + pow((locs[2].x - locs[3].x), 2));

		ThreeDGradient x1 = gradients[loc1[1]][loc1[0]];
		ThreeDGradient x2 = gradients[loc2[1]][loc2[0]];
		ThreeDGradient x3 = gradients[loc3[1]][loc3[0]];
		ThreeDGradient x4 = gradients[loc4[1]][loc4[0]];

		// TODO - experiment with these values
		sampleAngle1 = fmod((fabs(calculateGradientAngle(x1, x2))), 6.283);
		sampleAngle2 = fmod((fabs(calculateGradientAngle(x3, x4))), 6.283);

		sampleMag1 = fabs(magnitudes[loc1[1]][loc1[0]] - magnitudes[loc2[1]][loc2[0]]);
		sampleMag2 = fabs(magnitudes[loc3[1]][loc3[0]] - magnitudes[loc4[1]][loc4[0]]);

		ksum.magSum1 = sampleMag1;
		ksum.magSum2 = sampleMag2;

		distanceKernel = (1/dNorm) * exp((pow(sampleDistance1 - sampleDistance2, 2) / (-2 * pow(distanceBinWidth, 2))));
		angleKernel = (1/aNorm) * exp((pow(sampleAngle1 - sampleAngle2, 2) / (-2 * pow(angleBinWidth, 2))));

		if (sampleMag1 > 0 && sampleMag2 > 0 && distanceKernel > 0 && angleKernel > 0) {
			ksum.kernelSum = (sampleMag1 * sampleMag2 * distanceKernel * angleKernel);
		} else {
			ksum.kernelSum = 0;
		}
	} else {
		ksum.kernelSum = 0;
		ksum.magSum1 = 0;
		ksum.magSum2 = 0;
	}

	return ksum;

}



/**
 * Calculate the gradient angle between two 3D gradient vectors
 */
float calculateGradientAngle(ThreeDGradient a, ThreeDGradient b) {
	float dotProduct = (a.delX * b.delX) + (a.delY * b.delY) + (a.delT * b.delT);
	float lengthA = sqrt((a.delX * a.delX) + (a.delY * a.delY) + (a.delT * a.delT));
	float lengthB = sqrt((b.delX * b.delX) + (b.delY * b.delY) + (b.delT * b.delT));

	if (lengthA > 0 && lengthB > 0) {
		return (acos((double) (dotProduct / (lengthA * lengthB))));
	} else {
		return acos(0.0);
	}
}


/**
 * Determine if pixel locations are within the valid image
 * dimensions
 */
int inValidImageBounds(Location3D loc[4]) {
	int valid = 0;

	int i = 0;
	while (i < 4) {
		if (loc[i].y >= 0 && loc[i].y < height &&
				loc[i].x >= 0 && loc[i].x < width &&
				loc[i].t >= 0 && loc[i].t < slabFrames) {
			valid++;
		}
		i++;
	}

	if (valid == 4) {
		return 1;
	} else {
		return 0;
	}
}


/**
 * Gets the applicable region of the pixels forming two samples
 * that can be updated
 */
Bounds getApplicableBounds(Location3D locs[4], int neighborhood) {
	Bounds bounds;
	int L = neighborhood;		// Local neighborhood dimension
	int maxX = locs[0].x;
	int minX = locs[0].x;
	int maxY = locs[0].y;
	int minY = locs[0].y;
	int xDiff;		// Differences between the max / min x values
	int yDiff;		// Differences between the max / min y values
	int xDisp;		// The disparity in the x axis for forming an M x M neighborhood
	int yDisp;		// The disparity in the y axis for forming an M x M neighborhood

	/* Get the maximum / minimum x and y values of the pixels */
	for (int i = 1; i < 4; i++) {
		if (locs[i].x > maxX) {
			maxX = locs[i].x;
		}

		if (locs[i].y > maxY) {
			maxY = locs[i].y;
		}

		if (locs[i].x < minX) {
			minX = locs[i].x;
		}

		if (locs[i].y < minY) {
			minY = locs[i].y;
		}
	}

	/* Calculate the differences between the max / min values */
	xDiff = maxX - minX;
	yDiff = maxY - minY;

	/* Get the x and y disparity */
	xDisp = (L - xDiff) - 1;
	yDisp = (L - yDiff) - 1;

	/* Calculate the applicable bounds */
	bounds.topLeft.col = minX - xDisp;
	bounds.topLeft.row = minY - yDisp;

	bounds.topRight.col = maxX + xDisp;
	bounds.topRight.row = minY - yDisp;

	bounds.botLeft.col = minX - xDisp;
	bounds.botLeft.row = maxY + yDisp;

	bounds.botRight.col = maxX + xDisp;
	bounds.botRight.row = maxY + yDisp;

	return bounds;
}


/**
 * Used for updating the kernel sums within the applicable region
 * specified
 */
void updateApplicableRegion(Bounds bounds, KSumInfo ksum) {
	int sampleCountLimit = 1000; // change back to 30 - 100

	for (int i = bounds.topLeft.row; i <= bounds.botLeft.row; i++) {
		for (int j = bounds.topLeft.col; j <= bounds.topRight.col; j++) {
			/*
			 * Must be within the image dimensions
			 */
			if (i >= 0 && i < height && j >= 0 && j < width) {
				int ij = (i * width) + j;

				if (pixels[ij].sampleCount < sampleCountLimit) {
					pixels[ij].kernelSum += ksum.kernelSum;
					pixels[ij].magSum1 += ksum.magSum1;
					pixels[ij].magSum2 += ksum.magSum2;
					pixels[ij].sampleCount++;

					/*
					 * Update the pixel entropy every N iterations
					 *
					 * TODO changed from 32
					 */
					if (((pixels[ij].sampleCount + 1) % 64) == 0) {
						updateEntropy(&pixels[ij]);
					}
				}
			}
		}
	}
}


/**
 * This is used to update the entropy of a pixel in the image
 * after N updates. p is a ptr to a KDInfo object with the
 * required information for calculating the entropy.
 */
void updateEntropy(KDInfo *p) {
	if (p->sampleCount > 0) {
		double totalMagSum = p->magSum1 * p->magSum2;
		double estimation = 0;

		/* Special Case Handling */
		if (totalMagSum <= 0) {
			totalMagSum = p->sampleCount;
		}

		if (p->kernelSum < 0 || isnan(p->kernelSum)) {
			p->kernelSum = 0;
		}

		estimation = p->kernelSum / totalMagSum;

		if (estimation <= 0 || isnan(estimation)) {
			p->entropy = -1;
		} else {
			p->entropy = -1 * LOG2(estimation * estimation);
		}
	}
}


/**
 * Updates the saliency map to reflect the current calculations
 * This is called by computeIterativeSaliency after every
 * N iterations.
 */
void updateSaliencyMap() {
	if (pixels != NULL) {
		/* For proper visualization */
		float maxEntropy = -999.0;
		float minEntropy = 999.0;

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				if (pixels[ij].entropy > maxEntropy) {
					maxEntropy = pixels[ij].entropy;
				}
			}
		}

		/*
		 * Areas with 0 magnitude difference are assigned an entropy
		 * value of -1. These areas are now given the value of the
		 * maximum entropy estimated.
		 *
		 * We find the minimum entropy in tandem
		 */
		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				if (pixels[ij].entropy == -1) {
					pixels[ij].entropy = maxEntropy;
				}

				if (pixels[ij].entropy < minEntropy) {
					minEntropy = pixels[ij].entropy;
				}
			}
		}


		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;
				pixels[ij].entropy -= minEntropy;
			}
		}

		/* Also adjust the maximum entropy */
		 maxEntropy -= minEntropy;

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				int ij = (i * width) + j;

				if (maxEntropy > 0) {
					map(i, j) = (int) (255.0 - ((pixels[ij].entropy / maxEntropy) * 255.0));
				} else {
					map(i, j) = 0;
				}
			}
		}

	} else {
		errorMsg("Pixels array is NULL");
	}
}
