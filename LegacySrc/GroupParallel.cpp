/*
 * GroupParallel.c
 *
 *  Created on: Feb 12, 2010
 *      Author: kesterduncan
 */

#include "GroupParallel.h"


const int NUM_THREADS = 2;
int numberOfElements = 0;
pthread_t threads[NUM_THREADS];
ThreadInfo threadInfo[NUM_THREADS];

/**
 * Calculate the workload for each thread
 */
void getThreadWorkload(Group& group) {
	numberOfElements = group.members.size();
	int workload = (int) ceil(numberOfElements / NUM_THREADS);
	int current = 0;

	group.copyElementsToArray();

	for (int i = 0; i < NUM_THREADS && i < numberOfElements; i++) {
		threadInfo[i].id = i + 1;
		threadInfo[i].elemStart = current;
		threadInfo[i].elemEnd = current + workload;
		threadInfo[i].group = group;

		current += workload;
	}
}


/**
 * Builds the histogram in a parallel fashion using threads
 */
void * buildHistogramParallel(void * params) {
	ThreadInfo *tInfo = (ThreadInfo *) params;
	double distance, angle, magnitude;

	for (int i = tInfo->elemStart; i < tInfo->elemEnd; i++) {
		for (int j = tInfo->elemStart; j < tInfo->elemEnd; j++) {
			distance = sqrt(pow((tInfo->group.membersArray[i].x - tInfo->group.membersArray[j].x), 2) +
							pow((tInfo->group.membersArray[i].y - tInfo->group.membersArray[j].y), 2));
			angle = fmod(((2 * M_PI) + (tInfo->group.membersArray[i].angle - tInfo->group.membersArray[j].angle)), (2 * M_PI));
			if (angle > M_PI) {
				angle = (2 * M_PI) - angle;
			}

			magnitude = sqrt(tInfo->group.membersArray[i].magnitude * tInfo->group.membersArray[j].magnitude);
			tInfo->group.updateHistogram(distance, angle, magnitude);
		}
	}

	tInfo = NULL;
	return NULL;
}

/**
 * Issue threads to build up histogram
 */
void processHistogram() {
	int i;

	for (i = 0; i < NUM_THREADS && i < numberOfElements; i++) {
		int returnVal = pthread_create(&threads[i], NULL, buildHistogramParallel, (void *) &threadInfo[i]);
		if (returnVal) {
			printf("\tERROR; return code from pthread_create() is %d\n", returnVal);
			exit(-1);
		}
	}

	for (i = 0; i < NUM_THREADS && i < numberOfElements; i++) {
		pthread_join(threads[i], NULL);
	}

}
