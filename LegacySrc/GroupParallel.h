/*
 * GroupParallel.h
 *
 *  Created on: Feb 12, 2010
 *      Author: kesterduncan
 */

#ifndef GROUPPARALLEL_H_
#define GROUPPARALLEL_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

#include <list>
#include "Group.h"

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS
# define __END_DECLS
#endif

__BEGIN_DECLS

typedef struct {
	int id;
	int elemStart;
	int elemEnd;
	Group group;
} ThreadInfo;

void getThreadWorkload(Group& group);
void * buildHistogramParallel(void * params);
void processHistogram();


__END_DECLS

#endif /* GROUPPARALLEL_H_ */
