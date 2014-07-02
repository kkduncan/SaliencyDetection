/* histogram.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 *  Adopted on: Sep 26, 2010
 *  By: K. Duncan
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <stdlib.h>
#include <stdio.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/**
 * Structure
 */
typedef struct {
  size_t nx, ny ;
  double * xrange ;
  double * yrange ;
  double * bin ;
} histogram;

/**
 * Error Related...
 */
enum {
  HISTOGRAM_SUCCESS  = 0,                                                            //!< HISTOGRAM_SUCCESS
  HISTOGRAM_FAILURE  = -1,                                                           //!< HISTOGRAM_FAILURE
  HISTOGRAM_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */                   //!< HISTOGRAM_EDOM
  HISTOGRAM_EINVAL   = 2,   /* invalid argument supplied by user */                  //!< HISTOGRAM_EINVAL
  HISTOGRAM_ESANITY  = 3,   /* sanity check failed - shouldn't happen */             //!< HISTOGRAM_ESANITY
  HISTOGRAM_ENOMEM   = 4,   /* malloc failed */                                      //!< HISTOGRAM_ENOMEM
  HISTOGRAM_EFAILED  = 5
} ;

#define LINEAR_OPT 1

void histogram_error (const char * reason, const char * file, int line, int histogram_errno);

/* HISTOGRAM_ERROR: call the error handler, and return the error code */
#define HISTOGRAM_ERROR(reason, histogram_errno) \
       do { \
       histogram_error (reason, __FILE__, __LINE__, histogram_errno) ; \
       return histogram_errno ; \
       } while (0)

/* HISTOGRAM_ERROR_VAL: call the error handler, and return the given value */
#define HISTOGRAM_ERROR_VAL(reason, histogram_errno, value) \
       do { \
       histogram_error (reason, __FILE__, __LINE__, histogram_errno) ; \
       return value ; \
       } while (0)

/**
 * Histogram function prototypes
 */
histogram * histogram_alloc(const size_t nx, const size_t ny);
histogram * histogram_calloc(const size_t nx, const size_t ny);
histogram * histogram_calloc_uniform(const size_t nx, const size_t ny, const double xmin, const double xmax, const double ymin,	const double ymax);
void histogram_free(histogram * h);

//static int find(const size_t n, const double range[], const double x, size_t * i);
//static int find2d (const size_t nx, const double xrange[], 	const size_t ny, const double yrange[], const double x, const double y, size_t * i, size_t * j);
int histogram_equal_bins_p(const histogram * h1, const histogram * h2);

int histogram_increment(histogram * h, double x, double y);
int histogram_accumulate(histogram * h, double x, double y, double weight);
int histogram_find(const histogram * h, const double x, const double y,	size_t * i, size_t * j);

double histogram_get(const histogram * h, const size_t i, const size_t j);
int histogram_get_xrange(const histogram * h, const size_t i, double * xlower, double * xupper);
int histogram_get_yrange(const histogram * h, const size_t j, double * ylower, double * yupper);

double histogram_xmax(const histogram * h);
double histogram_xmin(const histogram * h);
size_t histogram_nx(const histogram * h);

double histogram_ymax(const histogram * h);
double histogram_ymin(const histogram * h);
size_t histogram_ny(const histogram * h);

void histogram_reset(histogram * h);

histogram * histogram_calloc_range(size_t nx, size_t ny, double *xrange, double *yrange);
int histogram_set_ranges_uniform(histogram * h, double xmin, double xmax, double ymin, double ymax);
int histogram_set_ranges(histogram * h, const double xrange[], size_t xsize, const double yrange[], size_t ysize);

int histogram_memcpy(histogram *dest, const histogram *source);

histogram * histogram_clone(const histogram * source);

double histogram_max_val(const histogram *h);
void histogram_max_bin(const histogram *h, size_t *i, size_t *j);
double histogram_min_val(const histogram *h);
void histogram_min_bin(const histogram *h, size_t *i, size_t *j);

double histogram_sum(const histogram *h);
int histogram_add(histogram *h1, const histogram *h2);

int histogram_fprintf(FILE * stream, const histogram * h, const char * range_format, const char * bin_format);
int histogram_fscanf(FILE * stream, histogram * h);

__END_DECLS


#endif /* HISTOGRAM_H_ */
