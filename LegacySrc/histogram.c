/*
 * hist.c
 *
 *  Created on: Sep 26, 2010
 *      Author: kesterduncan
 */
#include <stdlib.h>
#include "histogram.h"

#define LINEAR_OPT 1

static int find(const size_t n, const double range[], const double x,
		size_t * i) {
	size_t i_linear, lower, upper, mid;

	if (x < range[0]) {
		return -1;
	}

	if (x >= range[n]) {
		return +1;
	}

	/* optimize for linear case */

#ifdef LINEAR_OPT
	{
		double u = (x - range[0]) / (range[n] - range[0]);
		i_linear = (size_t) (u * n);
	}

	if (x >= range[i_linear] && x < range[i_linear + 1]) {
		*i = i_linear;
		return 0;
	}
#endif

	/* perform binary search */

	upper = n;
	lower = 0;

	while (upper - lower > 1) {
		mid = (upper + lower) / 2;

		if (x >= range[mid]) {
			lower = mid;
		} else {
			upper = mid;
		}
	}

	*i = lower;

	/* sanity check the result */

	if (x < range[lower] || x >= range[lower + 1]) {
		HISTOGRAM_ERROR("x not found in range", HISTOGRAM_ESANITY);
	}

	return 0;
}

static int find2d(const size_t nx, const double xrange[], const size_t ny,
		const double yrange[], const double x, const double y, size_t * i,
		size_t * j) {
	int status = find(nx, xrange, x, i);

	if (status) {
		return status;
	}

	status = find(ny, yrange, y, j);

	if (status) {
		return status;
	}

	return 0;
}

int histogram_equal_bins_p(const histogram * h1, const histogram * h2) {

	if ((h1->nx != h2->nx) || (h1->ny != h2->ny)) {
		return 0;
	}
	{
		size_t i;
		/* init ranges */
		for (i = 0; i <= (h1->nx); i++) {
			if (h1->xrange[i] != h2->xrange[i]) {
				return 0;
			}
		}
		for (i = 0; i <= (h1->ny); i++) {
			if (h1->yrange[i] != h2->yrange[i]) {
				return 0;
			}
		}
	}
	return 1;
}

/**
 * Print histogram errors
 */
void histogram_error(const char * reason, const char * file, int line,
		int histogram_errno) {
	fprintf(stderr, "histogram: %s:%d: ERROR: %s\n", file, line, reason);
	fflush(stderr);
	abort();
}

static void make_uniform(double range[], size_t n, double xmin, double xmax) {
	size_t i;

	for (i = 0; i <= n; i++) {
		double f1 = ((double) (n - i) / (double) n);
		double f2 = ((double) i / (double) n);
		range[i] = f1 * xmin + f2 * xmax;
	}
}

histogram * histogram_alloc(const size_t nx, const size_t ny) {
	histogram *h;

	if (nx == 0) {
		HISTOGRAM_ERROR_VAL("histogram2d length nx must be positive integer", HISTOGRAM_EDOM, 0);
	}

	if (ny == 0) {
		HISTOGRAM_ERROR_VAL("histogram2d length ny must be positive integer",
				HISTOGRAM_EDOM, 0);
	}

	h = (histogram *) malloc(sizeof(histogram));

	if (h == 0) {
		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram2d struct", HISTOGRAM_ENOMEM, 0);
	}

	h->xrange = (double *) malloc((nx + 1) * sizeof(double));

	if (h->xrange == 0) {
		free(h); /* exception in constructor, avoid memory leak */

		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram2d x ranges", HISTOGRAM_ENOMEM, 0);
	}

	h->yrange = (double *) malloc((ny + 1) * sizeof(double));

	if (h->yrange == 0) {
		free(h->xrange);
		free(h); /* exception in constructor, avoid memory leak */

		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram2d y ranges", HISTOGRAM_ENOMEM, 0);
	}

	h->bin = (double *) malloc(nx * ny * sizeof(double));

	if (h->bin == 0) {
		free(h->xrange);
		free(h->yrange);
		free(h); /* exception in constructor, avoid memory leak */

		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram bins", HISTOGRAM_ENOMEM, 0);
	}

	h->nx = nx;
	h->ny = ny;

	return h;
}

histogram * histogram_calloc(const size_t nx, const size_t ny) {
	histogram *h;

	if (nx == 0) {
		HISTOGRAM_ERROR_VAL("histogram2d length nx must be positive integer", HISTOGRAM_EDOM, 0);
	}

	if (ny == 0) {
		HISTOGRAM_ERROR_VAL("histogram2d length ny must be positive integer", HISTOGRAM_EDOM, 0);
	}

	h = (histogram *) malloc(sizeof(histogram));

	if (h == 0) {
		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram2d struct", HISTOGRAM_ENOMEM, 0);
	}

	h->xrange = (double *) malloc((nx + 1) * sizeof(double));

	if (h->xrange == 0) {
		free(h); /* exception in constructor, avoid memory leak */

		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram2d x ranges", HISTOGRAM_ENOMEM, 0);
	}

	h->yrange = (double *) malloc((ny + 1) * sizeof(double));

	if (h->yrange == 0) {
		free(h->xrange);
		free(h); /* exception in constructor, avoid memory leak */

		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram2d y ranges", HISTOGRAM_ENOMEM, 0);
	}

	h->bin = (double *) malloc(nx * ny * sizeof(double));

	if (h->bin == 0) {
		free(h->xrange);
		free(h->yrange);
		free(h); /* exception in constructor, avoid memory leak */

		HISTOGRAM_ERROR_VAL("failed to allocate space for histogram bins", HISTOGRAM_ENOMEM, 0);
	}

	{
		size_t i;

		for (i = 0; i < nx + 1; i++) {
			h->xrange[i] = i;
		}

		for (i = 0; i < ny + 1; i++) {
			h->yrange[i] = i;
		}

		for (i = 0; i < nx * ny; i++) {
			h->bin[i] = 0;
		}
	}

	h->nx = nx;
	h->ny = ny;

	return h;
}

void histogram_free(histogram * h) {
	free(h->xrange);
	free(h->yrange);
	free(h->bin);
	free(h);
}

int histogram_set_ranges_uniform(histogram * h, double xmin, double xmax,
		double ymin, double ymax) {
	size_t i;
	const size_t nx = h->nx, ny = h->ny;

	if (xmin >= xmax) {
		HISTOGRAM_ERROR_VAL("xmin must be less than xmax", HISTOGRAM_EINVAL, 0);
	}

	if (ymin >= ymax) {
		HISTOGRAM_ERROR_VAL("ymin must be less than ymax", HISTOGRAM_EINVAL, 0);
	}

	/* initialize ranges */

	make_uniform(h->xrange, nx, xmin, xmax);
	make_uniform(h->yrange, ny, ymin, ymax);

	/* clear contents */

	for (i = 0; i < nx * ny; i++) {
		h->bin[i] = 0;
	}

	return HISTOGRAM_SUCCESS;
}

int histogram_increment(histogram * h, double x, double y) {
	int status = histogram_accumulate(h, x, y, 1.0);
	return status;
}

int histogram_accumulate(histogram * h, double x, double y, double weight) {
	const size_t nx = h->nx;
	const size_t ny = h->ny;

	size_t i = 0, j = 0;

	int status = find2d(h->nx, h->xrange, h->ny, h->yrange, x, y, &i, &j);

	if (status) {
		return HISTOGRAM_EDOM;
	}

	if (i >= nx) {
		HISTOGRAM_ERROR("index lies outside valid range of 0 .. nx - 1",
				HISTOGRAM_ESANITY);
	}

	if (j >= ny) {
		HISTOGRAM_ERROR("index lies outside valid range of 0 .. ny - 1",
				HISTOGRAM_ESANITY);
	}

	h->bin[i * ny + j] += weight;

	return HISTOGRAM_SUCCESS;
}

double histogram_get(const histogram * h, const size_t i, const size_t j) {
	const size_t nx = h->nx;
	const size_t ny = h->ny;

	if (i >= nx) {
		HISTOGRAM_ERROR_VAL("index i lies outside valid range of 0 .. nx - 1",
				HISTOGRAM_EDOM, 0);
	}

	if (j >= ny) {
		HISTOGRAM_ERROR_VAL("index j lies outside valid range of 0 .. ny - 1",
				HISTOGRAM_EDOM, 0);
	}

	return h->bin[i * ny + j];
}

int histogram_add(histogram * h1, const histogram * h2) {
	size_t i;

	if (!histogram_equal_bins_p(h1, h2)) {
		HISTOGRAM_ERROR("histograms have different binning", HISTOGRAM_EINVAL);
	}

	for (i = 0; i < (h1->nx) * (h1->ny); i++) {
		h1->bin[i] += h2->bin[i];
	}

	return HISTOGRAM_SUCCESS;
}

double histogram_xmax(const histogram * h) {
	const int nx = h->nx;
	return h->xrange[nx];
}

double histogram_xmin(const histogram * h) {
	return h->xrange[0];
}

double histogram_ymax(const histogram * h) {
	const int ny = h->ny;
	return h->yrange[ny];
}

double histogram_ymin(const histogram * h) {
	return h->yrange[0];
}

size_t histogram_nx(const histogram * h) {
	return h->nx;
}

size_t histogram_ny(const histogram * h) {
	return h->ny;
}

void histogram_reset(histogram * h) {
	size_t i;
	const size_t nx = h->nx;
	const size_t ny = h->ny;

	for (i = 0; i < nx * ny; i++) {
		h->bin[i] = 0;
	}
}

double histogram_sum(const histogram * h) {
	const size_t n = h->nx * h->ny;
	double sum = 0;
	size_t i = 0;

	while (i < n)
		sum += h->bin[i++];

	return sum;
}

/*
 * Return the maximum contents value of a 2D histogram
 */
double histogram_max_val(const histogram * h) {
	const size_t nx = h->nx;
	const size_t ny = h->ny;
	size_t i;
	double max = h->bin[0 * ny + 0];

	for (i = 0; i < nx * ny; i++) {
		if (h->bin[i] > max) {
			max = h->bin[i];
		}
	}
	return max;
}

/*
 * Find the bin index for maximum value of a 2D histogram
 */

void histogram_max_bin(const histogram * h, size_t * imax_out,
		size_t * jmax_out) {
	const size_t nx = h->nx;
	const size_t ny = h->ny;
	size_t imax = 0, jmax = 0, i, j;
	double max = h->bin[0 * ny + 0];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			double x = h->bin[i * ny + j];

			if (x > max) {
				max = x;
				imax = i;
				jmax = j;
			}
		}
	}

	*imax_out = imax;
	*jmax_out = jmax;
}

/*
 * Return the minimum contents value of a 2D histogram
 */

double histogram_min_val(const histogram * h) {
	const size_t nx = h->nx;
	const size_t ny = h->ny;
	size_t i;
	double min = h->bin[0 * ny + 0];

	for (i = 0; i < nx * ny; i++) {
		if (h->bin[i] < min) {
			min = h->bin[i];
		}
	}

	return min;
}
/*
 * Find the bin index for minimum value of a 2D histogram
 */
void histogram_min_bin(const histogram * h, size_t * imin_out,
		size_t * jmin_out) {
	const size_t nx = h->nx;
	const size_t ny = h->ny;
	size_t imin = 0, jmin = 0, i, j;
	double min = h->bin[0 * ny + 0];

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			double x = h->bin[i * ny + j];

			if (x < min) {
				min = x;
				imin = i;
				jmin = j;
			}
		}
	}

	*imin_out = imin;
	*jmin_out = jmin;
}

int histogram_fprintf(FILE * stream, const histogram * h,
		const char *range_format, const char *bin_format) {
	size_t i, j;
	const size_t nx = h->nx;
	const size_t ny = h->ny;
	int status;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			status = fprintf(stream, range_format, h->xrange[i]);

			if (status < 0) {
				HISTOGRAM_ERROR("fprintf failed", HISTOGRAM_EFAILED);
			}

			status = putc(' ', stream);

			if (status == EOF) {
				HISTOGRAM_ERROR("putc failed", HISTOGRAM_EFAILED);
			}

			status = fprintf(stream, range_format, h->xrange[i + 1]);

			if (status < 0) {
				HISTOGRAM_ERROR("fprintf failed", HISTOGRAM_EFAILED);
			}

			status = putc(' ', stream);

			if (status == EOF) {
				HISTOGRAM_ERROR("putc failed", HISTOGRAM_EFAILED);
			}

			status = fprintf(stream, range_format, h->yrange[j]);

			if (status < 0) {
				HISTOGRAM_ERROR("fprintf failed", HISTOGRAM_EFAILED);
			}

			status = putc(' ', stream);

			if (status == EOF) {
				HISTOGRAM_ERROR("putc failed", HISTOGRAM_EFAILED);
			}

			status = fprintf(stream, range_format, h->yrange[j + 1]);

			if (status < 0) {
				HISTOGRAM_ERROR("fprintf failed", HISTOGRAM_EFAILED);
			}

			status = putc(' ', stream);

			if (status == EOF) {
				HISTOGRAM_ERROR("putc failed", HISTOGRAM_EFAILED);
			}

			status = fprintf(stream, bin_format, h->bin[i * ny + j]);

			if (status < 0) {
				HISTOGRAM_ERROR("fprintf failed", HISTOGRAM_EFAILED);
			}

			status = putc('\n', stream);

			if (status == EOF) {
				HISTOGRAM_ERROR("putc failed", HISTOGRAM_EFAILED);
			}
		}
		status = putc('\n', stream);

		if (status == EOF) {
			HISTOGRAM_ERROR("putc failed", HISTOGRAM_EFAILED);
		}
	}

	return HISTOGRAM_SUCCESS;
}

int histogram_fscanf(FILE * stream, histogram * h) {
	size_t i, j;
	const size_t nx = h->nx;
	const size_t ny = h->ny;
	double xupper, yupper;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			int status = fscanf(stream, "%lg %lg %lg %lg %lg", h->xrange + i,
					&xupper, h->yrange + j, &yupper, h->bin + i * ny + j);

			if (status != 5) {
				HISTOGRAM_ERROR("fscanf failed", HISTOGRAM_EFAILED);
			}
		}
		h->yrange[ny] = yupper;
	}

	h->xrange[nx] = xupper;

	return HISTOGRAM_SUCCESS;
}
