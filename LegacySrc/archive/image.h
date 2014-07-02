/**
 * Kester Duncan
 * Image.h
 *
 * Handles only ppm and pgm images
 *
 */
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef IMAGE_HPP_
#define IMAGE_HPP_

#define RED 0
#define GREEN 1
#define BLUE 2

typedef enum {
	NONE = 0,
	PGM = 1,
	PPM = 2,
	JPEG = 3,
	PNG  = 4

} ImageType;

class Image {

private:
	int *IM_r, *IM_g, *IM_b;
	int getint(FILE *fp);
	int initialize(int **IM, int nr, int nc);

public:
	int NR, NC; // Number of rows; Number of columns
	ImageType type;

	Image();
	~Image();
	int initialize(int nr, int nc);

	int save(char file[100]);
	int read(char file[100]);

	int inboundp(const int i, const int j);

	int *getChannel(int channel);
	int& operator()(const int i, const int j);
	int& operator()(const int i, const int j, const int rgb);
	int operator()(const int i, const int j) const;
	int operator()(const int i, const int j, const int rgb) const;
	const Image & operator=(const Image & src);
	void flushMemory();

};

#endif

