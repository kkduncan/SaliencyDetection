/*
 * RImage.h
 *
 * Header file fo the image representation used in R.E.M. (Relational
 * Entropy-Based Measure of Saliency). Currently, the image types that
 * can be read include PGM, PPM, and JPEG. In a future implementation,
 * it would be able to read and write PNG and BMP images also.
 *
 *  Created on: Dec 25, 2010
 *      Author: kesterduncan
 *
 *	DISCLAIMER:
 *	==============================================================
 *  This code was not designed well. It is not modular and doesn't
 *  obey the Object Oriented Principles. Therefore it may contain
 *  many bugs. Use at your own risk!
 *  ==============================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <math.h>
#include "include/jpeglib.h"


#ifndef RIMAGE_H_
#define RIMAGE_H_

/// Image channels
#define RED 0
#define GREEN 1
#define BLUE 2


/*
 * Types of images this class can read
 */
typedef enum {
	NONE,
	PGM,
	PPM,
	JPEG,
	PNG, 	/* Yet to be implemented */
	BMP		/* Yet to be implemented */

} ImageType;


class RImage {
private:
	int *r;
	int *g;
	int *b;

	int initialize (int **chan, int rows, int cols);
	int readPGM (char file[100]);
	int readPPM (char file[100]);
	int readJPEG (char file[100]);
	int readPNG (char file[100]);
	int savePGM (char file[100]);
	int savePPM (char file[100]);
	int saveJPEG (char file[100], int quality);
	int savePNG (char file[100]);
	int saveText (char file[100]);
	int readIntFromFile(FILE *fp);

	void convertJSAMPLEToImage(const JSAMPLE* jpegImageBuffer);
	void convertImageToJSAMPLE(JSAMPLE* jpegImageBuffer);

public:
	int numRows;
	int numCols;
	int isColorImage;
	ImageType type;

	RImage();
	RImage(const RImage& src);
	virtual ~RImage();
	int initialize (int rows, int cols);
	int initialize (int isColor, int rows, int cols);
	int save (char file[100]);
	int read (char file[100]);

	int inBounds(const int i, const int j);
	int *getChannel(int channel);
	void convertToGrayscale();
	void convertToColor();
	void makeThermalColorImage();
	void clearImage();

	int& operator()(const int i, const int j);
	int& operator()(const int i, const int j, const int rgb);
	int operator()(const int i, const int j) const;
	int operator()(const int i, const int j, const int rgb) const;
	const RImage & operator=(const RImage & src);

	void flushMemory();

};

#endif /* RIMAGE_H_ */
