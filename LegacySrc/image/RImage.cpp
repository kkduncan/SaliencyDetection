/*
 * RImage.cpp
 *
 * Author: Kester Duncan
 *
 */

#include "RImage.h"

/*
 * Constructor.
 * - initializes instance variables
 */
RImage::RImage() {
	numRows = numCols = 0;
	r = g = b = NULL;
	type = NONE;
	isColorImage = 0;
}


/*
 * Copy constructor
 */
RImage::RImage(const RImage& src) {
	numRows = src.numRows;
	numCols = src.numCols;
	type = src.type;
	initialize(src.numRows, src.numCols);

	for (int i = 0; i < src.numRows; i++) {
		for (int j = 0; j < src.numCols; j++) {
			int ij = (i * src.numCols) + j;
			*(r + ij) = src(i, j, RED);

			if (src.isColorImage) {
				*(g + ij) = src(i, j, GREEN);
				*(b + ij) = src(i, j, BLUE);
			}
		}
	}
}


RImage::~RImage() {
	flushMemory();
}


/*
 * Copy constructor.
 * - performs a deep copy of another RImage object
 */
const RImage & RImage::operator=(const RImage &src) {
	this->numRows = src.numRows;
	this->numCols = src.numCols;
	this->type = src.type;
	this->isColorImage = src.isColorImage;
	this->initialize(src.numRows, src.numCols);

	for (int i = 0; i < src.numRows; i++) {
		for (int j = 0; j < src.numCols; j++) {
			int ij = (i * src.numCols) + j;
			*(r + ij) = src(i, j, RED);

			if (src.isColorImage ==  1) {
				*(g + ij) = src(i, j, GREEN);
				*(b + ij) = src(i, j, BLUE);
			}
		}
	}

	return (*this);
}


/*
 * This function returns access to the value of the red channel
 * at the specified location.
 */
int& RImage::operator()(const int i, const int j) {
	return *(r + i * numCols + j);
}


int RImage::operator()(const int i, const int j) const {
	return *(r + i * numCols + j);
}


/*
 * This function returns access to the value of the specified
 * channel at the specified location.
 */
int& RImage::operator()(const int i, const int j, const int rgb) {
	switch (rgb) {
	case RED:
		return *(r + i * numCols + j);
		break;
	case GREEN:
		return *(g + i * numCols + j);
		break;
	case BLUE:
		return *(b + i * numCols + j);
		break;
	}

	return *(r + i * numCols + j);
}


int RImage::operator()(const int i, const int j, const int rgb) const {
	switch (rgb) {
	case RED:
		return *(r + i * numCols + j);
		break;
	case GREEN:
		return *(g + i * numCols + j);
		break;
	case BLUE:
		return *(b + i * numCols + j);
		break;
	}

	return *(r + i * numCols + j);
}


/*
 * Initialize.
 * - initializes an image channel
 */
int RImage::initialize(int **chan, int rows, int cols) {
	if (((*chan) == NULL) || (((*chan) != NULL) && (rows * cols >= numCols * numRows))) {
		if ((*chan) != NULL) {
			free(*chan);
		}

		if (((*chan) = (int*) malloc(sizeof(int) * rows * cols)) == NULL) {
			fprintf(stderr, "RImage Error: Could not allocate space(%d %d) at line %d in file %s\n",
					rows, cols, __LINE__, __FILE__);
			return(1);
		}
	}

	numRows = rows;
	numCols = cols;

	for (int i = 0; i < rows * cols; i++) {
		(*chan)[i] = 0;
	}

	return 0;
}


/*
 * Initialize
 * - initializes image channels.
 */
int RImage::initialize(int rows, int cols) {
	initialize(&(r), rows, cols);

	if (isColorImage == 1) {
		initialize(&(g), rows, cols);
		initialize(&(b), rows, cols);
	}

	return 0;
}


/*
 * Free allocated memory
 */
void RImage::flushMemory() {
	if (r != NULL) {
		free(r);
	}

	if (isColorImage == 1) {
		if (g != NULL) {
			free(g);
		}

		if (b != NULL) {
			free(b);
		}
	}

	numRows = numCols = 0;
	r = g = b = NULL;
}

/*
 * Clear the contents of an image
 */
void RImage::clearImage() {
	for (int i = 0; i < numRows * numCols; i++) {
		r[i] = 0;
		if (isColorImage) {
			g[i] = 0;
			b[i] = 0;
		}
	}
}


/*
 * Read an image
 */
int RImage::read(char file[100]) {
	if (strstr(file, ".ppm") != NULL) {
		type = PPM;
		readPPM(file);

	} else if (strstr(file, ".pgm") != NULL){
		type = PGM;
		readPGM(file);

	} else if (strstr(file, ".jpeg") != NULL ||
			strstr(file, ".jpg") != NULL) {
		type = JPEG;
		readJPEG(file);

	} else if (strstr(file, ".png") != NULL) {
		fprintf(stderr, "ERROR: Format not supported as yet.\n");
		type = PNG;
		readPNG(file);

	} else {
		fprintf(stderr, "ERROR: Cannot read image %s.\n", file);
		exit(1);
	}

	return 0;
}


/*
 * Read a PGM image
 */
int RImage::readPGM(char file[100]) {
	FILE *fp;
	char str[50];
	int i;

	type = PGM;
	isColorImage = 0;

	if ((fp = fopen(file, "r")) == NULL) {
		fprintf(stderr, "IMAGE ERROR: Cannot open image file %s\n", file);
		exit(1);
	}

	fscanf(fp, "%s", str);

	/* Check for PGM File */
	if ((strcmp(str, "P5") != 0)) {
		fprintf(stderr, "IMAGE ERROR: Image file %s not in PGM format\n", file);
		return (1);
	}

	// Get column and row values
	numCols = readIntFromFile(fp);
	numRows = readIntFromFile(fp);

	fscanf(fp, "%d", &i);
	(getc(fp));

	initialize(&this->r, numRows, numCols);

	for (i = 0; i < numRows * numCols; i++) {
		*(r + i) = (unsigned char) (getc(fp));
	}

	fclose(fp);

	return (0);
}


/*
 * Read a PPM image
 */
int RImage::readPPM(char file[100]) {
	FILE *fp;
	char str[50];
	int i;

	type = PPM;
	isColorImage = 1;

	if ((fp = fopen(file, "r")) == NULL) {
		fprintf(stderr, "Cannot open image file %s\n", file);
		exit(1);
	}

	fscanf(fp, "%s", str);

	/* Check for PGM File */
	if ((strcmp(str, "P6") != 0)) {
		fprintf(stderr, "Image file %s not in PPM format\n", file);
		return (1);
	}

	// Get column and row values
	numCols = readIntFromFile(fp);
	numRows = readIntFromFile(fp);

	fscanf(fp, "%d", &i);
	(getc(fp));

	initialize(numRows, numCols);

	for (i = 0; i < numRows * numCols; i++) {
		*(r + i) = (unsigned char) (getc(fp));
		*(g + i) = (unsigned char) (getc(fp));
		*(b + i) = (unsigned char) (getc(fp));
	}

	fclose(fp);
	return (0);
}


/**
 * Items relevant to readJPEG and saveJPEG
 */
struct my_error_mgr {
	struct jpeg_error_mgr pub;
	jmp_buf setjmp_buffer;
};

typedef struct my_error_mgr * my_error_ptr;

static void my_error_exit (j_common_ptr cinfo) {
	my_error_ptr myerr = (my_error_ptr) cinfo->err;
	(*cinfo->err->output_message) (cinfo);
	longjmp(myerr->setjmp_buffer, 1);
}

/*
 * Read a JPEG image
 */
int RImage::readJPEG(char file[100]) {
	int jpeg_height;
	int jpeg_width;
	int numComponents = 3; //TODO: changed from 1

	/* This struct contains the JPEG decompression parameters and pointers to
	 * working space (which is allocated as needed by the JPEG library).
	 */
	struct jpeg_decompress_struct cinfo;
	/* We use our private extension JPEG error handler.
	 * Note that this struct must live as long as the main JPEG parameter
	 * struct, to avoid dangling-pointer problems.
	 */
	struct my_error_mgr jerr;
	/* More stuff */
	FILE * infile; /* source file */
	JSAMPARRAY buffer; /* Output row buffer */
	int row_stride; /* physical row width in output buffer */

	/* In this example we want to open the input file before doing anything else,
	 * so that the setjmp() error recovery below can assume the file is open.
	 * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
	 * requires it in order to read binary files.
	 */
	if ((infile = fopen(file, "rb")) == NULL) {
		fprintf(stderr, "can't open %s\n", file);
		return 0;
	}

	/* Step 1: allocate and initialize JPEG decompression object */

	/* We set up the normal JPEG error routines, then override error_exit. */
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;
	/* Establish the setjmp return context for my_error_exit to use. */
	if (setjmp(jerr.setjmp_buffer)) {
		/* If we get here, the JPEG code has signaled an error.
		 * We need to clean up the JPEG object, close the input file, and return.
		 */
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return 0;
	}
	/* Now we can initialize the JPEG decompression object. */
	jpeg_create_decompress(&cinfo);

	/* Step 2: specify data source (eg, a file) */
	jpeg_stdio_src(&cinfo, infile);

	/* Step 3: read file parameters with jpeg_read_header() */
	(void) jpeg_read_header(&cinfo, TRUE);
	/* We can ignore the return value from jpeg_read_header since
	 *   (a) suspension is not possible with the stdio data source, and
	 *   (b) we passed TRUE to reject a tables-only JPEG file as an error.
	 * See libjpeg.txt for more info.
	 */

	/* Step 4: set parameters for decompression */

	/* In this example, we don't need to change any of the defaults set by
	 * jpeg_read_header(), so we do nothing here.
	 */

	/* Step 5: Start decompressor */

	(void) jpeg_start_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	 * with the stdio data source.
	 */

	/* We may need to do some setup of our own at this point before reading
	 * the data.  After jpeg_start_decompress() we have the correct scaled
	 * output image dimensions available, as well as the output colormap
	 * if we asked for color quantization.
	 * In this example, we need to make an output work buffer of the right size.
	 */
	/* JSAMPLEs per row in output buffer */
	jpeg_height = cinfo.output_height;
	jpeg_width = cinfo.output_width;

	if (cinfo.jpeg_color_space == JCS_RGB || cinfo.jpeg_color_space == JCS_YCbCr) {
		numComponents = 3;
	}
	if (numComponents != cinfo.output_components) {
		printf("ERROR: Number of color components do not match.\n");
		exit(1);
	}

	JSAMPLE* jpegImageBuffer = (JSAMPLE*) malloc(sizeof(JSAMPLE) * jpeg_height * jpeg_width * numComponents);


	row_stride = cinfo.output_width * cinfo.output_components;
	/* Make a one-row-high sample array that will go away when done with image */
	buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

	/* Step 6: while (scan lines remain to be read) */
	/*           jpeg_read_scanlines(...); */

	/* Here we use the library's state variable cinfo.output_scanline as the
	 * loop counter, so that we don't have to keep track ourselves.
	 */
	int row_ctr = 0;
	while (cinfo.output_scanline < cinfo.output_height) {
		/* jpeg_read_scanlines expects an array of pointers to scanlines.
		 * Here the array is only one element long, but you could ask for
		 * more than one scanline at a time if that's more convenient.
		 */
		(void) jpeg_read_scanlines(&cinfo, buffer, 1);
		/* Assume put_scanline_someplace wants a pointer and sample count. */
		int i;
		for (i = 0; i < jpeg_width * numComponents; i++) {
			int pos = (row_ctr * (jpeg_width * numComponents)) + i;
			jpegImageBuffer[pos] = buffer[0][i];
		}
		row_ctr++;
	}

	/* Step 7: Finish decompression */
	(void) jpeg_finish_decompress(&cinfo);
	/* We can ignore the return value since suspension is not possible
	 * with the stdio data source.
	 */

	/* Step 8: Release JPEG decompression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_decompress(&cinfo);

	/* After finish_decompress, we can close the input file.
	 * Here we postpone it until after no more JPEG errors are possible,
	 * so as to simplify the setjmp error logic above.  (Actually, I don't
	 * think that jpeg_destroy can do an error exit, but why assume anything...)
	 */
	fclose(infile);

	/*
	 * At this point you may want to check to see whether any corrupt-data
	 * warnings occurred (test whether jerr.pub.num_warnings is nonzero).
	 */

	/*
	 * Now, RImage customized processing
	 */
	this->numRows = jpeg_height;
	this->numCols = jpeg_width;

	if (cinfo.jpeg_color_space == JCS_GRAYSCALE) {
		isColorImage = 0;
		initialize(&(this->r), this->numRows, this->numCols);

	} else if (cinfo.jpeg_color_space == JCS_RGB || cinfo.jpeg_color_space == JCS_YCbCr) {
		isColorImage = 1;
		initialize(this->numRows, this->numCols);

	} else {
		fprintf(stderr, "ERROR: Unknown or unsupported image type being read.\n");
		exit(1);
	}

	convertJSAMPLEToImage(jpegImageBuffer);
	this->type = JPEG;

	free(jpegImageBuffer);

	/* And we're done! */
	return 1;

}


/*
 * Read a PNG image .. to be implemented
 */
int RImage::readPNG(char file[100]) {
	/**
	 * TO BE IMPLEMENTED
	 */
	return 0;
}


/*
 * Save an image to one of the formats implemented by this class
 */
int RImage::save(char file[100]) {
	if (strstr(file, ".ppm") != NULL) {
		savePPM(file);

	} else if (strstr(file, ".pgm") != NULL){
		savePGM(file);

	} else if (strstr(file, ".jpeg") != NULL ||
			strstr(file, ".jpg") != NULL) {
		saveJPEG(file, 100);

	} else if (strstr(file, ".png") != NULL) {
		fprintf(stderr, "ERROR: Format not supported as yet.\n");
		savePNG(file);

	} else if (strstr(file, ".txt") != NULL) {
		saveText(file);

	} else {
		fprintf(stderr, "ERROR: Cannot save image %s. Unsupported format.\n", file);
		exit(1);
	}

}


/*
 * Save the image as a PGM image
 */
int RImage::savePGM(char file[100]) {
	FILE *outfile1;
	int i, j;

	fprintf(stderr, "Saving PGM image in %s...\n", file);
	if ((outfile1 = fopen(file, "w")) == NULL) {
		fprintf(stderr, "Cannot open file, %s for writing.\n", file);
		return (0);
	}

	fprintf(outfile1, "P5\n%d %d\n255\n", numCols, numRows);

	for (i = 0; i < numRows; i++) {
		for (j = 0; j < numCols; j++) {
			int loc = (i * numCols) + j;
			putc((unsigned char) (((*(r + loc)) > 255) ? 255
					: (*(r + loc))), outfile1);
		}
	}

	fclose(outfile1);

	return 0;
}


/*
 * Save the image as a PPM image
 */
int RImage::savePPM(char file[100]) {
	FILE *outfile1;
	int i, j;

	fprintf(stderr, "Saving PPM image in %s...\n", file);
	if ((outfile1 = fopen(file, "w")) == NULL) {
		fprintf(stderr, "Cannot open file, %s for writing.\n", file);
		return (0);
	}

	fprintf(outfile1, "P6\n%d %d\n255\n", numCols, numRows);

	for (i = 0; i < numRows; i++) {
		for (j = 0; j < numCols; j++) {
			int loc = (i * numCols) + j;
			putc((unsigned char) (((*(r + loc)) > 255) ? 255
					: (*(r + loc))), outfile1);
			putc((unsigned char) (((*(g + loc)) > 255) ? 255
					: (*(g + loc))), outfile1);
			putc((unsigned char) (((*(b + loc)) > 255) ? 255
					: (*(b + loc))), outfile1);
		}
	}

	fclose(outfile1);

	return 0;
}


/*
 * Save the image as a JPEG
 */
int RImage::saveJPEG (char * filename, int quality) {
	int numChannels = 1;

	if (this->isColorImage == 1|| (g != NULL && b != NULL)) {
		numChannels = 3;
	}

	JSAMPLE* jpegImageBuffer = (JSAMPLE*) malloc(sizeof(JSAMPLE) * this->numRows * this->numCols * numChannels);
	convertImageToJSAMPLE(jpegImageBuffer);

	/*
	 * This struct contains the JPEG compression parameters and pointers to
	 * working space (which is allocated as needed by the JPEG library).
	 * It is possible to have several such structures, representing multiple
	 * compression/decompression processes, in existence at once.  We refer
	 * to any one struct (and its associated working data) as a "JPEG object".
	 */
	struct jpeg_compress_struct cinfo;
	/*
	 * This struct represents a JPEG error handler.  It is declared separately
	 * because applications often want to supply a specialized error handler
	 * (see the second half of this file for an example).  But here we just
	 * take the easy way out and use the standard error handler, which will
	 * print a message on stderr and call exit() if compression fails.
	 * Note that this struct must live as long as the main JPEG parameter
	 * struct, to avoid dangling-pointer problems.
	 */
	struct jpeg_error_mgr jerr;
	/* More stuff */
	FILE * outfile; /* target file */
	JSAMPROW row_pointer[1]; /* pointer to JSAMPLE row[s] */
	int row_stride; /* physical row width in image buffer */

	/* Step 1: allocate and initialize JPEG compression object */

	/* We have to set up the error handler first, in case the initialization
	 * step fails.  (Unlikely, but it could happen if you are out of memory.)
	 * This routine fills in the contents of struct jerr, and returns jerr's
	 * address which we place into the link field in cinfo.
	 */
	cinfo.err = jpeg_std_error(&jerr);
	/* Now we can initialize the JPEG compression object. */
	jpeg_create_compress(&cinfo);

	/* Step 2: specify data destination (eg, a file) */
	/* Note: steps 2 and 3 can be done in either order. */

	/* Here we use the library-supplied code to send compressed data to a
	 * stdio stream.  You can also write your own code to do something else.
	 * VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
	 * requires it in order to write binary files.
	 */
	if ((outfile = fopen(filename, "wb")) == NULL) {
		fprintf(stderr, "can't open %s\n", filename);
		exit(1);
	}
	jpeg_stdio_dest(&cinfo, outfile);

	/* Step 3: set parameters for compression */

	/* First we supply a description of the input image.
	 * Four fields of the cinfo struct must be filled in:
	 */
	cinfo.image_width = this->numCols; /* image width and height, in pixels */
	cinfo.image_height = this->numRows;

	if (this->isColorImage == 1) {
		cinfo.input_components = 3;
		cinfo.in_color_space = JCS_RGB; /* colorspace of input image */
	} else {
		cinfo.input_components = 1; /* # of color components per pixel */
		cinfo.in_color_space = JCS_GRAYSCALE; /* colorspace of input image */
	}

	/* Now use the library's routine to set default compression parameters.
	 * (You must set at least cinfo.in_color_space before calling this,
	 * since the defaults depend on the source color space.)
	 */
	jpeg_set_defaults(&cinfo);
	/* Now you can set any non-default parameters you wish to.
	 * Here we just illustrate the use of quality (quantization table) scaling:
	 */
	jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

	/* Step 4: Start compressor */

	/* TRUE ensures that we will write a complete interchange-JPEG file.
	 * Pass TRUE unless you are very sure of what you're doing.
	 */
	jpeg_start_compress(&cinfo, TRUE);

	/* Step 5: while (scan lines remain to be written) */
	/*           jpeg_write_scanlines(...); */

	/* Here we use the library's state variable cinfo.next_scanline as the
	 * loop counter, so that we don't have to keep track ourselves.
	 * To keep things simple, we pass one scanline per call; you can pass
	 * more if you wish, though.
	 */
	if (this->isColorImage == 1) {
		row_stride = this->numCols * 3; /* JSAMPLEs per row in jpegImageBuffer */
	} else {
		row_stride = this->numCols * 1;
	}

	while (cinfo.next_scanline < cinfo.image_height) {
		/* jpeg_write_scanlines expects an array of pointers to scanlines.
		 * Here the array is only one element long, but you could pass
		 * more than one scanline at a time if that's more convenient.
		 */
		row_pointer[0] = &jpegImageBuffer[cinfo.next_scanline * row_stride];
		(void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	/* Step 6: Finish compression */

	jpeg_finish_compress(&cinfo);
	/* After finish_compress, we can close the output file. */
	fclose(outfile);

	/* Step 7: release JPEG compression object */

	/* This is an important step since it will release a good deal of memory. */
	jpeg_destroy_compress(&cinfo);

	/* And we're done! */
	free(jpegImageBuffer);

}


/*
 * Save the image as a Text File Matrix
 */
int RImage::saveText(char file[100]) {
	FILE *outfile1;
	int i, j;

	if (this->type != PGM) {
		fprintf(stderr, "ERROR: This image must be grayscale in order to print to a text file.\n");
		exit(1);
	}

	fprintf(stderr, "Saving image in %s...\n", file);
	if ((outfile1 = fopen(file, "w")) == NULL) {
		fprintf(stderr, "Cannot open file, %s for writing.\n", file);
		return (0);
	}

	for (i = 0; i < numRows; i++) {
		for (j = 0; j < numCols; j++) {
			int loc = (i * numCols) + j;
			fprintf(outfile1, "%d ", *(r + loc));
		}
		fprintf(outfile1, "\n");
	}

	fclose(outfile1);

	return 0;
}


/*
 * Save the image as PNG ... to be implemented
 */
int RImage::savePNG(char file[100]) {
	/**
	 * TO BE IMPLEMENTED...
	 */
	return 0;
}


/*
 * Determine if a pixel location is within the image bounds
 */
int RImage::inBounds(int i, int j) {
	if ((i < 0) || (j < 0) || (j >= numCols) || (i >= numRows)) {
		return (0);
	} else {
		return (1);
	}
}


/*
 * Convert the image to grayscale
 */
void RImage::convertToGrayscale() {
	if (isColorImage && g != NULL && b != NULL) {
		for (int i = 0; i < (numRows * numCols); i++) {
			r[i] = (r[i] + g[i] + b[i]) / 3;
		}
		delete(g);
		delete(b);
	} else {
		fprintf(stderr, "ERROR: Cannot convert this image to a grayscale representation.\n");
		exit(1);
	}
}


/*
 * Convert the image to color
 */
void RImage::convertToColor() {
	if (!isColorImage) {
		if (g == NULL) {
			initialize(&g, numRows, numCols);
		}
		if (b == NULL) {
			initialize(&b, numRows, numCols);
		}

		for (int i = 0; i < (numRows * numCols); i++) {
			g[i] = r[i];
			b[i] = r[i];
		}

		isColorImage = 1;

	} else {
		fprintf(stderr, "ERROR: Cannot convert this image to a color representation.\n");
		exit(1);
	}
}


/*
 * Make a pretty thermal colored image similar to what is done
 * in MatLab
 */
void RImage::makeThermalColorImage() {
	double mn = 999, mx = -999;

	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			int ij = (i * numCols) + j;

			if (isColorImage == 1) {
				double value = (r[ij] + g[ij] + b[ij]) / 3.0;
				if (value > mx) {
					mx = value;
				}
				if (value < mn) {
					mn = value;
				}
			} else {
				if (r[ij] > mx) {
					mx = r[ij];
				}
				if (r[ij] < mn) {
					mn = r[ij];
				}
			}
		}
	}

	if (isColorImage == 0) {
		initialize(&g, numRows, numCols);
		initialize(&b, numRows, numCols);
		this->isColorImage = 1;
	}


	for (int i = 0; i < numRows; i++) {
		for (int j = 0; j < numCols; j++) {
			int ij = (i * numCols) + j;
			double value;

			g[ij] = b[ij] = r[ij];

			if (isColorImage == 1) {
				double intensity = (r[ij] + g[ij] + b[ij]) / 3.0;
				value = (intensity - mn) / (mx - mn);
			} else {
				value = (r[ij] - mn) / (mx - mn);
			}

			if (value < 0.60) {
				r[ij] = 0;
				g[ij] = 0;
				b[ij] = (int) (exp(-50 * ((value-0.60) * (value-0.60))) * 255);

			} else if (value < 0.70) {
				r[ij] = 0;
				g[ij] = (int) (((value - .60) / .10) * 255);
				b[ij] = 255;

			} else if (value < 0.80) {
				r[ij] = 0;
				g[ij] = 255;
				b[ij] = (int)(((0.80 - value) / .10) * 255);

			} else if (value < 0.90) {
				r[ij] = (int)(((value - .80) / .10) * 255);
				b[ij] = 255;
				b[ij] = 0;

			} else if (value < 1.0) {
				r[ij] = 255;
				g[ij] = (int)(((1.0 - value) / .10) * 255);
				b[ij] = 0;

			} else {
				r[ij] = 255;
				g[ij] = 0;
				b[ij] = 0;

			}
		}
	}
}


/*
 * Allow access to an image channel
 */
int * RImage::getChannel(int channel) {
	switch (channel) {
	case RED:
		return r;
		break;
	case GREEN:
		return g;
		break;
	case BLUE:
		return b;
		break;
	default:
		return r;
	}
}


/*
 * Convert from the JSAMPLE data-structure used by libjpeg to the format
 * used by this class
 */
void RImage::convertJSAMPLEToImage(const JSAMPLE* jpegImageBuffer) {
	if (numCols > 0 && numRows > 0) {
		if (isColorImage == 1) {
			for (int i = 0, k = 0; i < numRows * numCols; i++) {
				*(r + i) = (unsigned int) jpegImageBuffer[k];
				*(g + i) = (unsigned int) jpegImageBuffer[k+1];
				*(b + i) = (unsigned int) jpegImageBuffer[k+2];

				k += 3;

			}
		} else {
			for (int i = 0; i < numRows * numCols; i++) {
				*(r + i) = (unsigned int) jpegImageBuffer[i];
			}
		}

	} else {
		fprintf(stderr, "ERROR: Image dimensions are invalid.\n");
		exit(1);
	}
}


/*
 * Convert to the JSAMPLE data-structure used by libjpeg from the format
 * used by this class
 */
void RImage::convertImageToJSAMPLE(JSAMPLE* jpegImageBuffer) {
	if (numCols > 0 && numRows > 0) {
		if (isColorImage == 1) {
			for (int i = 0, k = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					int ij = (i * numCols) + j;
					jpegImageBuffer[k] = (JSAMPLE) this->r[ij];
					jpegImageBuffer[k+1] = (JSAMPLE) this->g[ij];
					jpegImageBuffer[k+2] = (JSAMPLE) this->b[ij];

					k += 3;
				}
			}

		} else {
			for (int i = 0; i < numRows; i++) {
				for (int j = 0; j < numCols; j++) {
					int ij = (i * numCols) + j;
					jpegImageBuffer[ij] = (JSAMPLE) this->r[ij];
				}
			}
		}

	} else {
		fprintf(stderr, "ERROR: Image dimensions are invalid.\n");
		exit(1);
	}

}


/*
 * Read an integer from a NETBPM formatted file
 */
int RImage::readIntFromFile(FILE *fp) {
	int item, i, flag;

	/* skip forward to start of next number */
	item = getc(fp);
	flag = 1;

	do {

		if (item == '#') { /* comment*/
			while (item != '\n' && item != EOF) {
				item = getc(fp);
			}
		}

		if (item == EOF) {
			return 0;
		}

		if (item >= '0' && item <= '9') {
			flag = 0;
			break;
		}

		/* illegal values */
		if (item != '\t' && item != '\r' && item != '\n' && item != ',') {
			return (-1);
		}

		item = getc(fp);

	} while (flag == 1);

	/* get the number*/
	i = 0;
	flag = 1;

	do {
		i = (i * 10) + (item - '0');
		item = getc(fp);

		if (item < '0' || item > '9') {
			flag = 0;
			break;
		}

		if (item == EOF) {
			break;
		}

	} while (flag == 1);

	return i;

}
