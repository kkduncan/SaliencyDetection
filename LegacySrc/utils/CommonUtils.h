/*
 * CommonUtils.h
 *
 *  Created on: Mar 2, 2011
 *      Author: kesterduncan
 *
 * 	Common utlity functions used for various miscellaneous tasks.
 */

#ifndef COMMONUTILS_H_
#define COMMONUTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>


/*
 * For checking processing time
 */
extern time_t start, stop;

#define STARTTIMER \
	{ \
		start = time(NULL); \
	}

#define ENDTIMER \
	{ \
		stop = time(NULL); \
		time_t total = (stop - start); \
		int h, m, s; \
		h = (long int) (total / 3600); \
		m = (long int) ((total / 60) % 60); \
		s = (long int) (total % 60); \
		infoMsg("TOTAL PROCESSING TIME: %d hour(s), %d minute(s), %d second(s)\n", h, m, s); \
	}


/*
 * Messages for feedback
 */
void infoMsg (const char *format, ...);
void errorMsg (const char *format, ...);
void usageMsg (const char *format, ...);
void writeProgress (char *what, int value, int maximum);

/*
 * File name handling
 */
void stripDirectory(char pathname[255], char output[255]);
char *getFileStem(const char *filename);


#endif /* COMMONUTILS_H_ */
