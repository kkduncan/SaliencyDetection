/*
 * CommonUtils.c
 *
 *  Created on: Mar 2, 2011
 *      Author: kesterduncan
 */

#include "CommonUtils.h"


void infoMsg(const char *format, ...) {
	va_list args;

	fprintf(stdout, "[INFO] - " );

	va_start(args, format);
	vfprintf(stdout, format, args );

	va_end( args );

	fprintf(stdout, "\n" );
	fflush(stdout);

}

void errorMsg(const char *format, ...) {
	va_list args;

	fprintf(stderr, "[ERROR] - " );

	va_start(args, format);
	vfprintf(stderr, format, args );

	va_end( args );
	fprintf(stderr, "\n" );
	fflush(stderr);
	exit(1);

}

void usageMsg (const char *format, ...) {
	va_list args;

	fprintf(stderr, "\n\tUSAGE: ");

	va_start(args, format);
	vfprintf(stderr, format, args );

	va_end( args );
	fprintf(stderr, "\n" );
	fflush(stderr);
	exit(1);

}

void writeProgress (char *what, int value, int maximum) {
	char *symbols = "\\|/-";
	static char *sym = "";

	if (*sym == '\0') {
		sym = symbols;
	}

	char msg[255];

	if (maximum) {
		sprintf(msg, "%s: %c processed: %11d (%3d%%)    ", what, *sym++, value, 100 * value / maximum);
		fprintf(stderr, "%s\r", msg);
		strcpy(msg, "");
	} else {
		sprintf(msg, "%s: %c processed: %11d           ", what, *sym++, value);
		fprintf(stderr, "%s\r", msg);
		strcpy(msg, "");
	}

	fflush(stdout);
}


void stripDirectory(char pathname[255], char output[255]) {
	int i = 0, j = 0;

	while (pathname[i] != '\0') {
		if (pathname[i] == '/')
			j = 0;
		else
			output[j++] = pathname[i];
		i++;
	}
	output[j++] = '\0';
}


char *getFileStem(char *filename) {
	char *str = (char*) malloc((strlen(filename) - 3) * sizeof(char));
	char *chPtr = strrchr(filename, '.');
	int index = chPtr - filename;
	int i = 0;;

	while (i < index) {
		str[i] = filename[i];
		i++;
	}
	str[i] = '\0';

	return str;
}
