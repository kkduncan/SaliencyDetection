# Kester Duncan
# University of South Florida - Computer Science and Engineering
#
# Makefile for Saliency computations.
#
# (c) Copyright 2011 - last Updated 05/15/11
#
#
#  DISCLAIMER:
#  ==============================================================
#  This code was not designed well. It is not modular and doesn't
#  obey the Object Oriented Principles. Therefore it may contain
#  many bugs. Use at your own risk!
#  ==============================================================
#
#

# !!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Ensure that jpeg library is built first before entering the
# image directory and running ./setupJPEG.sh
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
LIBS 	= image/lib/libjpeg.a
EXEC1	= bin/ComputeVideoSaliency
VID 	= ImageSequence.cpp
CC 		= g++
C 		= gcc
FLAGS 	= -g -Wall


CommonUtils.o: CommonUtils.c
	$(CC) -c CommonUtils.c
	
canny.o: canny.c canny.h
	$(CC) -c canny.c canny.h

RImage.o: image/RImage.cpp image/RImage.h
	$(CC) -c image/RImage.cpp image/RImage.h
	
gaussian_filter.o: gaussian_filter.c gaussian_filter.h
	$(CC) -c gaussian_filter.c gaussian_filter.h
						
IntegralImage.o: image/IntegralImage.cpp image/IntegralImage.h RImage.o
	$(CC) -c image/IntegralImage.cpp image/IntegralImage.h
	
OBJ = gaussian_filter.o canny.o RImage.o CommonUtils.o
		

video: ComputeVideoSaliency.cpp IntegralImage.o $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC1) $(FLAGS) ComputeVideoSaliency.cpp gaussian_filter.o RImage.o CommonUtils.o IntegralImage.o $(VID) $(LIBS)
	@echo ' '
	
all: video

clean:
	rm *.o; rm *.h.gch
