# Kester Duncan
# University of South Florida - Computer Science and Engineering
#
# Makefile for Saliency computations.
#
# (c) Copyright 2011 - last Updated 05/15/11

EXEC1= ../bin/ComputeSaliency
EXEC2= ../bin/ComputeIterativeSaliency
EXEC3= ../bin/ComputeVideoSaliency
EXEC4= ../bin/ComputePyramid
EXEC5= ../bin/CalculatePerformance
EXEC6= ../bin/GetSaliency
EXEC7= ../bin/GetPyramidSaliency

LIBS = image/lib/libjpeg.a
#VID = video/XYT.cpp video/Sequence.cpp video/Array.cpp video/StdAfx.cpp
VID = sequences/ImageSequence.cpp
CC = g++
C = gcc
FLAGS = -g -Wall

CommonUtils.o: utils/CommonUtils.c
	$(CC) -c utils/CommonUtils.c
	
canny.o: edge/canny.c edge/canny.h
	$(CC) -c edge/canny.c edge/canny.h

RImage.o: image/RImage.cpp image/RImage.h
	$(CC) -c image/RImage.cpp image/RImage.h
	
histogram.o: histogram/histogram.c
	$(CC) -c histogram/histogram.c histogram/histogram.h
	
gaussian_filter.o: gaussian/gaussian_filter.c gaussian/gaussian_filter.h
	$(CC) -c gaussian/gaussian_filter.c gaussian/gaussian_filter.h
	
LBP.o: image/LBP.cpp image/LBP.h
	$(CC) -c image/LBP.cpp image/LBP.h RImage.o
	
GaussianPyramid.o: gaussian/GaussianPyramid.cpp gaussian/GaussianPyramid.h
	$(CC) -c gaussian/GaussianPyramid.cpp gaussian/GaussianPyramid.h
	
IntegralImage.o: image/IntegralImage.cpp image/IntegralImage.h RImage.o
	$(CC) -c image/IntegralImage.cpp image/IntegralImage.h
	
SauvolasThreshold.o: image/SauvolasThreshold.cpp image/SauvolasThreshold.h
	$(CC) -c image/SauvolasThreshold.cpp image/SauvolasThreshold.h

REM.o: REM.cpp REM.h RImage.o histogram.o
	$(CC) -c REM.cpp REM.h
		
OBJ = gaussian_filter.o canny.o RImage.o histogram.o REM.o CommonUtils.o
		
saliency: ComputeSaliency.cpp $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC1) ComputeSaliency.cpp $(OBJ) $(LIBS)
	@echo ' '
	
iterative: ComputeIterativeSaliency.cpp $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC2) ComputeIterativeSaliency.cpp $(OBJ) $(LIBS)
	@echo ' '
	
video: ComputeVideoSaliency.cpp IntegralImage.o SauvolasThreshold.o $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC3) $(FLAGS) ComputeVideoSaliency.cpp gaussian_filter.o RImage.o histogram.o REM.o CommonUtils.o image/SauvolasThreshold.cpp image/IntegralImage.cpp $(VID) $(LIBS)
	@echo ' '
	
pyramid: ComputePyramidSaliency.cpp IntegralImage.o SauvolasThreshold.o GaussianPyramid.o $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC4) ComputePyramidSaliency.cpp gaussian/GaussianPyramid.cpp image/SauvolasThreshold.cpp image/IntegralImage.cpp $(OBJ) $(LIBS)
	@echo ' '
	
general: GetSaliency.cpp LBP.o $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC6) GetSaliency.cpp LBP.o $(OBJ) $(LIBS)
	
pgeneral: GetPyramidSaliency.cpp LBP.o GaussianPyramid.o $(OBJ) $(LIBS)
	@echo 'Building: $@'
	$(CC) -o $(EXEC7) GetPyramidSaliency.cpp gaussian/GaussianPyramid.cpp LBP.o $(OBJ) $(LIBS)
	
performance: CalculatePerformance.cpp RImage.o
	@echo 'Building: $@'
	$(CC) -o $(EXEC5) CalculatePerformance.cpp RImage.o $(LIBS)
		
all: saliency iterative video pyramid

clean:
	rm *.o; rm *.h.gch
