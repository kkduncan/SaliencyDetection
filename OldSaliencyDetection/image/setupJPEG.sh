#!/bin/bash
# This bash script simply sets up the IJPEG library to be used
# by the RImage class for reading and writing JPEG images
#
# K. Duncan
# University of South Florida

# !!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!
# The following variable must reflect the current absolute
# location of the 'image' directory on the machine this code 
# is running on.
IMAGE_DIR="/home/carrt/workspace/SaliencyDetection/OldSaliencyDetection/image";
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cd jpeg-8c/
./configure --prefix="$IMAGE_DIR" --libdir="$IMAGE_DIR/lib"

make
make install
