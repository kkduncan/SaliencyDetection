/*
 * Software License Agreement (BSD License)
 *
 *  Saliency Detection
 *  Copyright (c) 2011, Kester Duncan
 *
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *	\file	SaliencyDetector.h
 *	\brief	Defines the image and video saliency detection methods
 *	\author	Kester Duncan
 */

#ifndef SALIENCYDETECTOR_H_
#define SALIENCYDETECTOR_H_

#include <opencv2/core/core.hpp>



/**
 * Saliency namespace
 */
namespace sal {

/**
 * \brief Abstract class for saliency mechanisms following the theory
 * found in the paper "Relational Entropy-Based Saliency Detection in Images and Video"
 * by Kester Duncan and Sudeep Sarkar in IEEE International Conference on Image Processing (2012)
 *
 * For details on the algorithm, the interested reader is advised to read
 * the article, found here:
 *
 * http://www.cse.usf.edu/~kkduncan/research/DuncanICIP2012.pdf
 *
 */
class SaliencyDetector {
public:
	SaliencyDetector();
	virtual ~SaliencyDetector();

	/// Compute the saliency of an image or a 3D volume of frames
	void compute() = 0;

	/// Quantize the magnitude values in an effort to suppress lows and boost highs
	void quantizeMagnitudes() = 0;

	/// Process resultant saliency map to make it smooth and pretty
	void postProcessSaliencyMap(cv::Mat1f& salMap, const float& sigma = 1.20f);


protected:
	/// Percentage of image / frame pixels to use for calculations
	float samplingPercentage;

	/// Size of one dimension of the pixel neighborhood
	int neighborhoodSize;
};



/**
 * \brief Image saliency detector
 */
class ImageSaliencyDetector : public SaliencyDetector {


};



/**
 * \brief Video saliency detector
 */
class VideoSaliencyDetector : public SaliencyDetector {

};




} /* namespace sal */









#endif /* SALIENCYDETECTOR_H_ */
