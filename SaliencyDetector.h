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
	/*!
	 * Information structure to store the information required
	 * for performing iterative kernel density estimation for local
	 * pixel neighborhoods.
	 *
	 * Therefore, for each pixel, this information is calculated and
	 * stored. See Section 3 of the aforementioned paper for details.
	 */
	struct KernelDensityInfo {
		/// Constructor
		KernelDensityInfo() :
			kernelSum(0.f),
			entropy(0.f),
			firstWeight(0.f),
			secondWeight(0.f),
			sampleCount(0) {

		}

		/// Copy constructor
		KernelDensityInfo(const KernelDensityInfo& other) {
			this->kernelSum = other.kernelSum;
			this->entropy = other.entropy;
			this->firstWeight = other.firstWeight;
			this->secondWeight = other.secondWeight;
			this->sampleCount = other.sampleCount;
		}

		/// Assignment operator
		KernelDensityInfo& operator=(const KernelDensityInfo& other) {
			this->kernelSum = other.kernelSum;
			this->entropy = other.entropy;
			this->firstWeight = other.firstWeight;
			this->secondWeight = other.secondWeight;
			this->sampleCount = other.sampleCount;

			return *this;
		}

		/// The calculated kernel density sum
		float kernelSum;

		/// The entropy info for a pixel in relation to its neighbors
		float entropy;

		/// The value of the first weight used in the kernel density calculations
		float firstWeight;

		/// The value of the first weight used in the kernel density calculations
		float secondWeight;

		/// The number of samples used to estimate the kernel density
		int sampleCount;
	};

	/// Alias for 2D Pixel location
	typedef cv::Point2i Location2D;

	/// Alias for 3D Pixel location
	typedef cv::Point3i Location3D;

	/*!
	 * Stores the applicable bounding region of the 4 pixels forming
	 * the two samples used for incremental kernel density estimation
	 */
	struct BoundingBox2D {
		Location2D topLeft;
		Location2D topRight;
		Location2D botLeft;
		Location2D botRight;
	};


public:
	SaliencyDetector() : samplingPercentage(0.25f), neighborhoodSize(5) { };
	virtual ~SaliencyDetector() { };

	/// Compute the saliency of an image or a 3D volume of frames
	virtual void compute() = 0;

	void setSamplingPercentage(const float& p) {
		if (p > 0) samplingPercentage = p;
	}

	void setNeighborhoodSize(const int& m) {
		if (m > 0) neighborhoodSize = m;
	}

	float getSamplingPercentage() const {
		return samplingPercentage;
	}

	int getNeighborhoodSize() const {
		return neighborhoodSize;
	}


protected:
	/// Used to handle special cases which can produce errors (such as numbers below 1e-15)
	enum {ERROR_FLAG = -1};

	/// Quantize the magnitude values in an effort to suppress lows and boost highs
	virtual void quantizeMagnitudes() = 0;

	/// Process resultant saliency map to make it smooth and pretty
	void postProcessSaliencyMap(cv::Mat1f& salMap, const float& sigma = 18.0f);

	/// Percentage of image / frame pixels to use for calculations
	float samplingPercentage;

	/// Size of one dimension of the pixel neighborhood
	int neighborhoodSize;
};


/// Convenience aliases
typedef SaliencyDetector::KernelDensityInfo KernelDensityInfo;
typedef SaliencyDetector::BoundingBox2D BoundingBox2D;




/**
 * \brief 2D Image saliency detector
 */
class ImageSaliencyDetector : public SaliencyDetector {
public:
	/// Constructor that takes a source image as input
	ImageSaliencyDetector (const cv::Mat1f& src);

	/// Destructor
	~ImageSaliencyDetector();

	/*!
	 * Computes the saliency of pixels based on the entropy of distance
	 * and gradient orientation relationships
	 */
	void compute();

	/*!
	 * Quantize the 2D magnitude values by suppressing lows and boosting highs
	 */
	void quantizeMagnitudes();

	/*!
	 * Perform class specific post-processing of the saliency map
	 */
	void performPostProcessing();


public:
	void setMagnitudes(const cv::Mat1f& other) {
		magnitudes = other;
	}

	void setOrientations(const cv::Mat1f& other) {
		orientations = other;
	}

	void setSourceImage(const cv::Mat1f& theSrc) {
		srcImage = theSrc;
	}

	cv::Mat1f getMagnitudes() const {
		return magnitudes;
	}

	cv::Mat1f getOrientations() const {
		return orientations;
	}

	cv::Mat1f getSourceImage() const {
		return srcImage;
	}

	cv::Mat1f getSaliencyMap() const {
		return saliencyMap;
	}


private:
	/*!
	 * Calculate the intermediate kernel sum from the contribution of the
	 * pixels given by samples
	 */
	KernelDensityInfo calculateKernelSum(const std::vector<Location2D>& samples);

	/*!
	 * Gets the applicable region of the pixels forming two samples
	 * that can be updated
	 */
	BoundingBox2D getApplicableBounds(const std::vector<Location2D>& samples);

	/*!
	 * Updating the kernel sums of the pixels within the applicable region given
	 */
	void updateApplicableRegion(const BoundingBox2D& bounds, const KernelDensityInfo& kernelSum);

	/*!
	 * Updates the saliency map using the most recent density estimates
	 */
	void updateSaliencyMap();

	/*!
	 * Determine if sample locations are within the valid image boundaries
	 */
	bool inValidImageBounds(const std::vector<Location2D>& samples);

	/*!
	 * Updates the entropy of a pixel given the iteratively-estimated distribution of
	 * the distance and orientation relationships around it
	 */
	void updatePixelEntropy(KernelDensityInfo& kernelInfo);


private:
	cv::Mat1f magnitudes;
	cv::Mat1f orientations;
	cv::Mat1f srcImage;
	cv::Mat1f saliencyMap;
	std::vector< std::vector<KernelDensityInfo> > densityEstimates;

};


/**
 * \brief Video saliency detector
 */
class VideoSaliencyDetector : public SaliencyDetector {
	// TODO: re-implement
};




} /* namespace sal */









#endif /* SALIENCYDETECTOR_H_ */
