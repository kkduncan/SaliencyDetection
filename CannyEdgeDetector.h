/*
 * CannyEdgeDetector.h
 *
 *  Created on: Jun 27, 2014
 *      Author: carrt
 */

#ifndef CANNYEDGEDETECTOR_H_
#define CANNEDGEDETECTOR_H_

#include <opencv2/core/core.hpp>


/// Saliency Detection namespace
namespace sal {


/*
 * \brief Class to perform Canny edge detection
 *
 * Acquires the gradient direction and gradient magnitudes values of an
 * input image provided via the compute function
 */
class CannyEdgeDetector {
public:
	explicit CannyEdgeDetector(const float& sig = 1.20f);
	virtual ~CannyEdgeDetector();

	/// Perform Canny edge detection on the given image
	void compute(const cv::Mat1f& src);

	/// Smooth the input image and store the result in smoothedImg
	void smoothImage(const cv::Mat1f& src, cv::Mat1f& smoothedImg);


	/// Get sigma value
	float getSigma() const {
		return sigma;
	}

	/// Get the Gaussian window size
	int getWindowSize() const {
		return windowSize;
	}

	/// Get the gradient magnitudes
	cv::Mat1f getGradMagnitudes() const {
		return gradMagnitudes;
	}

	/// Get the gradientOrientations
	cv::Mat1f getGradOrientations() const {
		return gradOrientations;
	}

	/// Set the sigma value (should be greater than 0)
	void setSigma(const float& value) {
		if (value > 0) {
			sigma = value;
		}
	}


private:
	/// Create a 1-D Gaussian kernel
	void createGaussianKernel(cv::Mat1f& kernel);

	/// Calculate derivatives in the x and y directions
	void calculateDerivatives(const cv::Mat1f& smoothedImg, cv::Mat1f& deltaX, cv::Mat1f& deltaY);

	/// Calculate the gradient directions
	void calculateGradientDirections(const cv::Mat1f& deltaX, const cv::Mat1f& deltaY);

	/// Calculate the gradient magnitudes
	void calculateGradientMagnitudes(const cv::Mat1f& deltaX, const cv::Mat1f& deltaY);

	/// Calculate the angle of the vector represented by x and y
	float calculateVectorAngle(const float& x, const float& y);


	/// Standard deviation for the Gaussian smoothing filter
	float sigma;

	/// The Gaussian window size based on sigma
	int windowSize;

	/// Gradient magnitudes
	cv::Mat1f gradMagnitudes;

	/// Gradient orientations / directions
	cv::Mat1f gradOrientations;
};

} /* namespace sal */

#endif /* CANNYEDGEDETECTOR_H_ */
