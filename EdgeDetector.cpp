#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cmath>
#include "EdgeDetector.h"

#ifndef M_PI
#	define M_PI 3.14159265358979323846f
#endif

namespace sal 
{

EdgeDetector::EdgeDetector(const float& sig) {
	if (sig > 0) {
		this->sigma = sig;
	} else {
		this->sigma = 1.20f; // Default
	}

	windowSize = static_cast<int> (1 + (2 * ceil(2.5 * sigma)));

}


EdgeDetector::~EdgeDetector() {

}


void EdgeDetector::compute(const cv::Mat1f& src) {
	assert(!src.empty() && src.cols >= 3 && src.rows >= 3);

	cv::Mat1f smoothedImg;
	cv::Mat1f deltaX, deltaY;

	/*
	 * Perform Gaussian smoothing
	 */
	smoothImage(src, smoothedImg);

	/*
	 * Compute the first derivative in the x and y directions
	 */
	calculateDerivatives(smoothedImg, deltaX, deltaY);

	/*
	 * Compute the direction of the gradient, in radians that
	 * are specified counteclockwise from the positive x-axis
	 */
	calculateGradientDirections(deltaX, deltaY);

	/*
	 * Compute the magnitude of the gradient.
	 */
	calculateGradientMagnitudes(deltaX, deltaY);

}


void EdgeDetector::smoothImage(const cv::Mat1f& src, cv::Mat1f& smoothedImg) {
	assert(!src.empty() && src.cols >= 3 && src.rows >= 3);

	cv::Mat1f tempSmoothedImg(src.rows, src.cols);
	cv::Mat1f gaussianKernel(1, this->windowSize);

	if (!smoothedImg.empty()) {
		std::cout << "EdgeDetector WARNING: Expected an empty matrix for [smoothedImg]\n";
		std::cout << "\t... may have unexpected behavior!\n";
	}

	createGaussianKernel(gaussianKernel);
	smoothedImg.create(src.rows, src.cols);

	/// Blur in the x - direction ///
	int center = this->windowSize / 2;

	for (int r = 0; r < src.rows; r++) {
		for (int c = 0; c < src.cols; c++) {
			float dotProdSum = 0.f;
			float kernelWeightSum = 0.f;

			for (int cc = -center; cc <= center; cc++) {
				if (((c + cc) >= 0) && ((c + cc) < src.cols)) {
					dotProdSum += static_cast<float>(src(r, c + cc) * gaussianKernel(0, center + cc));
					kernelWeightSum += gaussianKernel(0, center + cc);
				}
			}

			if (kernelWeightSum > 0) {
				tempSmoothedImg(r, c) = dotProdSum / kernelWeightSum;
			} else {
				tempSmoothedImg(r, c) = dotProdSum;
			}
		}
	}

	/// Blur in the y - direction ///
	for (int c = 0; c < src.cols; c++) {
		for (int r = 0; r < src.rows; r++) {
			float dotProdSum = 0.f;
			float kernelWeightSum = 0.f;

			for (int rr = -center; rr <= center; rr++) {
				if (((r + rr) >= 0) && ((r + rr) < src.rows)) {
					dotProdSum += tempSmoothedImg(r + rr, c) * gaussianKernel(0, center + rr);
					kernelWeightSum += gaussianKernel(0, center + rr);
				}
			}

			if (kernelWeightSum > 0) {
				smoothedImg(r, c) = dotProdSum / kernelWeightSum + 0.5f;
			} else {
				smoothedImg(r, c) = dotProdSum;
			}
		}
	}
}


void EdgeDetector::createGaussianKernel(cv::Mat1f& kernel) {
	float sum = 0.0f;
	int windowCtr = static_cast<int> (windowSize / 2);

	kernel.create(windowSize, 1);

	for (int i = 0; i < windowSize; i++) {
		float loc = static_cast<float>(i - windowCtr);
		float val = powf(2.71828f, -0.5f * loc * loc / (sigma * sigma)) / (sigma * sqrt(6.2831853f));
		kernel(0, i) = val;
		sum += val;
	}

	assert(sum > 0);
	for (int i = 0; i < windowSize; i++) {
		kernel(0, i) /= sum;
	}

}


void EdgeDetector::calculateDerivatives(const cv::Mat1f& smoothedImg, cv::Mat1f& deltaX, cv::Mat1f& deltaY) {
	deltaX.create(smoothedImg.rows, smoothedImg.cols);
	deltaY.create(smoothedImg.rows, smoothedImg.cols);

	/*
	 * Adjusted to handle pixels at borders
	 */
	for (int r = 0; r < smoothedImg.rows; r++) {
		for (int c = 0; c < smoothedImg.cols; c++) {
			if (c == 0) {
				deltaX(r, c) = smoothedImg(r, c + 1) - smoothedImg(r, c);

			} else if (c == (smoothedImg.cols - 1)) {
				deltaX(r, c) = smoothedImg(r, c) - smoothedImg(r, c - 1);

			} else {
				deltaX(r, c) = smoothedImg(r, c + 1) - smoothedImg(r, c - 1);
			}
		}
	}

	for (int r = 0; r < smoothedImg.rows; r++) {
		for (int c = 0; c < smoothedImg.cols; c++) {
			if (r == 0) {
				deltaY(r, c) = smoothedImg(r + 1, c) - smoothedImg(r, c);

			} else if (r == (smoothedImg.rows - 1)) {
				deltaY(r, c) = smoothedImg(r, c) - smoothedImg(r - 1, c);

			} else {
				deltaY(r, c) = smoothedImg(r + 1, c) - smoothedImg(r - 1, c);
			}
		}
	}
}


void EdgeDetector::calculateGradientDirections(const cv::Mat1f& deltaX, const cv::Mat1f& deltaY) {
	assert(deltaX.size == deltaY.size);
	gradOrientations.create(deltaX.rows, deltaX.cols);

	for (int r = 0; r < deltaX.rows; r++) {
		for (int c = 0; c < deltaY.cols; c++) {
			float dx = deltaX(r, c);
			float dy = -deltaY(r, c);
			gradOrientations(r, c) = calculateVectorAngle(dx, dy);
		}
	}
}


void EdgeDetector::calculateGradientMagnitudes(const cv::Mat1f& deltaX, const cv::Mat1f& deltaY) {
	assert(deltaX.size == deltaY.size);

	gradMagnitudes.create(deltaX.rows, deltaX.cols);

	for (int r = 0; r < deltaX.rows; r++) {
		for (int c = 0; c < deltaX.cols; c++) {
			float sqX = deltaX(r, c) * deltaX(r, c);
			float sqY = deltaY(r, c) * deltaY(r, c);
			gradMagnitudes(r, c) = sqrt(sqX + sqY);
		}
	}
}


float EdgeDetector::calculateVectorAngle(const float& x, const float& y) {
	float xu = 0.f, yu = 0.f, angle = 0.f;
	yu = sin(x - y);
	xu = cos(x - y);

	if ((xu == 0) && (yu == 0)) {
		return (0.f);
	}

	angle = atan2(yu, xu);

	if (x >= 0) {
		if (y >= 0) {
			return (angle);
		} else {
			return (2.f * M_PI - angle);
		}

	} else {
		if (y >= 0) {
			return (M_PI - angle);
		} else {
			return (M_PI + angle);
		}
	}
}







} /* namespace sal */
