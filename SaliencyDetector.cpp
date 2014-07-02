#include <cstdlib>
#include "CannyEdgeDetector.h"
#include "SaliencyDetector.h"


namespace sal {

SaliencyDetector::SaliencyDetector() : samplingPercentage(0.25f), neighborhoodSize(7) {

}


SaliencyDetector::~SaliencyDetector() {

}


void SaliencyDetector::postProcessSaliencyMap(cv::Mat1f& salMap, const float& sigma){
	CannyEdgeDetector canny(sigma);
	cv::Mat1f filteredMap;

	/*
	 * First smoothing operation
	 */
	canny.smoothImage(salMap, filteredMap);
	filteredMap.copyTo(salMap);

	/*
	 * Smoothing another time effectively gets rid of
	 * unsightly rings caused by the Gaussian
	 */
	filteredMap.release();
	canny.setSigma(0.6f * sigma);
	canny.smoothImage(salMap, filteredMap);
	filteredMap.copyTo(salMap);

	/*
	 *  Normalization
	 */
	double maxVal = 0, minVal = 0;
	cv::minMaxLoc(salMap, &minVal, &maxVal);

	for (int i = 0; i < salMap.rows; i++) {
		for (int j = 0; j < salMap.cols; j++) {
			salMap(i, j) = ((((salMap(i, j) - static_cast<float>(minVal))
					/ static_cast<float>(maxVal - minVal)) * 255.0));
		}
	}
}


} /* namespace sal */
