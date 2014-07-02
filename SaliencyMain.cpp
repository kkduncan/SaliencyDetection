/*
 * SaliencyMain.cpp
 *
 */
#include <cstdlib>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "CannyEdgeDetector.h"

using namespace std;

int main(int argc, char *argv[]) {
	sal::CannyEdgeDetector ced;

	cv::Mat1f img = cv::imread("ads.ppm", CV_LOAD_IMAGE_GRAYSCALE);

	ced.compute(img);
	cout << "Writing gradient directions and magnitude images to file\n";

	cv::normalize(ced.getGradOrientations(), img, 0, 255, cv::NORM_MINMAX);
	cv::imwrite("GradDir.jpg", img);
	return 0;
}



