/*
 * SaliencyMain.cpp
 *
 */
#include <cstdlib>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "SaliencyDetector.h"
#include "EdgeDetector.h"

using namespace std;

int main(int argc, char *argv[]) {
	cv::Mat1f img = cv::imread("bike.jpeg", CV_LOAD_IMAGE_GRAYSCALE);

	if (img.empty()) {
		cout << "No image loaded. Press any key to exit." << endl;
		cin.get();
		return (EXIT_FAILURE);
	}

	sal::ImageSaliencyDetector detector(img);
	detector.setSamplingPercentage(0.10f);

	cout << "Computing..." << endl;
	detector.compute();

	cout << "Post-processing..." << endl;
	detector.performPostProcessing();

	cv::imwrite("SaliencyTestOutput.jpg", detector.getSaliencyMap());

	cout << "Done. Press any key to exit." << endl;

	std::cin.get();

	return 0;
}



