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

int main(int argc, char *argv[]) 
{
	try
	{
		if (argc == 2)
		{
			std::string imgName(argv[1]);
			if (!imgName.empty())
			{
				cv::Mat1f img = cv::imread(imgName, CV_LOAD_IMAGE_GRAYSCALE);
				std::string outputName(imgName);
				outputName.append("_SAL.png");

				if (img.empty()) 
				{
					cout << "Error: No image loaded. Press any key to exit." << endl;
					cin.get();
					return (EXIT_FAILURE);
				}

				sal::ImageSaliencyDetector detector(img);
				detector.setSamplingPercentage(0.25f);
				detector.setNeighborhoodSize(5);
				
				cout << "Computing..." << endl;
				detector.compute();

				cout << "Post-processing..." << endl;
				//detector.performPostProcessing();

				cv::imwrite(outputName, detector.getSaliencyMap());
				cout << "Output written to --> [ " << outputName << " ]" << std::endl;
			}
		}
		else
		{
			cout << "Error: No image provided.\n";
			cout << "Usage: ./SaliencyDetector.exe [IMAGE_PATH]\n";
		}
	}
	catch (const std::exception &ex)
	{
		cout << "Error: An unexpected error occurred during detection: " << ex.what() << endl;		
	}
	catch (...)
	{
		cout << "Error: An unknown error occurred during detection.\n";
	}

	//cout << "Done. Press any key to exit." << endl;
	//std::cin.get();

	return 0;
}



