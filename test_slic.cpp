#include "slic.h"
#include <cstdlib>
#include <cstdio>

using namespace std;

int main(int argc, char** argv)
{
	if (argc != 3) {
		printf("usage: test_slic <filename> <number of superpixels>\n");
		exit(-1);
	}

	cv::Mat img, result;
	
	img = imread(argv[1]);
	int numSuperpixel = atoi(argv[2]);

	SLIC slic;
	slic.GenerateSuperpixels(img, numSuperpixel);
	if (img.channels() == 3) 
		result = slic.GetImgWithContours(cv::Scalar(0, 0, 255));
	else
		result = slic.GetImgWithContours(cv::Scalar(128));

	cv::imwrite("result.jpg", result);
}