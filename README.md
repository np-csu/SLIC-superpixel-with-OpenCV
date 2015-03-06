# SLIC Superpixel Modification

We modified the [SLICO version of slic algorithm](http://ivrg.epfl.ch/research/superpixels) to work with OpenCV2, and make it easy to use.

## Class explanation

* `void SLIC::GenerateSuperpixels(cv::Mat& img, UINT numSuperpixels)`

	Perform SLIC algorithm on the given image with the given number of superpixels.


* `cv::Mat SLIC::GetImgWithContours(cv::Scalar color)`

	Get the result image with contours on the given color.

* `int* SLIC::GetLabel()`

	Get label on each pixel which shows the number of superpixel it belongs to.


## Example

original image

![bird_color.jpg](bird_color.jpg)

200 superpixels

![bird_color.jpg](bird_color_200.jpg)

500 superpixels

![bird_color.jpg](bird_color_500.jpg)

