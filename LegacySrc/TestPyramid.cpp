#include <stdlib.h>
#include <stdio.h>
#include "GaussianPyramid.h"
#include "image.h"

using namespace std;

int main() {
	Image mainImg;

	mainImg.read("105.pgm");

	GaussianPyramid pyr(mainImg);

	pyr.processPyramid();


	return 0;
}
