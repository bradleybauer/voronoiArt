#include <iostream>
#include <png++/png.hpp>
#include "util.hpp"

int main() {
	// Possibly accept arguments from user using gflags?

	// Load image using png++
	png::image< png::rgb_pixel_16 > image("input.png");

	// Compute transformations
	Process process;

	// Get intensities
	MatrixXd data(image.get_height(), image.get_width());
	process.intensity(image, data);
	std::cout << data.maxCoeff() << std::endl;

	// Compute Blur
	process.gaussianConvolve(11, data);

	// Compute gradient
	process.gradient(data);
	std::cout << data.maxCoeff() << std::endl;

	data /= data.maxCoeff();

	// Find vertices
	Matrix<int, Dynamic, 2> vertices(image.get_height()*image.get_width()/2,2);
	process.computeVertices(vertices, data);

	// Color image based on vertices
	process.colorImageWithVertices(image, vertices);

	// Write image to output file
	image.write("output.png");

	return 0;
}
