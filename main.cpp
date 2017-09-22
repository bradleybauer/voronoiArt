#include <iostream>
#include <png++/png.hpp>
#include "util.hpp"

int main() {
	// Possibly accept arguments from user using gflags?

	// Load image using png++
	png::image< png::rgb_pixel > image("input.png");

	// Compute transformations
	Process process;

	// Get intensities
	MatrixXd data(image.get_height(), image.get_width());
	process.intensity(image, data);

	// Compute gradient
	process.gradient(data);

	// Blur intensity gradient
	process.gaussianConvolve(25, data);

	// Find vertices
	Matrix<double, Dynamic, 2> vertices;
	process.computeVertices(vertices, data);

	// Color image based on vertices
	process.colorImageWithVertices(image, vertices);

	// Write image to output file
	image.write("output.png");

	return 0;
}
