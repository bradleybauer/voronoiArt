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
	std::cout << "Grayscale done" << std::endl;

	// Compute blur
	process.gaussianConvolve(11, data);
	std::cout << "Blur done" << std::endl;

	// Compute sobel filter
	process.sobel(data);
	std::cout << "Sobel done" << std::endl;

	// data = data.array().sqrt();
	data /= data.maxCoeff();

	// Find vertices
	const int max_vertices = 5000;
	Matrix<int, Dynamic, 2> vertices(image.get_height()*image.get_width()/2,2);
	process.computeVertices(vertices, data, max_vertices);
	std::cout << "Vertices done" << std::endl;

	// Color image based on vertices
	process.colorImageWithVertices(image, vertices);

	// Write image to output file
	image.write("output.png");

	return 0;
}
