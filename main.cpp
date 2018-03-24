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
    std::cout << "Converting to grayscale" << std::endl;
    process.intensity(image, data);

    // Compute blur
    std::cout << "Applying blur" << std::endl;
    process.gaussianConvolve(11, data);
    process.gaussianConvolve(11, data);
    process.gaussianConvolve(11, data);

    // Compute sobel filter
    std::cout << "Applying sobel" << std::endl;
    process.sobel(data);

    // data = data.array().sqrt();
    data /= data.maxCoeff() * .5;

    // Find vertices
    const int max_vertices = 5000;
    // each row in vertices contains (P, C, c)
    // P = 2D position of vertex
    // C = 3D color data
    // c = 1D number of pixels belonging to vertex
    Matrix<long, Dynamic, 6> vertices(image.get_height()*image.get_width() / 2, 6);
    vertices.setZero();
    std::cout << "Finding vertices" << std::endl;
    process.computeVertices(vertices, data, max_vertices);

    // Color image based on vertices
    process.colorImageWithVertices(image, vertices);

    // Write image to output file
    image.write("output.png");

    return 0;
}
