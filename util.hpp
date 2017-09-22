#include <eigen3/Eigen/Eigen>
#include <png++/png.hpp>
class Process {
public:
	Process(){};

	void gaussianConvolve(int kernel_size, Eigen::MatrixXd& data) {
		// compute gaussian kernel
		// compute zero padded convolution of data with kernel
		//    store convolution in temporary matrix
		// copy temporary matrix to data
	}

	void gradient(Eigen::MatrixXd& data) {
		// compute zero padded convolution of data with [-1, 0, 1]
		//    store convolution in temporary matrix
		// copy temporary matrix to data
	}

	void intensity(png::image<png::rgb_pixel>& image, Eigen::MatrixXd& data) {
		// set data_ij = length(image_ij)^2
	}

	void computeVertices(Eigen::Matrix<double, Dynamic, 2> vertices,
	                     Eigen::MatrixXd& data) {
    // if distribution(0,1) < data_ij
    //     add vertex(i, j) to vertices
	}

	void colorImageWithVertices(png::image<png::rgb_pixel>& image,
	                            Eigen::Matrix<double, Dynamic, 2> vertices) {
		// loop through image pixels
		//    find vertex that is closest to pixel
		//    set pixel to average color of image around point vertex
	}

};

