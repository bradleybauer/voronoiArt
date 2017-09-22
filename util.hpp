#include <eigen3/Eigen/Eigen>
#include <png++/png.hpp>
using namespace Eigen;
class Process {
public:
	Process(){};

	void gaussianConvolve(int kernel_size, MatrixXd& data) {
		// compute gaussian kernel
		// compute zero padded convolution of data with kernel
		//    store convolution in temporary matrix
		// copy temporary matrix to data
	}

	void gradient(MatrixXd& data) {
		// compute zero padded convolution of data with [-1, 0, 1]
		//    store convolution in temporary matrix
		// copy temporary matrix to data
	}

	void intensity(png::image<png::rgb_pixel>& image, MatrixXd& data) {
		// set data_ij = length(image_ij)^2

		// for (int i = 0; i < data_mat.cols(); ++i)
		// 	for (int j = 0; j < data_mat.rows(); ++j)
		// 		data_mat(j,i) = image[j][i];
	}

	void computeVertices(Matrix<double, Dynamic, 2> vertices,
	                     MatrixXd& data) {
    // if distribution(0,1) < data_ij
    //     add vertex(i, j) to vertices
	}

	void colorImageWithVertices(png::image<png::rgb_pixel>& image,
	                            Matrix<double, Dynamic, 2> vertices) {
		// loop through image pixels
		//    find vertex that is closest to pixel
		//    set pixel to average color of image around point vertex
	}

};

