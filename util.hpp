#include <eigen3/Eigen/Eigen>
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

};

