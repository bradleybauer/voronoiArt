#include <eigen3/Eigen/Eigen>
#include <png++/png.hpp>
#include <cmath>
#include <random>
#include <ctime>
using namespace Eigen;
class Process {
public:
	Process(){};

	void gaussianConvolve(int kernel_size, MatrixXd& data) {
		if (kernel_size % 2 == 0) {
			std::cout << "Please use odd sized kernels" << std::endl;
			exit(1);
		}
		// compute gaussian kernel
		MatrixXd kernel(kernel_size,kernel_size);
		for (int i = 0; i < kernel_size; ++i) {
			for (int j = 0; j < kernel_size; ++j) {
				// Generate coordinates
				Vector2d p = {i, j};
				p /= double(kernel_size-1);
				p = p.array() * 2. - 1.;
				p = p.array() * p.array();
				p*=6.;
				const double sigmaSqrd = 20;
				kernel(i,j) = 1./(2.*3.141592*sigmaSqrd)*std::exp(-p.sum()/(2.*sigmaSqrd));
			}
		}
		convolve_mirror_pad(data, kernel);
	}

	void gradient(MatrixXd& data) {
		MatrixXd D(3,3);
		D << -1,0,1,
		     -2,0,2,
		     -1,0,1;
		MatrixXd temp = data;
		convolve_mirror_pad(temp, D);
		D = D.transpose().eval();
		convolve_mirror_pad(data, D);
		data = Eigen::sqrt(data.array().pow(2) + temp.array().pow(2)).eval();
	}

	void intensity(png::image<png::rgb_pixel_16>& image, MatrixXd& data) {
		for (int i = 0; i < data.rows(); ++i) {
			for (int j = 0; j < data.cols(); ++j) {
				double r = image[i][j].red/65535.;
				double g = image[i][j].green/65535.;
				double b = image[i][j].blue/65535.;
				data(i,j) = sqrt((r*r + g*g + b*b)/3.);
			}
		}
	}

	// Thoughts on making the image coloring faster
	//     Keep a list of list of distance between vertexes
	//     Each vertex has a row in this matrix/listOfLists
	//     Each element in a row represents the distance of that vertex to a number of other vertexes
	//     To compute the color of the image, first choose a row in this list of lists
	//     Making that choice is equivelent to choosing a vertex
	//     Lets call that vertex the vertex X
	//     Now do a recursive search for pixels whose closest distanace out of all vertices in X's list is X
	//     This reduces the amount of vertexes we have to compare distances to by a lot
	//     I imagine that each row/sublist in this list of lists will have about 6 members
	//     So for each pixel in the image we have to loop about 6 times, we have to compute 6 distances
	//     So the number of distance computations, the num of loops, is 6 * width * height
	//     instead of being Num_Vertices * width * height.
	//     However, to compute the list of lists requires Num_Vertices^2 loops?
	//     That might make this optimization less efficient
	void computeVertices(Matrix<int, Dynamic, 2>& vertices,
	                     MatrixXd& data) {
		auto eng = std::default_random_engine(std::time(0));
		auto dist = std::uniform_real_distribution<double>(0, 1);
		int count = 0;
		for (int i = 0; i < data.rows(); ++i) {
			for (int j = 0; j < data.cols(); ++j) {
				if (dist(eng) < data(i,j)) {
					vertices(count, 0) = i;
					vertices(count, 1) = j;
					count++;
				}
			}
		}
		vertices.conservativeResize(count, 2);
	}

	void colorImageWithVertices(png::image<png::rgb_pixel_16>& image,
	                            Matrix<int, Dynamic, 2>& vertices) {
		// loop through image pixels
		//    find vertex that is closest to pixel
		//    set pixel to average color of image around point vertex
		int height = image.get_height();
		int width = image.get_width();
		for (int i = 0; i < height; ++i) {
			std::cout << "Progress:" << i/float(image.get_height()) <<std::endl;
			for (int j = 0; j < width; ++j) {
				int indx_of_min = 0;
				double d = std::numeric_limits<double>::infinity();
				for (int k = 0; k < vertices.rows(); ++k) {
					double x = j - vertices(k, 1);
					double y = i - vertices(k, 0);
					double dd = std::sqrt((x*x + y*y));
					// double dd = std::fabs(x) + std::fabs(y);
					if (dd < d) {
						indx_of_min = k;
						d = dd;
					}
				}
				png::rgb_pixel_16 pix = image[vertices(indx_of_min,0)][vertices(indx_of_min,1)];
				double r = pix.red;
				double g = pix.green;
				double b = pix.blue;
				image[i][j] = png::rgb_pixel_16(r,g,b);
			}
		}
	}

private:
	void convolve_mirror_pad(MatrixXd& data, MatrixXd& kernel) {
		if (kernel.rows() < 3 || kernel.cols() < 3) {
			std::cout << "Please use a bigger kernel" << std::endl;
			exit(1);
		}
		// compute mirror padded convolution of data with kernel
		//    store convolution in temporary matrix
		MatrixXd temp(data.rows(), data.cols());
		int kern_spread_col = kernel.cols() - kernel.cols()/2-1;
		int kern_spread_row = kernel.rows() - kernel.rows()/2-1;
		int rows = data.rows();
		int cols = data.cols();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				double sum = 0;
				for (int ii = -kern_spread_row; ii <= kern_spread_row; ++ii) {
					for (int jj = -kern_spread_col; jj <= kern_spread_col; ++jj) {
						int x = jj+j;
						int y = ii+i;
						if (x < 0) x = -x;
						else if (x >= cols) x = 2*cols - x-1;
						if (y < 0) y = -y;
						else if (y >= rows) y = 2*rows - y-1;
						sum += kernel(ii+kern_spread_row,jj+kern_spread_col) * data(y,x);
					}
				}
				temp(i,j) = sum;
			}
		}
		// copy temporary matrix to data
		data = temp;
	}

};

