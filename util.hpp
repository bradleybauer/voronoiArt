#include <eigen3/Eigen/Eigen>
#include <png++/png.hpp>
#include <cmath>
#include <random>
#include <ctime>
using namespace Eigen;
class Process {
public:
	Process(){};

	void gaussianConvolve(int kernel_width, MatrixXd& data) {
		kernel_width += kernel_width % 2 == 0; 
		// compute gaussian kernel
		MatrixXd kernel(kernel_width,kernel_width);
		for (int i = 0; i < kernel_width; ++i) {
			for (int j = 0; j < kernel_width; ++j) {
				// Generate coordinates
				Vector2d p = {i, j};
				p /= double(kernel_width-1);
				p = p.array() * 2. - 1.;
				p = p.array() * p.array();
				p*=6.;

				const double sigmaSqrd = 20;
				const double pi = 3.14159268;
				const double X = p.sum();
				kernel(i,j) = 1./(2.*pi*sigmaSqrd)*std::exp(-X/(2.*sigmaSqrd));
			}
		}
		convolve_mirror_pad(data, kernel);
	}

	void sobel(MatrixXd& data) {
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
				const double r = .33*image[i][j].red/65535.;
				const double g = .58*image[i][j].green/65535.;
				const double b = .11*image[i][j].blue/65535.;
				data(i,j) = r + g + b;
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
	//     However, to compute the list of lists requires about .5 * Num_Vertices^2 loops?
	//     That might make this optimization less efficient
	//
	//     I think the coloring computation could benefit greatly from running on the graphics card
	void computeVertices(Matrix<int, Dynamic, 2>& vertices,
	                     MatrixXd& data, const int max_vertices) {
		auto eng = std::default_random_engine(std::time(0));
		// A larger minimum value here causes larger sections of the image to be
		// covered by a single voronoi cell.
		const double minimum = .04;
		const double maximum = 1.;
		auto dist = std::uniform_real_distribution<double>(minimum, maximum);
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
		auto dist2 = std::uniform_real_distribution<double>(0,1);
		int rows = vertices.rows();
		while (rows > max_vertices) {
			int rowToRemove = int(dist2(eng) * vertices.rows());
			rows--;
			if( rowToRemove < rows )
					vertices.block(rowToRemove, 0, rows-rowToRemove, 2) =
						vertices.block(rowToRemove+1, 0, rows-rowToRemove,2);
			vertices.conservativeResize(rows, 2);
		}
	}

	void colorImageWithVertices(png::image<png::rgb_pixel_16>& image,
	                            Matrix<int, Dynamic, 2>& vertices) {
		// loop through image pixels
		//    find vertex that is closest to pixel
		//    set pixel to average color of image around point vertex
		const int height = image.get_height();
		const int width = image.get_width();
		int indx_of_min;
		double d, dd, x, y;
		for (int i = 0; i < height; ++i) {
			printf("Progress: %f\r", i/float(height)); fflush(stdout);
			for (int j = 0; j < width; ++j) {
				indx_of_min = 0;
				d = std::numeric_limits<double>::infinity();
				for (int k = 0; k < vertices.rows(); ++k) {
					x = j - vertices(k, 1);
					y = i - vertices(k, 0);
					dd = x*x + y*y;
					// dd = std::fabs(x) + std::fabs(y);
					if (dd < d) {
						indx_of_min = k;
						d = dd;
					}
				}
				image[i][j] = image[vertices(indx_of_min,0)][vertices(indx_of_min,1)];
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
		double sum = 0.;
		int x, y;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				for (int ii = -kern_spread_row; ii <= kern_spread_row; ++ii) {
					for (int jj = -kern_spread_col; jj <= kern_spread_col; ++jj) {
						x = fabs(jj+j);
						y = fabs(ii+i);
						if (x >= cols) x = 2*cols - x-1;
						if (y >= rows) y = 2*rows - y-1;
						sum += kernel(ii+kern_spread_row,jj+kern_spread_col) * data(y,x);
					}
				}
				temp(i,j) = sum;
				sum = 0;
			}
		}
		// copy temporary matrix to data
		data = temp;
	}

};

