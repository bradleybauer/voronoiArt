#include <eigen3/Eigen/Eigen>
#include <png++/png.hpp>
#include <cmath>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iterator>
#include <utility>
using namespace Eigen;
class Process {
public:
    void gaussianConvolve(int kernel_width, MatrixXd& data) {
        kernel_width += kernel_width % 2 == 0;
        // compute gaussian kernel
        MatrixXd kernel(kernel_width, kernel_width);
        for (int i = 0; i < kernel_width; ++i) {
            for (int j = 0; j < kernel_width; ++j) {
                // Generate coordinates
                Vector2d p = { i, j };
                p /= double(kernel_width - 1);
                p = p.array() * 2. - 1.;
                p = p.array() * p.array();
                p *= 6.;

                const double sigmaSqrd = 20;
                const double pi = 3.14159268;
                const double X = p.sum();
                kernel(i, j) = 1. / (2.*pi*sigmaSqrd)*std::exp(-X / (2.*sigmaSqrd));
            }
        }
        convolve_mirror_pad(data, kernel);
    }

    void sobel(MatrixXd& data) {
        MatrixXd D(3, 3);
        D << -1, 0, 1,
             -2, 0, 2,
             -1, 0, 1;
        MatrixXd temp = data;
        convolve_mirror_pad(temp, D);
        D = D.transpose().eval();
        convolve_mirror_pad(data, D);
        data = Eigen::sqrt(data.array().pow(2) + temp.array().pow(2)).eval();
    }

    void intensity(png::image<png::rgb_pixel_16>& image, MatrixXd& data) {
        for (int i = 0; i < data.rows(); ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                const double r = .33*image[i][j].red / 65535.;
                const double g = .58*image[i][j].green / 65535.;
                const double b = .11*image[i][j].blue / 65535.;
                data(i, j) = r + g + b;
            }
        }
    }

    void computeVertices(Matrix<long, Dynamic, 6>& vertices,
        MatrixXd& data, int max_vertices) {
        auto eng = std::default_random_engine(std::time(0));
        // A larger minimum value here causes larger sections of the image to be
        // covered by a single voronoi cell.
        const double minimum = .04;
        const double maximum = 1.;
        auto dist = std::uniform_real_distribution<double>(minimum, maximum);
        int count = 0;
        for (int i = 0; i < data.rows(); ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                if (dist(eng) < data(i, j)) {
                    vertices(count, 0) = i;
                    vertices(count, 1) = j;
                    count++;
                }
            }
        }

        max_vertices = std::min(count, max_vertices);

        std::cout << "Reducing number of vertices" << std::endl;
        // Reduce the number of vertices to max_vertices
        vertices.conservativeResize(count, NoChange);
        auto dist2 = std::uniform_int_distribution<int>(0, vertices.rows()-1);
        std::vector<int> verts_to_keep;
        std::generate_n(std::back_inserter(verts_to_keep), max_vertices, [&]() { return dist2(eng); });
        std::sort(verts_to_keep.begin(), verts_to_keep.end());
        const int cols = vertices.cols();
        for (int i = 0; i < max_vertices; ++i) {
            vertices.block(i, 0, 1, cols) = vertices.block(verts_to_keep[i], 0, 1, cols);
        }
        vertices.conservativeResize(max_vertices, NoChange);
    }

    void colorImageWithVertices(png::image<png::rgb_pixel_16>& image,
        Matrix<long, Dynamic, 6>& vertices) {
        // loop through image pixels
        //    find vertex that is closest to pixel
        //    set pixel to average color of image around point vertex
        std::vector<std::vector<std::pair<int, int>>> vertex_to_pixel_map;
        vertex_to_pixel_map.resize(vertices.rows());
        const int height = image.get_height();
        const int width = image.get_width();
        int indx_of_min;
        double d, dd, x, y;
        for (int i = 0; i < height; ++i) {
            printf("Coloring %f\r", i / float(height)); fflush(stdout);
            for (int j = 0; j < width; ++j) {

                // Find vertex of min distance to pixel at {i,j}
                indx_of_min = 0;
                d = std::numeric_limits<double>::infinity();
                for (int k = 0; k < vertices.rows(); ++k) {
                    y = i - vertices(k, 0);
                    x = j - vertices(k, 1);
                    dd = x * x + y * y;
                    // dd = std::fabs(x) + std::fabs(y);
                    if (dd < d) {
                        indx_of_min = k;
                        d = dd;
                    }
                }

                // Add pixel_ij color to vertex and increment count of pixels beloning
                // to vertex
                vertex_to_pixel_map[indx_of_min].push_back({ i,j });
                png::rgb_pixel_16 pix = image[i][j];
                vertices(indx_of_min, 2) += pix.red;
                vertices(indx_of_min, 3) += pix.green;
                vertices(indx_of_min, 4) += pix.blue;
                vertices(indx_of_min, 5) += 1;
            }
        }

        // Give pixel_ij the color of vertex_k
        for (int k = 0; k < vertex_to_pixel_map.size(); ++k) {
            for (int p = 0; p < vertex_to_pixel_map[k].size(); ++p) {
                const int i = vertex_to_pixel_map[k][p].first;
                const int j = vertex_to_pixel_map[k][p].second;
                const double r = vertices(k, 2) / double(vertices(k, 5));
                const double g = vertices(k, 3) / double(vertices(k, 5));
                const double b = vertices(k, 4) / double(vertices(k, 5));
                image[i][j] = png::rgb_pixel_16(r, g, b);
            }
        }
        std::cout << "Coloring 1.000000" << std::endl;
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
        const int kern_spread_col = kernel.cols() - kernel.cols() / 2 - 1;
        const int kern_spread_row = kernel.rows() - kernel.rows() / 2 - 1;
        const int rows = data.rows();
        const int cols = data.cols();
        double sum = 0.;
        int x, y;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                for (int ii = -kern_spread_row; ii <= kern_spread_row; ++ii) {
                    for (int jj = -kern_spread_col; jj <= kern_spread_col; ++jj) {
                        x = fabs(jj + j);
                        y = fabs(ii + i);
                        if (x >= cols) x = 2 * cols - x - 1;
                        if (y >= rows) y = 2 * rows - y - 1;
                        sum += kernel(ii + kern_spread_row, jj + kern_spread_col) * data(y, x);
                    }
                }
                temp(i, j) = sum;
                sum = 0;
            }
        }
        // copy temporary matrix to data
        data = temp;
    }

};

