//
// Gaussian Blur Generator
//
// Example Usage: ./GaussianBlurGen 7 0.84089642
//
// Sigma 1.8, Kernel: 9
//
// References:
//   https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
//   http://dev.theomader.com/gaussian-kernel-calculator/
//

#include <cmath>    // exp
#include <cstdio>   // printf
#include <cstdlib>  // atoi
#include <memory>   // unique_ptr

int main(int argc, char* argv[])
{
  static constexpr double k_TwoPI = 6.28318530718;

  if (argc < 3)
  {
    std::printf("usage: %s <blur-size> <blur-variance>\n", argv[0]);
    return 1;
  }

  const int blur_size = std::atoi(argv[1]);

  if (blur_size <= 0)
  {
    std::printf("ERROR: <blur-size> should be greater than 0 not `%i`.\n", blur_size);
    return 2;
  }

  if ((blur_size & 1) == 0)
  {
    std::printf("ERROR: <blur-size> must be an odd number not `%i`.\n", blur_size);
    return 2;
  }

  const double variance        = std::strtod(argv[2], nullptr);
  const double variance_sq     = variance * variance;
  const double two_variance_sq = 2.0 * variance_sq;
  const double denom           = k_TwoPI * variance_sq;
  const int    num_elements    = blur_size * blur_size;
  const int    half_size       = blur_size / 2;

  const std::unique_ptr<double[]> blur_kernel = std::make_unique<double[]>(num_elements);
  double                          total_value = 0.0f;
  int                             kernel_y    = 0;

  for (int y = -half_size; y <= half_size; ++y)
  {
    int kernel_x = 0;

    for (int x = -half_size; x <= half_size; ++x)
    {
      const double x_sq      = double(x * x);
      const double y_sq      = double(y * y);
      const double numerator = std::exp(-(x_sq + y_sq) / two_variance_sq);
      const double value     = numerator / denom;

      blur_kernel[kernel_x + kernel_y * blur_size] = value;
      total_value += value;

      ++kernel_x;
    }

    ++kernel_y;
  }

  const auto printMatrix = [blur_size](const std::unique_ptr<double[]>& matrix) {
    for (int y = 0; y < blur_size; ++y)
    {
      std::printf("| ");

      for (int x = 0; x < blur_size; ++x)
      {
        std::printf("%.8f ", matrix[x + y * blur_size]);
      }

      std::printf("|\n");
    }
  };

  const std::unique_ptr<double[]> normalized_blur_kernel = std::make_unique<double[]>(num_elements);

  for (int i = 0; i < num_elements; ++i)
  {
    normalized_blur_kernel[i] = blur_kernel[i] / total_value;
  }

  std::printf("Unnormalized Matrix(Total:%f):\n", total_value);
  printMatrix(blur_kernel);

  std::printf("\nNormalized Matrix:\n");
  printMatrix(normalized_blur_kernel);

  std::printf("\nWeights:\n");

  const double denom1D = std::sqrt(k_TwoPI * variance);
  int          xx      = -half_size;

  for (int i = 0; i < blur_size; ++i)
  {
    const double x_sq           = double(xx * xx);
    const double numerator      = std::exp(-x_sq / two_variance_sq);
    const double one_d_gaussian = numerator / denom1D;

    std::printf("[%3i] = %.15f => %.15f (%.15f)\n",
                xx,
                normalized_blur_kernel[i + i * blur_size],
                std::sqrt(normalized_blur_kernel[i + i * blur_size]),
                one_d_gaussian);

    ++xx;
  }

  return 0;
}
