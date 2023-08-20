//
// Gaussian Blur Generator
// 
// Example Usage: ./GaussianBlurGen 9 1.7573
//
// References:
//   https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm
//

#include <cmath>    // exp, acos, sqrt
#include <cstdio>   // printf
#include <cstdlib>  // atoi
#include <memory>   // unique_ptr

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    std::printf("usage: %s <blur-radius> <blur-variance>\n", argv[0]);
    return 1;
  }

  const int blur_radius = std::atoi(argv[1]);

  if (blur_radius <= 0)
  {
    std::printf("ERROR: <blur-radius> should be greater than 0 not `%i`.\n", blur_radius);
    return 2;
  }

  const double                    variance               = std::strtod(argv[2], nullptr);
  const int                       blur_size              = blur_radius * 2 - 1;
  const double                    two_pi                 = 2.0 * std::acos(-1.0);
  const double                    variance_sq            = variance * variance;
  const double                    two_variance_sq        = 2.0 * variance_sq;
  const double                    denom2D                = two_pi * variance_sq;
  const double                    denom1D                = std::sqrt(two_pi) * variance;
  const int                       num_elements           = blur_size * blur_size;
  const int                       half_size              = blur_size / 2;
  const std::unique_ptr<double[]> kernel_memory          = std::make_unique<double[]>(num_elements + num_elements + blur_size);
  double* const                   blur_kernel            = kernel_memory.get();
  double* const                   normalized_blur_kernel = blur_kernel + num_elements;
  double* const                   oneD_kernel            = normalized_blur_kernel + num_elements;
  double                          total_value            = 0.0f;

  int kernel_y = 0;
  for (int y = -half_size; y <= half_size; ++y)
  {
    int kernel_x = 0;

    for (int x = -half_size; x <= half_size; ++x)
    {
      const double x_sq      = double(x * x);
      const double y_sq      = double(y * y);
      const double numerator = std::exp(-(x_sq + y_sq) / two_variance_sq);
      const double value     = numerator / denom2D;

      blur_kernel[kernel_x + kernel_y * blur_size] = value;
      total_value += value;

      ++kernel_x;
    }

    ++kernel_y;
  }

  double normalized_total = 0.0;
  for (int i = 0; i < num_elements; ++i)
  {
    normalized_blur_kernel[i] = blur_kernel[i] / total_value;
    normalized_total += normalized_blur_kernel[i];
  }

  double oneDTotal = 0.0;
  int xx = -half_size;
  for (int i = 0; i < blur_size; ++i)
  {
    const double x_sq           = double(xx * xx);
    const double numerator      = std::exp(-x_sq / two_variance_sq);
    const double one_d_gaussian = numerator / denom1D;

    oneD_kernel[i] = one_d_gaussian;
    oneDTotal += one_d_gaussian;

    ++xx;
  }

  const auto PrintMatrix2D = [&](const char* title, const double total, const double* const matrix) {

    std::printf("%s(%ix%i, sigma = %f)(Total:%f):\n", title, blur_size, blur_size, variance, total);
    for (int y = 0; y < blur_size; ++y)
    {
      std::printf("  | ");

      for (int x = 0; x < blur_size; ++x)
      {
        std::printf("%.8f ", matrix[x + y * blur_size]);
      }

      std::printf("|\n");
    }
    };

  PrintMatrix2D("Gaussian2D", total_value, blur_kernel);
  std::printf("\n");
  PrintMatrix2D("Gaussian2D_normalized", normalized_total, normalized_blur_kernel);
  std::printf("\n");

  std::printf("Gaussian1D(Total:%f):\n", oneDTotal);
  {
    int xx = -half_size;

    for (int i = 0; i < blur_size; ++i)
    {
      const int index = i + i * blur_size;

      std::printf("  [%3i] = {orig(%.15f), orig_sqrt(%.15f), norm(%.15f), norm_sqrt(sqrt = %.15f), 1d(%.15f), 1d_norm(%.15f) }\n", 
        xx, 
        blur_kernel[index], 
        std::sqrt(blur_kernel[index]),
        normalized_blur_kernel[index], 
        std::sqrt(normalized_blur_kernel[index]), 
        oneD_kernel[i],
        oneD_kernel[i] / oneDTotal);
      ++xx;
    }

  }
  return 0;
}
