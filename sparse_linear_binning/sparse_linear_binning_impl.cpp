// #include <ctime>
// #include <cstdlib>
#include <vector>
// #include <iostream>

#include "boost/multiprecision/cpp_int.hpp"

#include "sparse_linear_binning_impl.hpp"

using boost::multiprecision::cpp_int;
using std::vector;


// Driver function for the linear binning method
void sparse_linear_binning_impl (const unsigned long D, // space dimension
                                 const double *X, // n_samples x D
                                 const double *weights, // n_samples
                                 const unsigned long n_samples,
                                 const double *extents, // D x 2
                                 const double *bin_sizes, // D
                                 const unsigned long *sizes, // D
                                 double *&result_X, // grid points x D
                                 double *&result_weights,  // grid points
                                 unsigned long &n_result) // # of grid points
{
  cpp_int prod = 1;
  for (unsigned long i = 0; i < D; ++i)
    prod *= sizes[i];

  vector<unsigned long> sizes_to_pass(D);
  vector<double>   extents_to_pass(D*2),
                   bin_sizes_to_pass(D);

  for (unsigned long i = 0; i < D; ++i) {
    sizes_to_pass [i] = sizes[i];
    extents_to_pass[2*i] = extents[2*i];
    extents_to_pass[2*i + 1] = extents[2*i + 1];
    bin_sizes_to_pass[i] = bin_sizes[i];
  }

  if (prod > std::numeric_limits<BigFiniteNum_t>::max()) {
    if (static_cast<unsigned long>(std::numeric_limits<FiniteNum_t>::digits) > D)
      make_linear_binning<arbitrary_precision_global_indices_D_fin<FiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
    else if (static_cast<unsigned long>(std::numeric_limits<BigFiniteNum_t>::digits) > D)
      make_linear_binning<arbitrary_precision_global_indices_D_fin<BigFiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
    else
      make_linear_binning<arbitrary_precision_global_indices_D_arb>(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
  }
  else if (prod > std::numeric_limits<FiniteNum_t>::max()) {
    if (static_cast<unsigned long>(std::numeric_limits<FiniteNum_t>::digits) > D)
      make_linear_binning<finite_global_indices_D_fin<BigFiniteNum_t, FiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
    else if (static_cast<unsigned long>(std::numeric_limits<BigFiniteNum_t>::digits) > D)
      make_linear_binning<finite_global_indices_D_fin<BigFiniteNum_t, BigFiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
    else
      make_linear_binning<finite_global_indices_D_arb<BigFiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
  }
  else {
    if (static_cast<unsigned long>(std::numeric_limits<FiniteNum_t>::digits) > D)
      make_linear_binning<finite_global_indices_D_fin<FiniteNum_t, FiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
    else if (static_cast<unsigned long>(std::numeric_limits<BigFiniteNum_t>::digits) > D)
      make_linear_binning<finite_global_indices_D_fin<FiniteNum_t, BigFiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
    else
      make_linear_binning<finite_global_indices_D_arb<FiniteNum_t> >(
                                  D, X, weights, n_samples,
                                  extents_to_pass, bin_sizes_to_pass,
                                  sizes_to_pass,
                                  result_X, result_weights, n_result);
  }

}


// using namespace std;
//
// int main(int /*argc*/, char * /*argv*/[])
// {
//   static const unsigned long D = 2;
//
//   const bool      print_points = false;
//   unsigned long   n_samples = 100000;
//   vector<unsigned long> sizes(D);
//   vector<double>   X(n_samples*D),
//                    weights(n_samples),
//                    extents(D*2),
//                    bin_sizes(D);
//   double          *result_X = NULL,
//                   *result_weights = NULL;
//   unsigned long    n_result;
//   double           weight_sum = 0.0,
//                    original_weight_sum = 0.0;
//
//   cout << "dimension: " << D << endl;
//   cout << "samples:   " << n_samples << endl;
//   srand(1);
//   if (print_points) cout << "\nsource points:" << endl;
//   for (unsigned long i = 0; i < n_samples; ++i) {
//     for (unsigned long j = 0; j < D; ++j) {
//       double num = static_cast<double>(rand() % 100000)/100000.;
//       X[i*D + j] = num;
//       if (print_points) {
//         cout << num;
//         if (j + 1 < D) cout << ", ";
//       }
//     }
//     double num = static_cast<double>(rand() % 100000)/100000.;
//     weights[i] = num;
//     original_weight_sum += num;
//     if (print_points) cout << ": " << num << endl;
//   }
//   if (print_points) cout << endl;
//
//   for (unsigned long i = 0; i < D; ++i) {
//     sizes[i] = 51;
//     extents[2*i] = 0.000001;
//     extents[2*i + 1] = 0.999999;
//     bin_sizes[i] = (extents[2*i + 1] - extents[2*i])/static_cast<double>(sizes[i] - 1);
//   }
//
//   std::clock_t start;
//   double duration;
//   start = std::clock();
//
//   if (D == 2) {
//     cout << "copula" << endl;
//     make_2D_linear_binning_for_copula(&X[0], &weights[0], n_samples,
//                                       &sizes[0],
//                                       result_X, result_weights, n_result);
//   }
//   else {
//     // cout << "static" << endl;
//     // sparse_linear_binning<D>(&X[0], &weights[0], n_samples,
//     //                          &extents[0], &bin_sizes[0], &sizes[0],
//     //                          result_X, result_weights, n_result);
//     cout << "dynamic" << endl;
//     sparse_linear_binning(D, &X[0], &weights[0], n_samples,
//                           &extents[0], &bin_sizes[0], &sizes[0],
//                           result_X, result_weights, n_result);
//   }
//
//   duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
//   cout<<"creating the linear bins took "<< duration << "s." << endl;
//   cout<<"number of grid points "<< n_result << endl;
//   if (print_points) cout << "\ngrid points:" << endl;
//   for (size_t i = 0; i < n_result; ++i) {
//     if (print_points) {
//       cout << i << ": ";
//       for (unsigned long j = 0; j < D; ++j) {
//         cout << result_X[i*D + j];
//         if (j + 1 < D) cout << ", ";
//       }
//       cout << ": " << result_weights[i] << endl;
//     }
//     weight_sum += result_weights[i];
//   }
//   if (print_points) cout << endl;
//   cout << "original sum of weights = " << original_weight_sum << endl;
//   cout << "sum of weights =          " << weight_sum << endl;
//   delete[] result_X;
//   delete[] result_weights;
//
//   return 0;
// }
