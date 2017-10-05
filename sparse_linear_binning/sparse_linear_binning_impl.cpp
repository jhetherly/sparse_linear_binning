#include <vector>

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
