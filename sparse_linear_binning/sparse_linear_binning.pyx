# distutils: language = c++
# distutils: sources = sparse_linear_binning/sparse_linear_binning_impl.cpp
# distutils: libraries = ["m"]

import numpy as np
cimport numpy as np
#from libcpp cimport double

np.import_array()

ctypedef np.int32_t DTYPE_t

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef data_to_numpy_array_with_spec(void * ptr, np.npy_intp N):
    cdef np.ndarray[double, ndim=1] arr = np.PyArray_SimpleNewFromData(
                                                    1, &N, np.NPY_DOUBLE, ptr)
    # NOTE: transfer ownership of array
    PyArray_ENABLEFLAGS(arr, np.NPY_OWNDATA)
    return arr

cdef extern from "sparse_linear_binning_impl.hpp":
    void sparse_linear_binning_impl(const unsigned long D,
                                    const double *X,
                                    const double *weights,
                                    const unsigned long n_samples,
                                    const double *extents,
                                    const double *bin_sizes,
                                    const unsigned long *sizes,
                                    double *&result_X,
                                    double *&result_weights,
                                    unsigned long &n_result)

def sparse_linear_binning (np.ndarray[double, ndim=2, mode='c'] X,
                           np.ndarray[double, ndim=1, mode='c'] weights,
                           np.ndarray[double, ndim=2, mode='c'] extents,
                           np.ndarray[unsigned long, ndim=1, mode='c'] sizes):
    """Perform a linear binning (optimized for sparsity)

    This algorithm is able to opperate on N data point in D dimensions and
    produces a variable number of grid points, G, with their corresponding
    weights.
    The input arrays need to be c-contiguous and strictly adhear to the
    required types (double and unsigned long).

    Parameters
    ----------
    X : 2D numpy array with shape (N, D) of type double
        data coordinates
    weights : 1D numpy array with shape (N) of type double
        data weights
    extents : 2D numpy array with shape (D, 2) of type double
        limits of grid (all data outside this retangular region is under- or
        overflow)
    sizes : 1D numpy array with shape (D) of type unsigned long
        number of bin centers (grid points) in each dimension

    Returns
    -------
    2D numpy array with shape (G, D) of type double
        grid points with non-zero weights
    1D numpy array with shape (N) of type double
        weights at the corresponding grid points
    """
    assert X.shape[0] == weights.size   # these must be the number of data points
    assert X.shape[1] == extents.shape[0] == sizes.size  # these must be the dimension
    assert extents.shape[1] == 2

    cdef unsigned long n_samples = X.shape[0]
    cdef unsigned long D = X.shape[1]
    cdef double* result_X = NULL
    cdef double* result_weights = NULL
    cdef unsigned long n_result = 0
    cdef np.ndarray[double, ndim = 1, mode = 'c'] bin_sizes = \
        np.empty(D, np.dtype('d'))

    bin_sizes = (np.abs(np.diff(extents, axis=1).T)/(sizes - 1)).flatten(order='C')

    sparse_linear_binning_impl(D, &X[0, 0], &weights[0], n_samples,
                               &extents[0, 0], &bin_sizes[0], &sizes[0],
                               result_X, result_weights, n_result)

    result_x = data_to_numpy_array_with_spec(result_X, D*n_result)
    result_w = data_to_numpy_array_with_spec(result_weights, n_result)

    return result_x.reshape((n_result, D)), result_w
