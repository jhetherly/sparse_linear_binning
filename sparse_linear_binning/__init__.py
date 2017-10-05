import numpy as np
from .sparse_linear_binning import cython_sparse_linear_binning


def sparse_linear_binning(X, weights=None, extents=None, sizes=None):
    """Perform a linear binning on weighted data (optimized for sparsity)

    This algorithm is able to operate on N data point in D dimensions and
    produces a variable number of grid points, G, with their corresponding
    weights.
    The input arrays need to be c-contiguous.

    Parameters
    ----------
    X : 2D numpy array with shape (N, D) of type double
        data coordinates
    weights : 1D numpy array with shape (N) of type double
        data weights (default is equally-weighted data)
    extents : 2D numpy array with shape (D, 2) of type double
        limits of grid (all data outside this rectangular region is under- or
        overflow) (default is the minimum and maximum along each dimension)
    sizes : 1D numpy array with shape (D) indicating the number of
        bin centers (grid points) in each dimension (default is 11)

    Returns
    -------
    2D numpy array with shape (G, D) of type double
        grid points with non-zero weights
    1D numpy array with shape (N) of type double
        weights at the corresponding grid points
    """

    if weights is None:
        weights = np.ones(X.shape[0], dtype=np.dtype('double'))
    if extents is None:
        extents = np.empty((X.shape[1], 2), dtype=np.dtype('double'))
        for d in range(X.shape[1]):
            extents[d, 0] = X[:, d].min()
            extents[d, 1] = X[:, d].max()
    if sizes is None:
        sizes = np.full(X.shape[1], 11, dtype=np.dtype('uint64'))

    return cython_sparse_linear_binning(X.astype(np.dtype('double')),
                                        weights.astype(np.dtype('double')),
                                        extents.astype(np.dtype('double')),
                                        sizes.astype(np.dtype('uint64')))
