from sparse_linear_binning import sparse_linear_binning
import numpy as np
import logging
from timeit import default_timer as timer

logging.basicConfig(level=logging.INFO)

def generate_data(n_samples=100000, D=2):
    sample_coords = np.random.random(size=(n_samples, D))
    sample_weights = np.random.random(size=n_samples)
    # NOTE: purposely limiting the range to test over- and underflow bins
    extents = np.tile([0.02, 0.8999], D).reshape((D, 2))
    sizes = np.full(D, 51)

    return sample_coords, sample_weights, extents, sizes


def test_sum_of_weights():

    # tests that the sum of weights in the binned grid is preserved
    sample_coords, sample_weights, extents, sizes = generate_data(1000000)
    start = timer()
    coords, weights = sparse_linear_binning(sample_coords, sample_weights,
                                            extents, sizes)
    end = timer()

    logging.info('\n')
    logging.info('One million 2D points binned with sparse_linear_binning in {}s'.format(end - start))
    assert np.allclose(weights.sum(), sample_weights.sum())

    x = np.ascontiguousarray(sample_coords[:,0])
    y = np.ascontiguousarray(sample_coords[:,1])
    start = timer()
    np.histogram2d(x, y,
                   weights=sample_weights,
                   bins=sizes, range=extents)
    end = timer()
    logging.info('For comparison, np.histogram2d finished in {}s'.format(end - start))


    # tests specific values on the grid
    sample_coords = np.array([[0.2, 0.9], [0.5, 1.1], [-0.1, 0.7]])
    sample_weights = np.array([25, 50, 25])
    extents = np.array([[0.0, 1.0], [0.0, 1.0]])
    sizes = np.array([11, 11])
    coords, weights = sparse_linear_binning(sample_coords, sample_weights,
                                            extents, sizes)

    pass_value_test = True
    value_tests = 0
    for i in range(coords.shape[0]):
        if np.allclose(coords[i, 0], 0.0) and np.allclose(coords[i, 1], 0.7):
            pass_value_test &= np.allclose(weights[i], 25.0)
            value_tests += 1
        elif np.allclose(coords[i, 0], 0.2) and np.allclose(coords[i, 1], 0.9):
            pass_value_test &= np.allclose(weights[i], 25.0)
            value_tests += 1
        elif np.allclose(coords[i, 0], 0.5) and np.allclose(coords[i, 1], 1.0):
            pass_value_test &= np.allclose(weights[i], 50.0)
            value_tests += 1
        else:
            pass_value_test &= np.allclose(weights[i], 0.0)

    assert pass_value_test and value_tests == 3
