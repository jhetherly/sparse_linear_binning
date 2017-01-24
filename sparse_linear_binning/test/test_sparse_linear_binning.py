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
    sizes = np.full(D, 51, dtype=np.dtype('uint64'))

    return sample_coords, sample_weights, extents, sizes


def test_sum_of_weights():
    sample_coords, sample_weights, extents, sizes = generate_data(1000000)
    start = timer()
    coords, weights = sparse_linear_binning(sample_coords, sample_weights,
                                            extents, sizes)
    end = timer()

    logging.info('One million 2D points binned in {}s'.format(end - start))
    assert np.allclose(weights.sum(), sample_weights.sum())

    # sample_coords, sample_weights, extents, sizes = generate_data(2)
    # D = 2
    # extents = np.tile([-20., 0.], D).reshape((D, 2))
    # sample_coords = 40.*sample_coords - 20.
    # coords, weights = sparse_linear_binning(sample_coords, sample_weights,
    #                                         extents, sizes)
    # print(sample_coords, sample_weights)
    # print(coords, weights)
    # print(sample_weights.sum(), weights.sum())
