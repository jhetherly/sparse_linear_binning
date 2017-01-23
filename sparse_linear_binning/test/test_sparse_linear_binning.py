from sparse_linear_binning import sparse_linear_binning
import numpy as np
import logging
from timeit import default_timer as timer

logging.basicConfig(level=logging.INFO)

def generate_data(n_samples=100000, D=2):
    sample_coords = np.random.random(size=(n_samples, D))
    sample_weights = np.random.random(size=n_samples)
    extents = np.tile([0., 1.], D).reshape((D, 2))
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
