from sparse_linear_binning import sparse_linear_binning
import numpy as np


def generate_data(n_samples=100000, D=2):
    sample_coords = np.random.random(size=(n_samples, D))
    sample_weights = np.random.random(size=n_samples)
    extents = np.tile([0., 1.], D).reshape((D, 2))
    sizes = np.full(D, 51, dtype=np.dtype('uint64'))

    return sample_coords, sample_weights, extents, sizes


def test_sum_of_weights():
    sample_coords, sample_weights, extents, sizes = generate_data(1000000)
    coords, weights = sparse_linear_binning(sample_coords, sample_weights,
                                            extents, sizes)

    assert np.allclose(weights.sum(), sample_weights.sum())
