|Build Status|

`sparse\_linear\_binning: linear binning (optimized for sparsity) <https://github.com/jhetherly/sparse_linear_binning>`__
=========================================================================================================================

Performs a linear binning technique described in `Wand and
Jones <https://www.crcpress.com/Kernel-Smoothing/Wand-Jones/p/book/9780412552700>`__
on a regularly-spaced grid in an arbitrary number of dimensions. The
`asymptotic
behavior <http://www.tandfonline.com/doi/abs/10.1080/00949658308810650>`__
of this binning technique performs better than so-called simple binning
(i.e. as in histograms). Each data point in ``d``-dimensional space must
have an associated weight (for equally weighted points just use a weight
of ``1.0`` for each point).

For example, within a (segment of a) 2D grid with corners A, B, C, and D
and a 2D data point P with weight wP:

::

    A-----------------------------------B
    |        |                          |
    |                                   |
    |        |                          |
    |- - - - P- - - - - - - - - - - - - |
    |        |                          |
    D-----------------------------------C

-  Assign a weight to corner A of the proportion of area between P and C
   (times wP)
-  Assign a weight to corner B of the proportion of area between P and D
   (times wP)
-  Assign a weight to corner C of the proportion of area between P and A
   (times wP)
-  Assign a weight to corner D of the proportion of area between P and B
   (times wP)

Note that the under- and overflow bins need to be accounted for when
specifying the numbers of grid points in each dimension (grid points act
as bin centers). For instance, if you want grid points in steps of 0.1
in a range of [0,1] (i.e. (0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1)),
specify the number of grid points to be 11. Internally, the grid points
are stored in a high performance, C++-based hash table
(`sparsepp <https://github.com/greg7mdp/sparsepp>`__). This allows for
finer binning in some circumstances because the hash table doesn't
allocates memory for grid points with near-zero weight. To accommodate
arbitrary numbers of bins along each dimension, an arbitrary precision
numeric library (`boost
multiprecision <http://www.boost.org/doc/libs/1_63_0/libs/multiprecision/doc/html/boost_multiprecision/intro.html>`__)
may be used internally and will negatively impact performance. If this
degradation in performance is unacceptable, consider reducing the number
of grid points in such a way that the product of grid points in all
dimensions is less than the numeric maximum of "unsigned long long" on
your system. For instance, in 20 dimenisons with each dimension having
51 grid points gives a total of 14171098670753043575626125424226001
potential grid points at which point the arbitrary precision library
must take care of all arithmetic related to determining grid points.

Quickstart
----------

-  pip install sparse\_linear\_binning

or

-  git clone https://github.com/jhetherly/sparse\_linear\_binning
-  cd sparse\_linear\_binning
-  python setup.py install

Example
-------

This constructs one million random 2D points in the unit square with
random weights and constructs a grid of ``51`` by ``51`` (can be
different along different dimensions) linearly binned "bin centers." The
boundaries of the grid of bin centers are specified by ``extents`` and
can be thought of as the under- and overflow bins.

.. code:: python

    from sparse_linear_binning import sparse_linear_binning
    import numpy as np

    # generate one million random 2D points and weights
    # (should take less than a second to bin)
    n_samples=1000000
    D=2

    # coordinates, weights, and extents must be of type "double"
    sample_coords = np.random.random(size=(n_samples, D))
    sample_weights = np.random.random(size=n_samples)
    extents = np.tile([0., 1.], D).reshape((D, 2))
    # n_bins must be of type "unsigned long"
    n_bins = np.full(D, 51, dtype=np.dtype('uint64'))

    coords, weights = sparse_linear_binning(sample_coords, sample_weights,
                                            extents, n_bins)

    # check that weights on grid match original weights
    print(np.allclose(weights.sum(), sample_weights.sum()))

Dependencies
------------

-  numpy

.. |Build Status| image:: https://travis-ci.org/jhetherly/sparse_linear_binning.svg?branch=master
   :target: https://travis-ci.org/jhetherly/sparse_linear_binning
