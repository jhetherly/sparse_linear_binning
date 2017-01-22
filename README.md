
# sparse_linear_binning: performs a linear binning (optimized for sparsity)

Performs a linear binning technique described in [Wand and Jones](https://www.crcpress.com/Kernel-Smoothing/Wand-Jones/p/book/9780412552700).
The
[asymptotic behavior](http://www.tandfonline.com/doi/abs/10.1080/00949658308810650)
of this binning technique performs better than so-called
simple binning (i.e. as in histograms).

For example, within a 2D grid with corners A, B, C, and D and a 2D point P with
weight wP:

    A-----------------------------------B
    |        |                          |
    |                                   |
    |        |                          |
    |- - - - P- - - - - - - - - - - - - |
    |        |                          |
    D-----------------------------------C

* Assign a weight to corner A of the proportion of area (times wP) between P and C
* Assign a weight to corner B of the proportion of area (times wP) between P and D
* Assign a weight to corner C of the proportion of area (times wP) between P and A
* Assign a weight to corner D of the proportion of area (times wP) between P and B

Note that the under- and overflow bins need to be accounted for when specifying
the number of bins.
For instance, if you want grid points in steps of 0.1 in a range of \[0,1\]
(i.e. (0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1), specify the number of bins to be 11.
Internally, the grid points are stored in a high performance, C++-based hash
table ([sparsepp](https://github.com/greg7mdp/sparsepp)).
This allows for finer binning because it never allocates memory for grid points
with near-zero weight.
To accommodate arbitrary numbers of bins along each dimension, an arbitrary
precision numeric library ([boost multiprecision](http://www.boost.org/doc/libs/1_63_0/libs/multiprecision/doc/html/boost_multiprecision/intro.html))
may be used internally and will negatively impact performance.
If this degradation in performance is unacceptable, consider reducing the number
of bins in such a way that the product of bin sizes is less than the numeric
maximum of "unsigned long" or "unsigned long long" on your system.

## Quickstart

## Dependencies

* numpy
