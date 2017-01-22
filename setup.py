from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

# NOTE: to compile, run in the current directory
# python setup.py build_ext --inplace

ext_modules = \
    cythonize('sparse_linear_binning/sparse_linear_binning.pyx')

setup(
    name="sparse_linear_binning",
    ext_modules=ext_modules,
    include_dirs=[np.get_include(), 'sparse_linear_binning',
                  'sparse_linear_binning/sparsepp']
)
