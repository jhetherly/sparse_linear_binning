#! /usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

# NOTE: to compile, run in the current directory
# python setup.py build_ext --inplace

try:
    # from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = { }
ext_modules = [ ]

if use_cython:
    ext_modules += [
        Extension("sparse_linear_binning.sparse_linear_binning",
                  sources=["sparse_linear_binning/sparse_linear_binning.pyx",
                    "sparse_linear_binning/sparse_linear_binning_impl.cpp"],
                  language='c++'
                  ),
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("sparse_linear_binning.sparse_linear_binning",
                  sources=['sparse_linear_binning/sparse_linear_binning.cpp',
                    "sparse_linear_binning/sparse_linear_binning_impl.cpp"],
                  language='c++'
                  ),
    ]

long_description = read('README.md', 'CHANGES.txt')

setup(
    name="sparse_linear_binning",
    version='1.0.0',
    url='https://github.com/jhetherly/sparse_linear_binning',
    license='MIT',
    author='Jeff Hetherly',
    author_email='fredslacks@gmail.com',
    platforms='any',
    description='Python function for performing a linear binning and optimized for sparsity',
    long_description=long_description,
    # install_requires=['numpy>=10.0.0'],
    # tests_require=['pytest'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    # ext_modules=cythonize('sparse_linear_binning/sparse_linear_binning.pyx')
    include_dirs=[np.get_include(), 'sparse_linear_binning',
                  'sparse_linear_binning/sparsepp']
)
