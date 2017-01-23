#! /usr/bin/env python

import os
import io
from distutils.core import setup, Extension
# from setuptools import setup, Extension
import numpy as np

# NOTE: to compile, run in the current directory
# python setup.py build_ext --inplace
# python setup.py develop

try:
    from Cython.Build import cythonize
    # from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True


def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)


cmdclass = {}
ext_modules = []

here = os.path.abspath(os.path.dirname(__file__))
sourcefiles = ['sparse_linear_binning/sparse_linear_binning_impl.cpp',
               'sparse_linear_binning/sparse_linear_binning.pyx']
include_path = [np.get_include(), 'sparse_linear_binning',
                'sparse_linear_binning/sparsepp']

if use_cython:
    extensions = \
        Extension("sparse_linear_binning.sparse_linear_binning",
                  sources=sourcefiles,
                  include_dirs=include_path,
                  language='c++')
else:
    sourcefiles[1] = sourcefiles[1].replace('.pyx', '.cpp')
    extensions = \
        Extension("sparse_linear_binning.sparse_linear_binning",
                  sources=sourcefiles,
                  include_dirs=include_path,
                  language='c++')


if use_cython:
    ext_modules += \
        cythonize(extensions)
else:
    ext_modules += [extensions]

long_description = read('README.rst', 'CHANGES.txt')

setup(
    name="sparse_linear_binning",
    version='1.0.0',
    url='https://github.com/jhetherly/sparse_linear_binning',
    license='MIT',
    author='Jeff Hetherly',
    author_email='jeffrey.hetherly@gmail.com',
    platforms='any',
    description='Python function for performing a linear ' +
                'binning that is optimized for sparsity',
    long_description=long_description,
    install_requires=['numpy'],
    packages=['sparse_linear_binning'],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    include_dirs=[np.get_include(), 'sparse_linear_binning',
                  'sparse_linear_binning/sparsepp']
)
