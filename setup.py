from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import os.path

curdir = os.path.dirname(os.path.realpath(__file__))
ext_modules = [ Extension(
    name="adaptive_anderson_solver",
    sources=["adaptive_anderson_solver.pyx"],
    libraries=["adaptive_anderson_solver"] + np.__config__.blas_opt_info['libraries'],
    include_dirs=[np.get_include(), os.path.dirname(__file__)],
    library_dirs=["./lib/"] + np.__config__.blas_opt_info['library_dirs'],
    extra_compile_args=["-O3"],
    extra_link_args=["-Wl,-rpath,./lib:"+curdir]
)]

setup(
    ext_modules=cythonize(ext_modules)
)
