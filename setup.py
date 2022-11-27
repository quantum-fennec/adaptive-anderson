from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os.path

try:
  import numpy.distutils.system_info as np_config
except ImportError:
  """ In my system, this gives wrong blas linking info """
  import numpy.__config__ as np_config

get_config = np_config.get_info
config = get_config('lapack_opt') or get_config('lapack_ilp64_opt')

if not config:
  raise Exception('Cannot determine LAPACK libraries')

curdir = os.path.dirname(os.path.realpath(__file__))
ext_modules = [ Extension(
    name="adaptive_anderson_solver",
    sources=["adaptive_anderson_solver.pyx"],
    libraries=["adaptive_anderson_solver"] + config['libraries'],
    include_dirs=[np.get_include(), os.path.dirname(__file__)],
    library_dirs=["./lib/"] + config['library_dirs'],
    extra_compile_args=["-O3"],
    extra_link_args=["-Wl,-rpath,./lib:"+curdir]
)]

setup(
    ext_modules=cythonize(ext_modules)
)
