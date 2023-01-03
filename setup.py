from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np
import os.path

import sys
import site
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

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
    sources=["src/adaptive_anderson_solver/adaptive_anderson_solver.pyx"],
    libraries=["adaptive_anderson_solver"] + config['libraries'],
    include_dirs=[
      np.get_include(),
      os.path.join(curdir, 'src', 'fortran'),
      os.path.join(curdir),
      ],
    library_dirs=[
      "./lib/",
      os.path.join(curdir, 'src', 'adaptive_anderson_solver'),
      ] + config['library_dirs'],
    extra_compile_args=["-O3"],
    extra_link_args=["-Wl,-rpath,$ORIGIN/adaptive_anderson_solver/:$ORIGIN/../lib:"+curdir]
)]

setup(
    name='adaptive_anderson_solver',
    ext_modules=cythonize(ext_modules),
    packages = ['adaptive_anderson_solver'],
    package_dir= { '':'src' },
    package_data={ '': ['libadaptive_anderson_solver.so']},
    exclude_package_data= { '': ['*.c']},
    license='BSD',
    author='Matyáš Novák',
    version='1.0.0',
    author_email = 'novakmat@fzu.cz',
    include_package_data=True
)
