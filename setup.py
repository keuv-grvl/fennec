import os
from distutils.core import Extension, setup

from fennec import __version__

# Add $CONDA_PREFIX/{lib,include} if we are in a Conda environment
try:
  libdir = os.environ['CONDA_PREFIX']
except:
  libdir = '.'

# List requirement from `requirements.txt`
requirements = open('requirements.txt').read().splitlines()

# Use the C implementation of VBGMM from https://github.com/BinPro/CONCOCT/
module1 = Extension(
            name = 'vbgmm',
            library_dirs = [
              libdir+'/lib',
              libdir+'/include'
            ],
            runtime_library_dirs = [
              libdir+'/lib',
              libdir+'/include'
            ],
            libraries = ['pthread', 'gsl', 'gslcblas'],
            include_dirs = ['c-concoct'],
            sources = ['c-concoct/vbgmmmodule.c'],
            define_macros = [
              ('N_RTHREADS', os.cpu_count())  # use all CPU
            ]
          )

setup(name = 'fennec',
      version = __version__,
      description = 'DNA sequence modeling for machine learning',
      url = 'https://github.com/keuv-grvl/fennec',
      author = 'Kévin Gravouil',
      author_email = 'kevin.gravouil@uca.fr',
      license = 'MIT',
      packages = ['fennec'],
      ext_modules = [module1],
      install_requires = requirements,
      #see: https://pypi.python.org/pypi?%3Aaction = list_classifiers
      classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research'
      ],
      include_package_data = True,
    )
