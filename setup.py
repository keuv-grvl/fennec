import codecs
import os
import re
from setuptools.extension import Extension
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()

def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

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
              ('N_RTHREADS', "32")  # str(os.cpu_count() - 1))  # use all CPU
            ]
          )

setup(name = 'fennec',
      version = find_version("fennec", "__init__.py"),
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
