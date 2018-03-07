from setuptools import setup
from fennec import __version__

requirements = open('requirements.txt').read().splitlines()

setup(name='fennec',
      version=__version__,
      description='DNA sequence modeling for machine learning',
      url='https://github.com/keuv-grvl/fennec',
      author='Kévin Gravouil',
      author_email='kevin.gravouil@uca.fr',
      license='MIT',
      packages=['fennec'],
      install_requires=requirements,
      #see: https://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research'
        ]
      )
