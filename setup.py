from setuptools import setup
import fennec

setup(name='fennec',
      version=fennec.__version__,
      description='DNA sequence modeling for machine learning',
      url='https://github.com/keuv-grvl/fennec',
      author='Kévin Gravouil',
      author_email='kevin.gravouil@uca.fr',
      license='MIT',
      packages=['fennec'],
      install_requires=[
          'scikit-bio',
          'scikit-learn',
      ],
      zip_safe=False)
