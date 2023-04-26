#!/usr/bin/env python
from setuptools import setup, find_packages
setup(name='ler',
      version='0.1.3',
      description='Lensing Rates',
      author='Hemantakumar',
      license="MIT",
      author_email='hemantaphurailatpam@gmail.com',
      url='https://github.com/hemantaph/ler',
      packages=find_packages(),
      install_requires=[
        "setuptools>=61.1.0",
        "numpy>=1.18",
        "bilby>=1.0.2",
        "pycbc>=2.0.4",
        "quintet>=0.1",
        "scipy>=1.9.0",
        "lenstronomy>=1.10.4",
        "astropy>=5.1",
        "gwcosmo>=1.0.0",
      ]
     )
