#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name="ler",
    version="0.1.6",
    description="Lensing Rates",
    author="Hemantakumar",
    license="MIT",
    author_email="hemantaphurailatpam@gmail.com",
    url="https://github.com/hemantaph/ler",
    packages=find_packages(),
    install_requires=[
        "setuptools>=67.8.0",
        "numpy>=1.18",
        "numba==0.56.4",
        "bilby>=1.0.2",
        "pycbc<=2.0.6",
        "gwsnr>=0.1",
        "scipy>=1.9.0",
        "lenstronomy>=1.10.4",
        "astropy>=5.1",
        "tqdm>=4.64.1",
        "pointpats>=2.3",
        "shapely>=2.0.1",
    ],
)
