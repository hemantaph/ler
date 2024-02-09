#!/usr/bin/env python
from setuptools import setup, find_packages
import sys
# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# check that python version is 3.7 or above
python_version = sys.version_info
if python_version < (3, 10):
    sys.exit("Python < 3.10 is not supported, aborting setup")

setup(
    name="ler",
    version="0.3.3",
    description="Gravitational waves Lensing Rates",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Hemantakumar",
    license="MIT",
    author_email="hemantaphurailatpam@gmail.com",
    url="https://github.com/hemantaph/ler",
    packages=find_packages(),
    python_requires='>=3.10',
    install_requires=[
        "setuptools>=65.5.0",
        "matplotlib>=3.4.2",
        "pycbc>=1.18.0",
        "numpy>=1.18",
        "numba>=0.57.1",
        "bilby>=1.0.2",
        "gwsnr>=0.2.0",
        "scipy>=1.9.0",
        "lenstronomy>=1.10.4",
        "astropy>=5.1",
        "tqdm>=4.64.1",
        "pointpats>=2.3",
        "shapely>=2.0.1",
        "gwcosmo>=2.0.0",
    ],
)
