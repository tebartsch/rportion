from os import path
from setuptools import setup, find_packages
from codecs import open

with open(
    path.join(path.abspath(path.dirname(__file__)), "README.md"), encoding="utf-8"
) as f:
    long_description = f.read()

setup(
    name="rportion",
    version="0.1.0",
    license="MIT",
    author="Tilmann Bartsch",
    url="https://github.com/tilmann-bartsch/rportion",
    description="Python data structure and operations for 2-dimensional rectilinear polygons",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Developers",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    keywords="rectangle polygon interval-tree",
    packages=find_packages(include=["rportion"]),
    python_requires="~= 3.7",
    install_requires=[
        "sortedcontainers ~= 2.2", "portion ~= 2.2.0"
    ],
    extras_require={
        "test": ["pytest ~= 7.0",
                 "coverage ~= 6.0",
                 "numpy ~= 1.21.6"],
        "docu": ["numpy ~= 1.21.6",
                 "matplotlib ~= 3.5.2",
                 "imageio ~= 2.19.5",
                 "pillow ~= 9.2.0",
                 "tqdm ~= 4.64.0"]
    },
)
