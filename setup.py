#!/usr/bin/env python

import pathlib

from setuptools import find_packages, setup

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()


setup(
    name="pymol-labimm",
    version="0.8.1",
    description="Some PyMOL utilities",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/pslacerda/pymol-labimm",
    author="Pedro Sousa Lacerda",
    author_email="pslacerda@gmail.com",
    license="MIT",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: MIT License",
        "Environment :: Plugins",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(),
    # TODO remove unused requirements
    install_requires=[
        "lxml",
        "pandas",
        "scipy",
        "requests",
        "cached_property",
        "patool",
        "matplotlib",
        "seaborn",
    ],
)
