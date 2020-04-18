#!/usr/bin/env python

from setuptools import find_packages, setup

setup(
    name="pymol-labimm",
    packages=find_packages(),
    install_requires=["pandas", "scipy", "requests", "cached_property"],
)
