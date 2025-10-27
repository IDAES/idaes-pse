#!/usr/bin/env python
"""
Institute for the Design of Advanced Energy Systems
"""
from pathlib import Path
import os
import sys
from setuptools import setup, find_namespace_packages
from typing import List, Tuple


def get_version():
    code_file = os.path.join("idaes", "ver.py")
    code = open(code_file).read()
    local_namespace = {}
    exec(code, {}, local_namespace)
    return local_namespace["__version__"]


NAME = "idaes-pse"
VERSION = get_version()
README = open("README.md").read()
README = README[README.find("#") :]  # ignore everything before title


class ExtraDependencies:
    """
    A convenience shorthand to define and combine dependencies for ``extras_require``.

    >>> extras = ExtraDependencies()
    >>> extras["testing"]
    ['pytest', 'addheader', 'pyyaml']
    >>> list(extras)
    ['ui', 'omlt', 'grid', 'coolprop', 'testing']
    >>> dict(extras)
    {'testing': ['pytest', 'addheader', 'pyyaml'], 'omlt': ['omlt==1.1', 'tensorflow', ...], ...}
    """

    ui = [
        "idaes-ui",
    ]
    omlt = [
        "omlt==1.1",  # fix the version for now as package evolves
        "tensorflow",
        "onnx",
    ]
    grid = [
        "gridx-prescient>=2.2.1",  # idaes.tests.prescient
    ]
    coolprop = [
        "coolprop>=7.0",  # idaes.generic_models.properties.general.coolprop
    ]
    testing = ["pytest", "addheader", "pyyaml"]

    def __init__(self):
        self._data = dict(type(self).__dict__)

    def keys(self):
        for name, attr in self._data.items():
            if not name.startswith("_") and isinstance(attr, list):
                yield name

    def __getitem__(self, key):
        return self._data[key]


kwargs = dict(
    zip_safe=False,
    name=NAME,
    version=VERSION,
    packages=find_namespace_packages(
        include=[
            "idaes*",
        ]
    ),
    # Put abstract (non-versioned) deps here.
    # Concrete dependencies go in requirements[-dev].txt
    install_requires=[
        "pyomo == 6.9.4",  # Temporary pin to avoid Pylint issues TODO revert
        "pint >= 0.24.1",  # required to use Pyomo units. Pint 0.24.1 needed for Python 3.9 support
        "networkx",  # required to use Pyomo network
        "numpy>=1,<3",
        # pandas constraint added on 2023-08-30 b/c bug in v2.1
        # see IDAES/idaes-pse#1253
        "pandas!=2.1.0,<3",
        "scipy",
        "sympy",  # idaes.core.util.expr_doc
        "matplotlib",
        "click>=8",
    ],
    entry_points={
        "console_scripts": [
            "idaes = idaes.commands.base:command_base",
        ],
        "idaes.flowsheets": [
            "0D Fixed Bed TSA = idaes.models_extra.temperature_swing_adsorption.fixed_bed_tsa0d_ui",
        ],
    },
    # Only installed if [<key>] is added to package name
    extras_require=dict(ExtraDependencies()),
    package_data={
        # If any package contains these files, include them:
        "": [
            "*.template",
            "*.json",
            "*.yaml",
            "*.svg",
            "*.png",
            "*.jpg",
            "*.csv",
            "*.ipynb",
            "*.txt",
            "*.js",
            "*.css",
            "*.html",
            "*.json.gz",
            "*.dat",
            "*.h5",
            "*.pb",  # for Keras Surrogate folder
            "*.data-00000-of-00001",  # for Keras Surrogate folder
            "*.index",  # for Keras Surrogate folder
            "*.trc",
            "*.nl",
            "*.keras",  # idaes/core/surrogate/tests/data/keras_models
            "*.onnx",
        ]
    },
    include_package_data=True,
    maintainer="Keith Beattie",
    maintainer_email="ksbeattie@lbl.gov",
    url="https://idaes.org",
    license="BSD ",
    platforms=["any"],
    description="IDAES Process Systems Engineering Framework",
    long_description=README,
    long_description_content_type="text/markdown",
    keywords=[NAME, "energy systems", "chemical engineering", "process modeling"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)

setup(**kwargs)
