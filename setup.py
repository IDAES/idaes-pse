#!/usr/bin/env python
"""
Institute for the Design of Advanced Energy Systems
"""
from pathlib import Path
import os
import sys
from setuptools import setup, find_namespace_packages


def warn(s):
    sys.stderr.write("*** WARNING *** {}\n".format(s))


def get_version():
    code_file = os.path.join("idaes", "ver.py")
    code = open(code_file).read()
    local_namespace = {}
    exec(code, {}, local_namespace)
    return local_namespace["__version__"]


NAME = "idaes-pse"
VERSION = get_version()
README = open("README.md").read()
README = README[README.find("#"):]  # ignore everything before title


def rglob(path, glob):
    """Return list of paths from `path` matching `glob`.
    """
    p = Path(path)
    return list(map(str, p.rglob(glob)))


kwargs = dict(
    zip_safe=False,
    name=NAME,
    version=VERSION,
    packages=find_namespace_packages(),
    # Put abstract (non-versioned) deps here.
    # Concrete dependencies go in requirements[-dev].txt
    install_requires=[
        # idaes core / dmf
        "backports.shutil_get_terminal_size",
        "bunch",
        "click",
        "colorama",
        "flask",  # for ui/fsvis
        "flask-cors",
        "jupyter",
        "lxml",
        "matplotlib",
        "mock",
        "nbconvert",
        "nbformat",
        "numpy",
        "networkx",
        "pandas",
        "pint",
        "psutil",
        "pyutilib>=6.0.0",
        "pyomo>5.7.1,<6",
        "pytest",
        "pyyaml",
        "requests",  # for ui/fsvis
        "python-slugify", # for ui/fsvis
        "sympy",
        "tinydb",
        "rbfopt",
    ],
    entry_points={
        "console_scripts": [
            "dmf = idaes.dmf.cli:base_command",
            "idaes = idaes.commands.base:command_base",
        ]
    },
    extras_require={
        # For developers. Only installed if [dev] is added to package name
        "dev": [
            "alabaster>=0.7.7",
            # temporarily hold coverage version due to avoid bug in coveralls
            # -alee 12/20/2019
            "coverage==4.5.4",
            "flake8",
            "jsonschema",
            "jupyter_contrib_nbextensions",
            "mock",
            "pylint",
            "pytest-cov",
            "python-coveralls",
            "snowballstemmer==1.2.1",
            # temporarily hold sphinx version to avoid bug with 3.x
            # -dang 4/22/2020
            "sphinx>=2.4.4,<3.0.0",
            # note: 4/22/2020, removed the version requirement here
            "sphinx-rtd-theme",
            "sphinxcontrib-napoleon>=0.5.0",
            "sphinx-argparse",
        ]
    },
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
        ]
    },
    include_package_data=True,
    data_files=[],
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
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)

setup(**kwargs)
