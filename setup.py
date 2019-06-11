#!/usr/bin/env python
"""
Institute for the Design of Advanced Energy Systems
"""
import os
import sys
from setuptools import setup, find_packages


def warn(s):
    sys.stderr.write("*** WARNING *** {}\n".format(s))


def get_version(file, name="__version__"):
    """Get the version of the package from the given file by
    executing it and extracting the given `name`.
    """
    path = os.path.realpath(file)
    local_namespace = {}
    exec(open(path).read(), {}, local_namespace)
    return local_namespace[name]


NAME = "idaes"
VERSION = get_version(os.path.join(NAME, "ver.py"))

kwargs = dict(
    name=NAME,
    version=VERSION,
    packages=find_packages(),
    # Put abstract (non-versioned) deps here.
    # Concrete dependencies go in requirements[-dev].txt
    install_requires=[
        "backports.shutil_get_terminal_size",
        "blessings",
        "bokeh",
        "bunch",
        "click",
        "humanize",
        "jupyter",
        "lxml",
        "matplotlib",
        "pandas",
        "pendulum",
        "psutil",
        "pytest",
        "pyyaml",
        "sympy",
        "tinydb",
        "toml",
    ],
    entry_points="""
    [console_scripts]
    dmf=idaes.dmf.cli:base_command
    """,
    extras_require={
        # For developers. Only installed if [dev] is added to package name
        "dev": [
            "alabaster>=0.7.7",
            "coverage",
            "flake8",
            "flask>=1.0",
            "flask-bower",
            "flask-restful",
            "jsonschema",
            "jupyter_contrib_nbextensions",
            "mock",
            "networkx",
            "pytest-cov",
            "python-coveralls",
            "six",
            "snowballstemmer==1.2.1",
            "sphinx-rtd-theme>=0.1.9",
            "sphinxcontrib-napoleon>=0.5.0",
            "sphinx-argparse",
        ]
    },
    package_data={
        # If any package contains *.template, *.json files, *.dll files, or
        # *.so file, include them:
        "": ["*.template", "*.json", "*.dll", "*.so", "*.svg"]
    },
    author="IDAES Team",
    author_email="idaes-dev@idaes.org",
    maintainer="Keith Beattie",
    url="https://github.com/IDAES/idaes",
    license="BSD 3-clause",
    description="IDAES core framework",
    long_description=__doc__,
    data_files=[],
    keywords=[NAME, "energy systems"],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)

setup(**kwargs)
