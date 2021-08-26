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
        "click<=7.1.2", # problems with 8.x
        "colorama",
        "flask",  # for ui/fsvis
        "flask-cors",
        "jupyter",
        # pinning pywin32 to version 225 as a workaround for Python 3.8 compatibility issues
        # (ImportError: DLL load failed while importing ...)
        # for more information see e.g. https://stackoverflow.com/a/62249872
        "pywin32==225; sys_platform=='win32' and python_version>='3.8'",
        "lxml",
        "matplotlib",
        "nbconvert",
        "nbformat",
        "numpy",
        "networkx",
        "pandas",
        "pint",
        "psutil",
        "pyomo>=6.1.2",
        "pytest",
        "pyyaml",
        "requests",  # for ui/fsvis
        "python-slugify", # for ui/fsvis
        "scipy",
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
    # Only installed if [<key>] is added to package name
    extras_require={
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
            "*.json.gz",
            "*.dat"
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
