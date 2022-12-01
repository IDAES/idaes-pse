#!/usr/bin/env python
"""
Institute for the Design of Advanced Energy Systems
"""
from pathlib import Path
import os
import sys
from setuptools import setup, find_namespace_packages
from typing import List, Tuple


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
README = README[README.find("#") :]  # ignore everything before title


def rglob(path, glob):
    """Return list of paths from `path` matching `glob`."""
    p = Path(path)
    return list(map(str, p.rglob(glob)))


DEPENDENCIES_FOR_PRERELEASE_VERSION = []

# For included DMF data
DMF_DATA_ROOT = "data"


def dmf_data_files(root: str = DMF_DATA_ROOT) -> List[Tuple[str, List[str]]]:
    """Generate a list of pairs (directory, [files..]), covering all the DMF data
    files, for the `data_files` option to :func:`setup()`.
    """
    file_list = [
        (
            root,
            [f"{root}/config.yaml", f"{root}/resourcedb.json"],
        )
    ]
    files_root = Path(root) / "files"
    for files_subdir in files_root.glob("*"):
        file_names = [f.as_posix() for f in files_subdir.glob("*")]
        if file_names:  # empty for non-directories and empty directories
            file_list.append((files_subdir.as_posix(), file_names))
    return file_list


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
        "click>=8",
        "colorama",
        "distro",  # help identify linux distros for binary downloads
        "flask",  # for ui/fsvis
        "flask-cors",
        "jupyter",
        "lxml",
        "matplotlib",
        "nbconvert",
        "nbformat",
        "numpy",
        "networkx",
        "omlt==0.3.1",  # fix the version for now as package evolves
        "pandas<1.5",
        "pint",
        "psutil",
        "pyomo>=6.4.3",
        "pytest",
        "pyyaml",
        "requests",  # for ui/fsvis
        "python-slugify",  # for ui/fsvis
        "scipy",
        "sympy",
        "tinydb",
        "rbfopt",
        "xlrd",  # for DMF read of old .xls Excel files
        "openpyxl",  # for DMF read of new .xls Excel files
        # lbianchi-lbl: see https://github.com/IDAES/idaes-pse/issues/661
        "ipython<8.0.0",
    ],
    entry_points={
        "console_scripts": [
            "dmf = idaes.dmf.cli:base_command",
            "idaes = idaes.commands.base:command_base",
        ]
    },
    # Only installed if [<key>] is added to package name
    extras_require={
        "optional": [
            "tensorflow",  # idaes.core.surrogate.keras_surrogate
            "gridx-prescient>=2.1",  # idaes.tests.prescient
            # A Lee 11-Jan-22: no precompiled version of CoolProp available for Pyhton 3.9
            "coolprop; python_version < '3.9'",  # idaes.generic_models.properties.general.coolprop
        ],
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
            "*.dat",
            "*.h5",
            "*.pb",  # for Keras Surrogate folder
            "*.data-00000-of-00001",  # for Keras Surrogate folder
            "*.index",  # for Keras Surrogate folder
            "*.trc",
            "*.xlsx",  # idaes/dmf/tests/data_files - tabular import test files
        ]
    },
    include_package_data=True,
    data_files=dmf_data_files(),
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
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)

setup(**kwargs)
