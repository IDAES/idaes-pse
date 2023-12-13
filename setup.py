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


# For included DMF data
DMF_DATA_ROOT = "dmf_data"


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


class ExtraDependencies:
    """
    A convenience shorthand to define and combine dependencies for ``extras_require``.

    >>> extras = ExtraDependencies()
    >>> extras["ui"]
    ['requests', 'pint']
    >>> list(extras)
    ['ui', 'dmf', 'omlt', 'grid', 'coolprop', 'testing']
    >>> dict(extras)
    {'ui': ['requests', 'pint'], 'dmf': ['jsonschema', 'setuptools', 'traitlets', ...], ...}
    """

    ui = [
        "idaes-ui >= 0.23.12",
    ]
    _ipython = [
        'ipython <= 8.12; python_version == "3.8"',
    ]
    dmf = [
        # all modules relative to idaes.core.dmf
        "jsonschema",  # commands, resource, workspace
        "setuptools",  # provides pkg_resources?
        "traitlets",  # dmfbase
        "lxml",  # help
        "seaborn",  # model_data (optional^2)
        "PyPDF2",  # model_data (optional^2)
        "colorama",  # util
        *_ipython,  # magics
        "pyyaml",  # workspace
        "tinydb",  # resourcedb
        "xlrd",  # tables (implicitly by pandas.read_excel())
        "openpyxl",  # tables (implicitly by pandas.read_excel())
    ]
    omlt = [
        "omlt==1.1",  # fix the version for now as package evolves
        "tensorflow",
    ]
    grid = [
        "gridx-prescient>=2.2.1",  # idaes.tests.prescient
    ]
    coolprop = [
        "coolprop",  # idaes.generic_models.properties.general.coolprop
    ]
    testing = [
        "pytest",
        "addheader",
        "pyyaml",
    ]

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
            "dmf_data*",
        ]
    ),
    # Put abstract (non-versioned) deps here.
    # Concrete dependencies go in requirements[-dev].txt
    install_requires=[
        "pyomo >= 6.7.0",
        "pint",  # required to use Pyomo units
        "networkx",  # required to use Pyomo network
        "numpy",
        # pandas constraint added on 2023-08-30 b/c bug in v2.1
        # see IDAES/idaes-pse#1253
        "pandas!=2.1.0",
        "scipy",
        "sympy",  # idaes.core.util.expr_doc
        "matplotlib",
        "click>=8",
    ],
    entry_points={
        "console_scripts": [
            "dmf = idaes.dmf.cli:base_command",
            "idaes = idaes.commands.base:command_base",
        ]
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
            "*.xlsx",  # idaes/dmf/tests/data_files - tabular import test files
            "*.nl",
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
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)

setup(**kwargs)
