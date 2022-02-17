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


DEPENDENCIES_FOR_PRERELEASE_VERSION = [
    "pyomo @ https://github.com/IDAES/pyomo/archive/6.2.zip"
]


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
        "nbconvert",
        "nbformat",
        "numpy",
        "networkx",
        "omlt==0.3", # fix the version for now as package evolves
        "pandas",
        "pint",
        "psutil",
        "pyomo>=6.2",
        "pytest",
        "pyyaml",
        "requests",  # for ui/fsvis
        "python-slugify", # for ui/fsvis
        "scipy",
        "sympy",
        "tinydb",
        "rbfopt",
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
        "prerelease": DEPENDENCIES_FOR_PRERELEASE_VERSION,
        "optional": [
            "tensorflow",  # idaes.surrogate.keras_surrogate
            # A Lee 11-Jan-22: no precompiled version of CoolProp available for Pyhton 3.9
            "coolprop; python_version < '3.9'"  # idaes.generic_models.properties.general.coolprop
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
            "*.json.gz",
            "*.dat",
            "*.h5",
            "*.trc",
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
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)

setup(**kwargs)
