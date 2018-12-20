#!/usr/bin/env python
"""
Institute for the Design of Advanced Energy Systems
"""
import io
import os
import sys
from setuptools import setup, find_packages


def warn(s):
    sys.stderr.write('*** WARNING *** {}\n'.format(s))


def get_version(file, name='__version__'):
    """Get the version of the package from the given file by
    executing it and extracting the given `name`.
    """
    path = os.path.realpath(file)
    version_ns = {}
    with io.open(path, encoding="utf8") as f:
        exec(f.read(), {}, version_ns)
    return version_ns[name]


NAME = 'idaes'
VERSION = get_version(os.path.join(NAME, 'ver.py'))

kwargs = dict(
    name=NAME,
    version=VERSION,
    packages=find_packages(),
    install_requires=[],
    extras_require={},
    package_data={
        # If any package contains *.template, *.json files, *.dll files, or
        # *.so file, include them:
        '': ['*.template', '*.json', '*.dll', '*.so']
    },
    author='IDAES Team',
    author_email='idaes-dev@idaes.org',
    maintainer='Keith Beattie',
    url='https://github.com/IDAES/idaes',
    license='BSD 3-clause',
    description="IDAES core framework",
    long_description=__doc__,
    data_files=[],
    keywords=[NAME, "energy systems"],
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules"
    ],
)

setup(**kwargs)
