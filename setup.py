#!/usr/bin/env python
"""
Institute for the Design of Advanced Energy Systems
"""
import sys
from setuptools import setup, find_packages


def warn(s):
    sys.stderr.write('*** WARNING *** {}\n'.format(s))


kwargs = dict(
    name='idaes',
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
    url="https://github.com/IDAES/idaes-tmp",
    license='BSD 3-clause',
    description="IDAES core framework",
    long_description=__doc__,
    data_files=[],
    keywords=["idaes", "energy systems"],
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

try:
    setup(setup_requires=['setuptools_scm'], use_scm_version=True, **kwargs)
except (ImportError, LookupError):
    default_version = '1.0.0'
    warn('Cannot use .git version: package setuptools_scm not installed '
         'or .git directory not present.')
    print('Defaulting to version: {}'.format(default_version))
    setup(**kwargs)
