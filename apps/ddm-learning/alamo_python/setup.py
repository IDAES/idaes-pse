#!/usr/bin/env python
"""
ALAMO Python Interface
"""
import sys
from setuptools import setup, find_packages


def warn(s):
    sys.stderr.write('*** WARNING *** {}\n'.format(s))


kwargs = dict(
    name='alamopy',
    packages=find_packages(),
    install_requires=[],
    extras_require={},
    package_data={},
    scripts=[],
    author='ALAMO team',
    author_email='idaes-dev@idaes.org',
    maintainer='TBD',
    url="https://github.com/IDAES/Thermodynamics-and-kinetics/",
    license='BSD 3-clause',
    description="ALAMO and RIPE wrappers for IDAES project",
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
