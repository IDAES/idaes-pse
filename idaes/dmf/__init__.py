"""
IDAES Data Management Framework (DMF)

The DMF lets you save, search, and retrieve provenance related
to your models.

This package is documented with Sphinx. To build the documentation,
change to the 'docs' directory and run, e.g., 'make html'.

.. automodule:: idaes.dmf.resource

"""
import logging

from .dmfbase import DMF                      # noqa: F401
from .userapi import get_workspace            # noqa: F401
from .userapi import find_property_packages   # noqa: F401
from .userapi import index_property_packages  # noqa: F401

# Version info.
# Uses `setuptools_scm <https://pypi.python.org/pypi/setuptools_scm>`_
from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = '0.0.0'

# default log format
h = logging.StreamHandler()
h.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] '
                                 '%(name)s: %(message)s'))
logging.getLogger('idaes.dmf').addHandler(h)
logging.getLogger('idaes.dmf').propagate = False
