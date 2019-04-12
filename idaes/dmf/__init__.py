"""
IDAES Data Management Framework (DMF)

The DMF lets you save, search, and retrieve provenance related
to your models.
"""
__author__ = 'Dan Gunter'

import logging

from .dmfbase import DMF, DMFConfig           # noqa: F401
from .userapi import get_workspace            # noqa: F401
from .userapi import find_property_packages   # noqa: F401
from .userapi import index_property_packages  # noqa: F401
from . import resource                        # noqa: F401
# DMF version is the same as IDAES version
from idaes import __version__                 # noqa


# default log format
h = logging.StreamHandler()
h.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] '
                                 '%(name)s: %(message)s'))
logging.getLogger('idaes.dmf').addHandler(h)
logging.getLogger('idaes.dmf').propagate = False
