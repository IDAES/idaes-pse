#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
API for accessing core IDAES datasets

Usage, e.g., for the Pitzer(1984) data::

    from idaes.core.datasets import Pitzer
    pitzer = Pitzer()
    gibbs_data = pitzer.get_table("Standard G").data

"""
import re
from typing import Dict

# package
from idaes.core.dmf.datasets import Publication, AvailableResult

__authors__ = ["Dan Gunter (LBNL)"]
__author__ = __authors__[0]


class Pitzer(Publication):
    """Pitzer(1984) publication and related tables."""

    def __init__(self, workspace=None):
        super().__init__("Pitzer:1984", workspace=workspace)


def __f():
    pass


def available() -> Dict[str, AvailableResult]:
    """Find and return all available datasets.

    Returns:
        Mapping of dataset class name to its class and description
    """

    def doc_desc(cls):
        """Get description from docstring."""
        docstr = cls.__doc__
        # find first blank line or line ending in period (or whole string)
        m = re.search(r"[.]$|^$", docstr, flags=re.M)
        pos = m.start() if m else len(docstr)
        # return text as a single line
        return " ".join([x.strip() for x in docstr[:pos].strip().split("\n")])

    result = {}
    for key, obj in __f.__globals__.items():
        if (
            isinstance(obj, type)
            and issubclass(obj, Publication)
            and not obj.__name__ == "Publication"
        ):
            result[key] = AvailableResult(obj, doc_desc(obj))
    return result
