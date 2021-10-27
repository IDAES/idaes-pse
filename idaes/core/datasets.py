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
from idaes.dmf.datasets import Publication, AvailableResult

__authors__ = ["Dan Gunter (LBNL)"]
__author__ = __authors__[0]


class Pitzer(Publication):
    """Pitzer(1984) publication and related tables."""

    def __init__(self, workspace=None):
        super().__init__("Pitzer:1984", workspace=workspace)


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
    for obj in globals().copy():
        if isinstance(obj, type) and issubclass(obj, Publication):
            result[obj.__name__] = AvailableResult(obj, doc_desc(obj))
    return result
