"""
API for accessing core IDAES datasets

Usage, e.g., for the Pitzer(1984) data::

    from idaes.core.datasets import Pitzer
    pitzer = Pitzer()
    gibbs_data = pitzer.get_table("Standard G").data

"""
from typing import Union, List
from idaes.dmf import datasets
from idaes.dmf.tables import Table

__authors__ = ["Dan Gunter (LBNL)"]
__author__ = __authors__[0]


class _Publication:
    """Abstract superclass for all publication-derived datasets.

    Do not instantiate directly.
    """

    def __init__(self, name, workspace=None):
        self._ds = datasets.PublicationDataset(workspace=workspace)
        self._pub, self._tables = self._ds.retrieve(name)

    def list_tables(self) -> List[str]:
        return list(self._tables.keys())

    def get_table(self, name) -> Union[Table, None]:
        return self._tables.get(name, None)


class Pitzer(_Publication):
    """Pitzer(1984) publication and related tables.
    """
    def __init__(self, workspace=None):
        super().__init__("Pitzer:1984", workspace=workspace)
