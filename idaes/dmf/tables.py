###############################################################################
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
###############################################################################
"""
Table handling for DMF
"""
# stdlib
from typing import List, Tuple, Dict
import re

# ext
import pandas as pd

__author__ = "Dan Gunter"


class Table:
    """Represent a table stored in the DMF.
    """
    def __init__(self):
        self._data = pd.DataFrame({})
        self._units = {}
        self._filepath = None

    @property
    def data(self) -> pd.DataFrame:
        return self._data

    @property
    def units(self) -> Dict[str, str]:
        return self._units.copy()

    def read_csv(self, filepath, **kwargs) -> None:
        """Read the table from a CSV file using pandas' `read_csv()`.
        See `Pandas read_csv docs
        <https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html>`_
        for details.

        Existing table will be replaced.

        Args:
            filepath: Any valid first argument to pandas `read_csv`
            kwwargs: Keyword arguments passed to pandas `read_csv`

        Returns:
            None
        """
        self._data = pd.read_csv(filepath, **kwargs)
        self._extract_units()
        self._filepath = filepath

    def read_excel(self, filepath, **kwargs) -> None:
        """Read the table from a CSV file using pandas' `read_excel()`.
        See `Pandas read_excel docs
        <https://pandas.pydata.org/docs/reference/api/pandas.read_excel.html>`_
        for details.

        Existing table will be replaced.

        Args:
            filepath: Any valid first argument to pandas `read_excel`
            **kwargs: Keyword arguments passed to pandas `read_excel`

        Returns:
            None

        Raises:
            ValueError, if more than one Excel sheet is returned
        """
        # Workaround for older versions of Python/Pandas (python 3.6):
        # set engine explicitly to openpyxl for *.xlsx files
        v = [int(_) for _ in pd.__version__.split(".")]
        if v[0] <= 1 and v[1] <= 1:  # version < 1.2.0
            from io import BufferedIOBase, RawIOBase
            import os
            # if it's a file and has xlsx extension, set engine
            if not isinstance(filepath, (BufferedIOBase, RawIOBase)):
                ext = os.path.splitext(str(filepath))[-1]
                if ext == ".xlsx":
                    kwargs["engine"] = "openpyxl"

        data = pd.read_excel(filepath, **kwargs)
        if isinstance(data, dict):
            raise ValueError(f"Read from excel file must return a single sheet, "
                             f"but sheet_name='{kwargs.get('sheet_name', '?')}' "
                             f"returned {len(data)} sheets: {list(data.keys())}")
        self._data = data
        self._extract_units()
        self._filepath = filepath

    def _extract_units(self):
        new_names, units_dict = {}, {}

        for name in self._data.columns:
            base_name, units = self._split_units(name)
            new_names[name] = base_name
            units_dict[base_name] = units

        self._data.rename(columns=new_names, inplace=True)
        self._units = units_dict

    @staticmethod
    def _split_units(name) -> Tuple[str, str]:
        m = re.match(r"(.*)\s*(\(.*\)|\[.*])", name)
        if m is None:
            return name, ""
        else:
            new_name = m.group(1).strip()
            unit = m.group(2)[1:-1]  # strip parentheses
            return new_name, unit

    def add_to_resource(self, rsrc):
        rsrc.v["data"]["table"] = self.as_dict()

    @classmethod
    def get_from_resource(cls, rsrc) -> "Table":
        """Get an instance of this class from data in a Resource.

        Args:
            rsrc: A DMF resource instance

        Returns:

        """
        try:
            table_dict = rsrc.v["data"]["table"]
        except KeyError:
            raise KeyError("No table in resource")
        return cls.from_dict(table_dict)

    def as_dict(self) -> Dict:
        header = list(self._data.columns)
        d = {
            column: {
                "units": self._units[column],
                "values": list(self._data[column])
            }
            for column in header
        }
        return d

    @classmethod
    def from_dict(cls, data) -> "Table":
        tbl = Table()
        dataframe_dict = {}
        for column, info in data.items():
            dataframe_dict[column] = info["values"]
            tbl._units[column] = info["units"]
        tbl._data = pd.DataFrame(dataframe_dict)
        return tbl

