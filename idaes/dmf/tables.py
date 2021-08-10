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
Table handling for DMF.


To use the API, create an instance of the :class:`Table` class.
Then (or previously) use the DMF's :class:`idaes.dmf.resource.Resource` class
to create a resource object, and call :meth:`Table.add_to_resource()` to add
the table to the resource.
"""
# stdlib
from typing import List, Tuple, Dict
import re

# ext
import pandas as pd

__author__ = "Dan Gunter"


class DataFormatError(Exception):
    def __init__(self, source, problem):
        message = f"in {source}: {problem}"
        super().__init__(self, message)


class Table:
    """Represent a table stored in the DMF.

    Tables are expected to have a header row with optional units,
    encoded in [square brackets]. Whitespace is ignored between the column
    name and the units. For example::

            T [C], P [bar], G0/RT H2O [-], G0/RT NaCl [-], A phi [(kg/mol^0.5]
            0, 1, -23.4638, -13.836, 0.3767
    """
    def __init__(self):
        """Create new, empty, table.

        Use :meth:`read_csv` or :meth:`read_excel` to populate the table with data.
        """
        self._data = pd.DataFrame({})
        self._units = {}
        self._filepath = None

    @property
    def data(self) -> pd.DataFrame:
        """Pandas dataframe for data.
        """
        return self._data

    @property
    def units_dict(self) -> Dict[str, str]:
        """Units as a dict keyed by table column name.
        """
        return self._units.copy()

    @property
    def units_list(self) -> List[str]:
        """Units in order of table columns.
        """
        return [self._units[c] for c in self._data.columns]

    #: Shorthand for getting list of units
    units = units_list

    def read_csv(self, filepath, **kwargs) -> None:
        """Read the table from a CSV file using pandas' `read_csv()`.
        See `Pandas read_csv docs
        <https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html>`_
        for details.

        Existing table will be replaced.

        Args:
            filepath: Any valid first argument to pandas `read_csv`
            kwargs: Keyword arguments passed to pandas `read_csv`

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
            ValueError: if more than one Excel sheet is returned
            DataFormatError: if the input data or header is invalid
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

    UNITS_REGEX = r"""
        (?P<name>[^[]+) # column name
        (?:\s*\[        # start of [units] section
        (?P<units>.*?)  # column units
        \])?            # end of [units] section, which is optional
        """

    @classmethod
    def _split_units(cls, name) -> Tuple[str, str]:
        m = re.match(cls.UNITS_REGEX, name, flags=re.X)
        if m is None:
            raise DataFormatError(name, "No recognized column name. Expected format "
                                        "is 'name [units]', where [units] is optional")
        new_name = m.group("name").strip()
        unit = m.group("units")
        if unit == "-" or unit is None:
            unit = ""  # normalize empty units to empty string
        else:
            unit = unit.strip()
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

