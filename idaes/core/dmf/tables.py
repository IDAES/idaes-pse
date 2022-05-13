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
Table handling for DMF.

The main class defined here is :class:`Table`. It provides constructor methods
for reading from Excel and CSV files. There is a convention defined for
indicating units in column headers so that this code can split the unit from
the column name. Other methods are defined for adding and extracting tables
from DMF :class:`idaes.core.dmf.resource.Resource` objects.

In the simplest case, you would create a new DMF resource for a CSV table like this::

    from idaes.core.dmf.resource import Resource
    resource = Resource()
    resource.add_table("my_file.csv")
    # you can now save this resource in the DMF

Then you could retrieve and use that table like this::

    # retrieve resource from the DMF
    table = resource.tables["my_file.csv"]
    dataframe = table.data    # Pandas dataframe
    units = table.units       # Units extracted from header row (strings)

See also, on the DMF Resource class:

    * :meth:`idaes.core.dmf.resource.Resource.add_table`
    * :attr:`idaes.core.dmf.resource.Resource.tables`

"""
# stdlib
from typing import List, Tuple, Dict
import re

# ext
import pandas as pd

# Local
from idaes.core.dmf.resource import Resource

__authors__ = ["Dan Gunter (LBNL)"]
__author__ = __authors__[0]


class DataFormatError(Exception):
    def __init__(self, source, problem):
        message = f"in {source}: {problem}"
        super().__init__(self, message)


class Table:
    """Represent a table stored in the DMF.

    Tables are expected to have a header row with optional units, which if present
    are encoded in [square brackets]. Whitespace is ignored between the column
    name and the units. For example::

            T [C], P [bar], G0/RT H2O, G0/RT NaCl [-], A phi [(kg/mol^0.5]
            0, 1, -23.4638, -13.836, 0.3767
    """

    def __init__(self):
        """Create new, empty, table.

        Use :meth:`read_csv` or :meth:`read_excel` to populate the table with data.
        """
        self._data = pd.DataFrame({})
        self._units = {}
        self._filepath = None
        self._desc = ""

    @property
    def data(self) -> pd.DataFrame:
        """Pandas dataframe for data."""
        return self._data

    @property
    def units_dict(self) -> Dict[str, str]:
        """Units as a dict keyed by table column name."""
        return self._units.copy()

    @property
    def units_list(self) -> List[str]:
        """Units in order of table columns."""
        return [self._units[c] for c in self._data.columns]

    #: Shorthand for getting list of units
    units = units_list

    @property
    def description(self):
        return self._desc

    @description.setter
    def description(self, value):
        self._desc = value

    @staticmethod
    def read_table(filepath, inline: bool, file_format: str) -> "Table":
        """Determine the input file type, then construct a new Table object
        by calling one of :meth:`Table.read_csv` or :meth:`Table.read_excel`.

        Args:
            filepath: Any valid first argument to pandas `read_csv`
            inline: If True, read the whole table in; otherwise just get the
                    column names and units from the header row.
            file_format: One of 'infer', 'csv', or 'excel'. For 'infer',
                         use the file extension (and only the extension) to
                         determine if it's a CSV or Excel file.

        Returns:
            Constructed Table object

        Raises:
            IOError: If the input cannot be read or parsed
        """
        fmt = file_format.lower()
        name = filepath.name
        if fmt == "infer":
            if name.endswith(".csv"):
                fmt = "csv"
            elif name.endswith(".xls") or name.endswith(".xlsx"):
                fmt = "excel"
            else:
                raise ValueError(f"Cannot infer file format for '{name}'")
        elif fmt not in ("csv", "excel"):
            raise ValueError(f"Unknown file format '{fmt}'; must be csv or excel")

        # create a new table to work with
        table = Table()

        # set up keywords to read only header row if we are not including data inline
        kwargs = {}
        if not inline:
            kwargs["nrows"] = 0

        # read the table (or at least its header)
        try:
            if fmt == "csv":
                table.read_csv(filepath, **kwargs)
            elif fmt == "excel":
                table.read_excel(filepath, **kwargs)
        except Exception as err:
            raise IOError(f"Cannot read '{filepath}': {err}")

        return table

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
            raise ValueError(
                f"Read from excel file must return a single sheet, "
                f"but sheet_name='{kwargs.get('sheet_name', '?')}' "
                f"returned {len(data)} sheets: {list(data.keys())}"
            )
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

    #: Regular expression for extracting units from column names.
    #: In plain English, the following forms are expected for a
    #: column name: "Name", "Name[Units]", "Longer Name With $% Chars [ Units ]"
    #: For both the Name and the Units, any sequence of characters valid
    #: in the current encoding are acceptable (except, of course, a "["
    #: in the name, which means start-of-units)
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
            raise DataFormatError(
                name,
                "No recognized column name. Expected syntax is "
                "'name' or 'name [units]'",
            )
        new_name = m.group("name").strip()
        unit = m.group("units")
        if unit == "-" or unit is None:
            unit = ""  # normalize empty units to empty string
        else:
            unit = unit.strip()  # note: may also end up
        return new_name, unit

    def add_to_resource(self, rsrc: Resource):
        """Add the current table, inline, to the given resource.

        Args:
            rsrc: A DMF :class:`Resource` instance

        Returns:
            None
        """
        rsrc.data[Resource.TABLE_FIELD] = self.as_dict()

    @classmethod
    def from_resource(cls, rsrc: Resource) -> Dict[str, "Table"]:
        """Get an instance of this class from data in the given resource.

        Args:
            rsrc: A DMF :class:`Resource` instance

        Returns:
            Dictionary of tables in resource. If there is only one inline
            table, the dictionary is of length one with only key "" (empty string).
            If there are multiple tables referenced by file the dictionary
            keys are the (relative) file names.
            If there are no tables in this resource, raises KeyError.

        Raises:
            KeyError: if there are no tables in this resource
        """
        data = rsrc.v["data"]
        if Resource.TABLE_FIELD in data:
            # Single inline resource
            table_ = cls.from_dict(data[Resource.TABLE_FIELD])
            return {"": table_}
        elif Resource.TABLE_INFO_FIELD in data:
            # One or more files
            tables = {}
            for idx, path in enumerate(rsrc.get_datafiles()):
                table_ = cls.read_table(path, True, "infer")
                table_.description = rsrc.v[Resource.DATAFILES_FIELD][idx].get(
                    "desc", ""
                )
                tables[path.name] = table_
            return tables
        else:
            raise KeyError("No table in resource")

    def as_dict(self, values=True) -> Dict:
        """Get the representation of this table as a dict.

        Args:
            values: If True, include the values in the dict. Otherwise only
                    include the units for each column.

        Returns:
            Dictionary with the structure accepted by :meth:`from_dict`.
            If the "values" argument is False, that key will be missing from
            the dict for each column.
        """
        header = list(self._data.columns)
        d = {}
        for column in header:
            d[column] = {"units": self._units[column]}
            if values:
                d[column]["values"] = list(self._data[column])
        return d

    @classmethod
    def from_dict(cls, data: Dict) -> "Table":  # unquote in Py3.7+ see PEP563
        """Create a new Table object from a dictionary of data and units.

        Args:
            data: Dictionary with the following structure::

                {
                    'column-name-1': {
                        'units': 'unit',
                        'values': [ value, value, .. ]
                    },
                    'column-name-2': {
                        'units': 'unit',
                        'values': [ value, value, .. ]
                    },
                    ...etc...
                }

        Returns:
            :class:`Table` object
        """
        tbl = Table()
        dataframe_dict = {}
        for column, info in data.items():
            dataframe_dict[column] = info.get("values", [])
            tbl._units[column] = info.get("units", "")
        tbl._data = pd.DataFrame(dataframe_dict)
        return tbl
