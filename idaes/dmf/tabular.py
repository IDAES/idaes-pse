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
Tabular data handling
"""
# standard
import abc
import csv
import json
import logging
import re
# third-party
import jsonschema
try:
    import pandas as pd
    import numpy as np
except ImportError:
    np, pd = None, None
# local
from . import errors
from .util import get_file

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)


class Fields(object):
    """Constants for field names.
    """
    DATA, META = 'data', 'meta'
    DTYPE, AUTH, INFO, TITLE, DATE = \
        'datatype', 'authors', 'info', 'title', 'date'
    VALS, ROWS = 'values', 'rows'
    #: Keys for data mapping
    DATA_NAME = 'name'
    DATA_UNITS = 'units'
    DATA_VALUES = 'values'
    DATA_ERRORS = 'errors'
    DATA_ERRTYPE = 'error_type'
    # Used during parsing
    COLTYPE = 'type'

# --------------------------------------------------------------------------
# Schemas


DATA_SCHEMA_DEF = {
    "type": "object",
    "properties": {
        "name": {
            "type": "string",
            "examples": [
                "Density",
                "r"
            ]
        },
        "units": {
            "type": "string",
            "examples": [
                "mPa-s",
                "K"
            ]
        },
        "values": {
            "description": "Column of numeric values",
            "type": "array",
            "items": {
                "type": "number"
            },
            "examples": ["[2.6, 6.21]"]
        },
        "error_type": {
            "description": "Type for error values",
            "type": "string"
        },
        "errors": {
            "description": "Column of numeric errors",
            "type": "array",
            "items": {
                "type": "number"
            },
            "examples": [
                "[0.001, 0.035]"
            ]
        },
        "type": {
            "description": "Type of column",
            "enum": ["state", "property"]
        }
    },
    "required": [
        "name",
        "units",
        "values"
    ],
    "additionalProperties": False
}

METADATA_SCHEMA_DEF = {
    "type": "object",
    "properties": {
        "datatype": {
            "description": "name of the data type",
            "type": "string",
            "examples": ["MEA"]
        },
        "info": {
            "description": "Additional information about the source "
                           "(i.e. publication)",
            "type": "string",
            "examples": [
                "J. Chem. Eng. Data, 2009, Vol 54, pg. 3096-30100"]
        },
        "notes": {
            "description": "Free-form text with notes about the data",
            "type": "string",
            "examples": ["r is MEA weight fraction in aqueous soln."]
        },
        "authors": {
            "description": "Author list in format Last1, First1, Last2,"
                           " First2, etc.",
            "type": "string",
            "examples": ["Amundsen, T.G., Lars, E.O., Eimer, D.A."]
        },
        "title": {
            "description": "Title of the source (e.g. publication"
                           " title)",
            "type": "string",
            "examples": [
                "Density and Viscosity of Monoethanolamine + .etc."]
        },
        "date": {
            "description": "Date of source data",
            "type": "string",
            "examples": ["2009"]
        }
    },
    "required": [
        "datatype",
        "authors",
        "title",
        "date"
    ],
    "additionalProperties": True
}

TABLE_SCHEMA = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "id": "http://idaes.org/table",
    "type": "object",
    "properties": {
        "meta": {
            "description": "List of information about the data source",
            "type": "array",
            "items": METADATA_SCHEMA_DEF
        },
        "data": {
            "description": "Measured data columns",
            "type": "array",
            "items": DATA_SCHEMA_DEF
        }
    },
    "required": ["meta", "data"],
    "additionalProperties": False
}

COLUMN_SCHEMA = {
    "$schema": "http://json-schema.org/draft-04/schema#",
    "id": "http://idaes.org/table",
    "type": "array",
    "items": DATA_SCHEMA_DEF,
}

# Schemas
# --------------------------------------------------------------------------


class TabularObject(object):
    """Abstract Property data class.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def as_dict(self):
        """Return Python dict representation.
        """
        pass


class Table(TabularObject):
    """Tabular data and metadata together (at last!)
    """
    _validator = jsonschema.Draft4Validator(TABLE_SCHEMA)

    def __init__(self, data=None, metadata=None):
        """Wrapper object for data + metadata of properties.

        Args:
            data (list|TabularData): Raw data dictionaries
            metadata (dict|list|Metadata): Metadata dictionar(ies)
        """
        if isinstance(data, list):
            self._data = TabularData(data)
        elif isinstance(data, TabularData):
            self._data = data
        else:
            raise TypeError('Expected list or TabularData object for: {}'
                            .format(data))
        self._meta = []
        if metadata:
            if isinstance(metadata, list) or isinstance(metadata, tuple):
                for m in metadata:
                    self.add_metadata(m)
            else:
                self.add_metadata(metadata)

    def __iter__(self):
        yield 'data', self._data.as_list()
        yield 'meta', [m.as_dict() for m in self._meta]

    def add_metadata(self, m):
        if isinstance(m, dict):
            obj = Metadata(values=m)
        elif isinstance(m, Metadata):
            obj = m
        else:
            raise TypeError('Expected dict or Metadata object for: {}'
                            .format(m))
        self._meta.append(obj)

    @property
    def data(self):
        return self._data

    @property
    def metadata(self):
        return self._meta

    def as_dict(self):
        """Represent as a Python dictionary.

        Returns:
            (dict) Dictionary representation
        """
        return {k: v for k, v in self}

    def dump(self, fp, **kwargs):
        """Dump to file as JSON.
        Convenience method, equivalent to converting to a
        dict and calling :meth:`json.dump`.

        Args:
            fp (file): Write output to this file
            **kwargs: Keywords passed to json.dump()

        Returns:
            see json.dump()
        """
        return json.dump(self.as_dict(), fp, **kwargs)

    def dumps(self, **kwargs):
        """Dump to string as JSON.
        Convenience method, equivalent to converting to a
        dict and calling :meth:`json.dumps`.

        Args:
            **kwargs: Keywords passed to json.dumps()

        Returns:
            (str) JSON-formatted data
        """
        return json.dumps(self.as_dict(), **kwargs)

    def __str__(self):
        return self.dumps()

    @classmethod
    def load(cls, file_or_path, validate=True):
        """Create from JSON input.

        Args:
            file_or_path (file or str): Filename or file object
                from which to read the JSON-formatted data.
            validate (bool): If true, apply validation to input JSON data.

        Example input::

            {
                "meta": [{
                    "datatype": "MEA",
                    "info": "J. Chem. Eng. Data, 2009, Vol 54, pg. 3096-30100",
                    "notes": "r is MEA weight fraction in aqueous soln.",
                    "authors": "Amundsen, T.G., Lars, E.O., Eimer, D.A.",
                    "title": "Density and Viscosity of Monoethanolamine + etc."
                }],
                "data": [
                    {
                        "name": "Viscosity Value",
                        "units": "mPa-s",
                        "values": [2.6, 6.2],
                        "error_type": "absolute",
                        "errors": [0.06, 0.004],
                        "type": "property"
                    }
                ]
            }
        """
        fp = get_file(file_or_path)
        d = json.load(fp)
        if validate:
            cls._validate_json(d)
        metalist = d[Fields.META]
        meta = [Metadata(m) for m in metalist]
        data = TabularData(d[Fields.DATA])
        tbl = Table(data=data)
        for m in meta:
            tbl.add_metadata(m)
        return tbl

    @classmethod
    def _validate_json(cls, d):
        # print('@@ validating:\n----\n{}\n-----'.format(d))
        try:
            cls._validator.validate(d)
        except jsonschema.ValidationError as err:
            raise errors.DataFormatError('tabulardata', str(err))


class TabularData(object):
    """Class representing tabular data that knows how to
    construct itself from a CSV file.

    You can build objects from multiple CSV files as well.
    See the property database section of the API docs for
    details, or read the code in :meth:`add_csv` and the
    tests in :mod:`idaes_dmf.propdb.tests.test_mergecsv`.
    """
    embedded_units = r'(.*)\((.*)\)'

    _validator = jsonschema.Draft4Validator(COLUMN_SCHEMA)

    def __init__(self, data, error_column=False):
        """Construct from a list.

            [  {
                "name": "Density Data",
                "units": "g/cm^3",
                "values": [1.0053, 1.0188, .., ],
                "errors": [.00005, .., .00005],
                "error_type": "absolute"
              },
              ...etc...
            ]

        Args:
            data (list): Input dictionary
            error_column (bool): Whether there are error columns
        Returns:
            TabularData: New instance.
        Raises:
            TypeError: Bad type for `data`
            ValueError: Bad value for `data`
        """
        if not isinstance(data, list):
            raise TypeError('Expected list of dicts, got {}'
                            .format(type(data)))
        if len(data) == 0:
            raise ValueError('Input data must have at least one column')
        try:
            self._validator.validate(data)
        except jsonschema.ValidationError as err:
            raise ValueError(str(err))
        self._data = data
        self._nrows = self._get_nrows()
        self._errcol = error_column

    @property
    def columns(self):
        return self._data

    def __len__(self):
        return self._nrows

    def names(self):
        """Get column names.

        Returns:
            list[str]: List of column names.
        """
        return [v[Fields.DATA_NAME] for v in self.columns]

    @property
    def num_columns(self):
        """Number of columns in this table.

        A "column" is defined as data + error. So if there
        are two columns of data, each with an associated
        error column, then `num_columns` is 2 (not 4).

        Returns:
            int: Number of columns.
        """
        return len(self.columns)

    @property
    def num_rows(self):
        """Number of rows in this table.

        obj.num_rows is a synonym for len(obj)

        Returns:
            int: Number of rows.
        """
        return self._nrows

    def _get_nrows(self):
        n = 0
        for v in self.columns:
            vals = v[Fields.DATA_VALUES]
            if n == 0:
                n = len(vals)
            elif len(vals) != n:
                raise ValueError('Column "{}" length {} != {}'
                                 .format(v[Fields.DATA_NAME], len(vals), n))
        return n

    def get_column(self, key):
        """Get an object for the given named column.

        Args:
            key (str): Name of column

        Returns:
            (TabularColumn) Column object.

        Raises:
            KeyError: No column by that name.
        """
        result = None
        for v in self.columns:
            if v[Fields.DATA_NAME] == key:
                result = Column(key, v)
                break
        if result is None:
            name_list = ', '.join(self.names())
            raise KeyError('Bad column name "{}", not in ({})'.format(
                key, name_list))
        return result

    def get_column_index(self, key):
        """Get an index for the given named column.

        Args:
            key (str): Name of column

        Returns:
            (int) Column number.

        Raises:
            KeyError: No column by that name.
        """
        # print('@@ get column index for name: {}'.format(key))
        for i, v in enumerate(self.columns):
            if v[Fields.DATA_NAME] == key:
                return i
        raise KeyError('Bad column name "{}", not in ({})'.format(
            key, ', '.join(self.names())))

    def as_list(self):
        """Export the data as a list.

        Output will be in same form as data passed to constructor.

        Returns:
            (list) List of dicts
        """
        return self._data

    def as_arr(self):
        """Export property data as arrays.

        Returns:
            (values[M,N], errors[M,N]) Two arrays of floats,
            each with M columns having N values.

        Raises:
            ValueError if the columns are not all the same length
        """
        values, errvals = [], []
        # extract columns
        for v in self.columns:
            values.append(v[Fields.DATA_VALUES])
            errvals.append(v[Fields.DATA_ERRORS])
        return values, errvals

    def values_dataframe(self):
        """Get values as a dataframe.

        Returns:
            (pd.DataFrame) Pandas dataframe for values.

        Raises:
            ImportError: If `pandas` or `numpy` were never
                successfully imported.
        """
        return self._get_dataframe(Fields.DATA_VALUES)

    def errors_dataframe(self):
        """Get errors as a dataframe.

        Returns:
            pd.DataFrame: Pandas dataframe for values.

        Raises:
            ImportError: If `pandas` or `numpy` were never
                successfully imported.
        """
        return self._get_dataframe(Fields.DATA_ERRORS)

    def _get_dataframe(self, field):
        self._check_pandas_import()
        a1, names = [], []
        a1.extend([v[field] for v in self.columns])
        names.extend([v[Fields.DATA_NAME] for v in self.columns])
        a2 = np.array(a1).transpose()
        return pd.DataFrame(a2, columns=names)

    @staticmethod
    def _check_pandas_import():
        if pd is None:
            raise ImportError('Failed to import Pandas and/or Numpy packages '
                              'at module load. Cannot return a Pandas '
                              'Dataframe without Pandas.')

    @staticmethod
    def from_csv(file_or_path, error_column=False):
        """Import the CSV data.

        Expected format of the  files is a header plus data rows.

        Header: Index-column,  Column-name(1), Error-column(1),  \
                Column-name(2), Error-column(2), ..
        Data: <index>, <val>, <errval>, <val>, <errval>, ..

        Column-name is in the format "Name (units)"

        Error-column is in the format "<type> Error", where "<type>" is
        the error type.

        Args:
            file_or_path (file-like or str): Input file
            error_column (bool): If True, look for an error column after each
                           value column. Otherwise, all columns are
                           assumed to be values.

        Returns:
            TabularData: New table of data
        """
        input_file = get_file(file_or_path)
        csv_file = csv.reader(input_file)
        row = next(csv_file)
        names, data = TabularData._parse_csv_headers(row,
                                                     error_column=error_column)
        for row in csv_file:
            # print('@@ parse csv row: {}'.format(row))
            TabularData._parse_csv_row(data, row, error_column=error_column)
        obj = TabularData(data, error_column=error_column)
        return obj

    @classmethod
    def _parse_csv_headers(cls, headers, error_column=None):
        """Parse a row of CSV headers which are pairs
        of columns like "<name> [(units)], <error-type> Error".

        Returns:
             (names, data). Names is a list of all the column names.
                            Data is a list of property/state objects.
        """
        if error_column:
            if len(headers) < 3:
                raise ValueError('Less than 3 columns')
            if len(headers) % 2 != 1:
                raise ValueError('Number of columns must be odd')
        else:
            if len(headers) < 2:
                raise ValueError('Less than 2 columns')
        data = []
        all_names = []
        # Add new item for each value/error pair in column headers
        column_step = 2 if error_column else 1
        for i in range(1, len(headers), column_step):
            errhdr = ''
            if error_column:
                hdr, errhdr = headers[i], headers[i + 1]
            else:
                hdr = headers[i]
            m = re.match(cls.embedded_units, hdr)
            name, units = m.groups() if m else (hdr, '')
            name = name.strip()  # ignore extra ws for column names
            if error_column:
                errtype = errhdr.strip().split()[0].lower()
            else:
                errtype = 'none'
            item = {Fields.DATA_NAME: name,
                    Fields.DATA_UNITS: units,
                    Fields.DATA_VALUES: [],
                    Fields.DATA_ERRORS: [],
                    Fields.DATA_ERRTYPE: errtype}
            data.append(item)
            all_names.append(name)
        return all_names, data

    @classmethod
    def _parse_csv_row(cls, data, row, error_column=None):
        """Add data in row to dict in data, which has the form
        returned by `_parse_csv_headers`.

        Rows are laid out like this:

            id, state-value1, state-error1, state-value2, state-error2, ..., \
            prop-value1, prop-error1, prop-value2, prop-error2, ...

        The number of state-value/error column pairs is equal to `nstates`.
        """
        rowlen = len(row)
        # check that the row has the right number of columns
        if error_column:
            expected_rowlen = 2 * len(data) + 1
        else:
            expected_rowlen = len(data) + 1
        if rowlen != expected_rowlen:
            raise ValueError('CSV row, expected {:d} columns, got {:d}'
                             .format(expected_rowlen, rowlen))
        # iterate over the values in the row, adding each to
        # the appropriate values in data['states'] or data['properties']
        column_step = 2 if error_column else 1
        for i in range(1, rowlen, column_step):
            value = row[i]
            value = float(value) if value else float('nan')
            # add value and error to the property column
            c = data[(i - 1) // column_step]
            c[Fields.DATA_VALUES].append(value)
            if error_column:
                error = row[i + 1]
                error = float(error) if error else float('nan')
                c[Fields.DATA_ERRORS].append(error)


class Metadata(object):
    """Class to import metadata.
    """
    # generic line regular expression
    line_expr = re.compile(r'\s*(\w+)\s*,\s*(.*)\s*')
    # source line regular expression
    source_expr = re.compile(r'\s*(.*)\s*,\s*"(.*)"\s*,\s*(.*)\s*')

    def __init__(self, values=None):
        if values is None:
            values = {}
        self._meta = values

    @property
    def datatype(self):
        return self._meta.get(Fields.DTYPE, '')

    @datatype.setter
    def datatype(self, value):
        self._meta[Fields.DTYPE] = value

    @property
    def author(self):
        """Publication author(s)."""
        return self._meta[Fields.AUTH]

    @property
    def date(self):
        """Publication date"""
        return self._meta.get(Fields.DATE, '1970-01-01')

    @property
    def title(self):
        """Publication title."""
        return self._meta[Fields.TITLE]

    @property
    def info(self):
        """Publication venue, etc."""
        return self._meta[Fields.INFO]

    @property
    def source(self):
        """Full publication info."""
        return '{},"{}",{}'.format(self.author, self.title, self.info)

    def as_dict(self):
        return self._meta

    @staticmethod
    def from_csv(file_or_path):
        """Import metadata from simple text format.

        Example input::

            Source,Han, J., Jin, J., Eimer, D.A., Melaaen, M.C.,"Density of \
            Water(1) + Monoethanolamine(2) + CO2(3) from (298.15 to 413.15) K\
            and Surface Tension of Water(1) + Monethanolamine(2) from ( \
            303.15 to 333.15)K", J. Chem. Eng. Data, 2012, Vol. 57, \
            pg. 1095-1103"
            Retrieval,"J. Morgan, date unknown"
            Notes,r is MEA weight fraction in aqueous soln. (CO2-free basis)

        Args:
            file_or_path (str or file): Input file

        Returns:
            (PropertyMetadata) New instance
        """
        meta = {}
        f = get_file(file_or_path)
        got_source = False
        for line in f:
            line = line.strip()
            if len(line) == 0 or line[0] == '#':
                continue
            m = Metadata.line_expr.match(line)
            if m is None:
                raise ValueError('Cannot parse metadata from: {}'.format(line))
            label, info = m.groups()
            label = label.lower()
            if label == 'source':
                meta.update(Metadata._parse_source(info))
                got_source = True
            else:
                meta[label] = info
        if not got_source:
            raise ValueError('"source" is required')
        return Metadata(meta)

    @staticmethod
    def _parse_source(s):
        m = Metadata.source_expr.match(s)
        if m is None:
            raise ValueError('Cannot parse source: {}'.format(s))
        auth, title, rest = m.groups()
        date_match = re.search(r',\s*(\d\d\d\d(?:-\d\d(?:-\d\d)))', rest)
        date = date_match.group(1) if date_match else '1970-01-01'
        return {Fields.AUTH: auth,
                Fields.DATE: date,
                Fields.TITLE: title,
                Fields.INFO: rest}


class Column(object):
    """Generic, abstract column
    """
    type_name = 'generic'

    def __init__(self, name, data):
        self.name = name
        self.units = data[Fields.DATA_UNITS]
        self.values = data[Fields.DATA_VALUES]
        self.entity_type = self.type_name
        if Fields.DATA_ERRORS in data:
            self.errors = data[Fields.DATA_ERRORS]
            self.error_type = data[Fields.DATA_ERRTYPE]

    def data(self):
        return {
            Fields.DATA_UNITS: self.units,
            Fields.DATA_VALUES: self.values,
            Fields.DATA_ERRORS: self.errors,
            Fields.DATA_ERRTYPE: self.error_type
        }

    def __len__(self):
        return len(self.values)
