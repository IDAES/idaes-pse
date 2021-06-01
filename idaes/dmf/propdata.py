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
Property data types.

Ability to import, etc. from text files is
part of the methods in the type.

Import property database from textfile(s):
* See :meth:`PropertyData.from_csv`, for the expected format for data.
* See :meth:`PropertyMetadata()` for the expected format for metadata.

"""
# stdlib
import csv
import json
import logging

# third-party
try:
    import pandas as pd
    import numpy as np
except ImportError:
    np, pd = None, None

# local
from .util import get_file
from . import tabular

__author__ = 'Dan Gunter'

_log = logging.getLogger(__name__)


class AddedCSVColumnError(KeyError):
    """Error for :meth:PropertyData.add_csv()
    """

    def __init__(self, names, how_bad, column_type=''):
        ctype = column_type + ' ' if column_type else ''
        if len(names) == 1:
            msg = 'Added CSV data {} {}column "{}"'.format(
                how_bad, ctype, list(names)[0]
            )
        else:
            msg = 'Added CSV data {} {}columns: {}'.format(
                how_bad, ctype, ', '.join(list(names))
            )
        KeyError.__init__(self, msg)


class Fields(tabular.Fields):
    """Constants for fields.
    """

    # Values for "type" field
    C_STATE, C_PROP = 'state', 'property'


class PropertyTable(tabular.Table):
    """Property data and metadata together (at last!)
    """

    def __init__(self, data=None, **kwargs):
        """Constructor.
        """
        if isinstance(data, PropertyData):
            pdata = data
        elif isinstance(data, list):
            pdata = PropertyData(data)
        else:
            raise TypeError('list or PropertyData object required')
        super(PropertyTable, self).__init__(data=pdata, **kwargs)

    @classmethod
    def load(cls, file_or_path, validate=True):
        """Create PropertyTable from JSON input.

        Args:
            file_or_path (file or str): Filename or file object
                from which to read the JSON-formatted data.
            validate (bool): If true, apply validation to input JSON data.

        Example input::

            {
                "meta": [
                    {"datatype": "MEA",
                     "info": "J. Chem. Eng. Data, 2009, Vol 54, pg. 306-310",
                     "notes": "r is MEA weight fraction in aqueous soln.",
                     "authors": "Amundsen, T.G., Lars, E.O., Eimer, D.A.",
                     "title": "Density and Viscosity of ..."}
                ],
                "data": [
                    {"name": "Viscosity Value",
                     "units": "mPa-s",
                     "values": [2.6, 6.2],
                     "error_type": "absolute",
                     "errors": [0.06, 0.004],
                     "type": "property"},
                    {"name": "r",
                     "units": "",
                     "values": [0.2, 1000],
                     "type": "state"}
                ]
            }
        """
        fp = get_file(file_or_path)
        d = json.load(fp)
        PropertyTable._validate_json(d)
        metalist = d[Fields.META]
        meta = [PropertyMetadata(m) for m in metalist]
        data = PropertyData(d[Fields.DATA])
        tbl = PropertyTable(data=data)
        for m in meta:
            tbl.add_metadata(m)
        return tbl


class PropertyData(tabular.TabularData):
    """Class representing property data that knows how to
    construct itself from a CSV file.

    You can build objects from multiple CSV files as well.
    See the property database section of the API docs for
    details, or read the code in :meth:`add_csv` and the
    tests in :mod:`idaes_dmf.propdb.tests.test_mergecsv`.
    """

    embedded_units = r'(.*)\((.*)\)'

    def __init__(self, data):
        """Construct new object from input list.

        Example input::

            [{
                "name": "Density Data",
                "units": "g/cm^3",
                "values": [1.0053, 1.0188, .., ],
                "errors": [.00005, .., .00005],
                "error_type": "absolute",
                "type": "property"
            }, ...etc...]

        Args:
            data (list): Input data columns

        Returns:
            (PropertyData) New instance.
        """
        super(PropertyData, self).__init__(data, error_column=True)
        self._nstates = len(self.states)

    @property
    def states(self):
        return [c for c in self.columns if self._is_state(c)]

    @property
    def properties(self):
        return [c for c in self.columns if self._is_prop(c)]

    @staticmethod
    def _is_state(c):
        return c[Fields.COLTYPE] == Fields.C_STATE

    @staticmethod
    def _is_prop(c):
        return c[Fields.COLTYPE] == Fields.C_PROP

    def names(self, states=True, properties=True):
        """Get column names.

        Args:
            states (bool): If False, exclude "state" data, e.g. the
                          ambient temperature, and only
                          include measured property values.
            properties (bool): If False, excluse property data

        Returns:
            list[str]: List of column names.
        """
        result = []
        if states:
            result.extend([v[Fields.DATA_NAME] for v in self.states])
        if properties:
            result.extend([v[Fields.DATA_NAME] for v in self.properties])
        return result

    def is_state_column(self, index):
        """Whether given column is state.

        Args:
            index (int): Index of column
        Returns:
            (bool) State or property and the column number.
        Raises:
            IndexError: No column at that index.
        """
        col = self.columns[index]
        return self._is_state(col)

    def is_property_column(self, index):
        """Whether given column is a property. See :meth:`is_state_column`."""
        return not self.is_state_column(index)

    def as_arr(self, states=True):
        """Export property data as arrays.

        Args:
            states (bool): If False, exclude "state" data, e.g. the
                          ambient temperature, and only
                          include measured property values.
        Returns:
            (values[M,N], errors[M,N]) Two arrays of floats,
            each with M columns having N values.
        Raises:
            ValueError if the columns are not all the same length
        """
        n, values, errors = None, [], []
        # extract state columns
        if states:
            for v in self.states:
                vals = v[Fields.DATA_VALUES]
                if n is None:
                    n = len(vals)
                elif len(vals) != n:
                    raise ValueError(
                        'State values "{}" length {} != {}'.format(
                            v[Fields.DATA_NAME], len(vals), n
                        )
                    )
                values.append(vals)
                errors.append([0] * len(vals))
        # extract property columns
        for v in self.properties:
            vals = v[Fields.DATA_VALUES]
            if n is None:
                n = len(vals)
            elif len(vals) != n:
                raise ValueError(
                    'Property values "{}" length {} != {}'.format(
                        v[Fields.DATA_NAME], len(vals), n
                    )
                )
            values.append(v[Fields.DATA_VALUES])
            errors.append(v[Fields.DATA_ERRORS])
        return values, errors

    def values_dataframe(self, states=True):
        """Get values as a dataframe.

        Args:
            states (bool): see :meth:`names()`.

        Returns:
            (pd.DataFrame) Pandas dataframe for values.
 
        Raises:
            ImportError: If `pandas` or `numpy` were never
                successfully imported.
        """
        return self._get_prop_dataframe(Fields.DATA_VALUES, states)

    def errors_dataframe(self, states=False):
        """Get errors as a dataframe.

        Args:
            states (bool): If False, exclude state data.
                           This is the default, because states do not
                           normally have associated error information.
        Returns:
            pd.DataFrame: Pandas dataframe for values.

        Raises:
            ImportError: If `pandas` or `numpy` were never
                successfully imported.
        """
        return self._get_prop_dataframe(Fields.DATA_ERRORS, states)

    def _get_prop_dataframe(self, field, states):
        self._check_pandas_import()
        a1, names = [], []
        if states:
            a1 = [v[field] for v in self.states]
            names = [v[Fields.DATA_NAME] for v in self.states]
        a1.extend([v[field] for v in self.properties])
        names.extend([v[Fields.DATA_NAME] for v in self.properties])
        a2 = np.array(a1).transpose()
        return pd.DataFrame(a2, columns=names)

    @staticmethod
    def from_csv(file_or_path, nstates=0):
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
            nstates (int): Number of state columns, appearing
               first before property columns.

        Returns:
            PropertyData: New properties instance
        """
        input_file = get_file(file_or_path)
        csv_file = csv.reader(input_file)
        row = next(csv_file)
        names, data = PropertyData._prop_parse_csv_headers(nstates, row)
        for row in csv_file:
            # print('@@ parse csv row: {}'.format(row))
            PropertyData._parse_csv_row(data, row, error_column=True)
        obj = PropertyData(data)
        return obj

    def add_csv(self, file_or_path, strict=False):
        """Add to existing object from a new CSV file.

        Depending on the value of the `strict` argument (see
        below), the new file may or may not have the same 
        properties as the object -- but it always needs to have
        the same number of state columns, and in the same order.

        .. note:: Data that is "missing" because of property columns in
           one CSV and not the other will be filled with `float(nan)` values.
        
        Args:
            file_or_path (file or str): Input file. This should be in exactly
                the same format as expected by :meth:from_csv().
            strict (bool): If true, require that the columns in the input
                CSV match columns in this object. Otherwise, only require
                that *state* columns in input CSV match columns in this
                object. New property columns are added, and matches
                to existing property columns will append the data.
                
            
        Raises:
            AddedCSVColumnError: If the new CSV column headers are not the
                same as the ones in this object.
            
        Returns:
            (int) Number of added rows
        """
        nstates = self._nstates
        input_file = get_file(file_or_path)
        csv_file = csv.reader(input_file)

        # Parse the header
        row = next(csv_file)
        hdr_names, hdr_data = PropertyData._prop_parse_csv_headers(nstates, row)

        # print('@@ add_csv, column names = {}, data columns = {}'
        #      .format(hdr_names, self.names()))

        # Check that set of keys in new data is the same
        cur_keys = set(self.names())
        new_keys = set(hdr_names)
        # This is used to re-order input data
        rowmap = None
        if strict:
            if cur_keys > new_keys:
                missing = cur_keys - new_keys
                raise AddedCSVColumnError(missing, 'is missing')
            elif new_keys > cur_keys:
                extra = new_keys - cur_keys
                raise AddedCSVColumnError(extra, 'has extra')
            elif new_keys != cur_keys:
                extra = new_keys - cur_keys
                missing = cur_keys - new_keys
                namelist = (
                    '(' + ','.join(extra) + ')',
                    'instead of',
                    '(' + ','.join(missing) + ')',
                )
                raise AddedCSVColumnError(namelist, 'has different')
        else:
            # check that all states are in common
            hdr_states = filter(self._is_state, hdr_data)
            new_states = [s[Fields.DATA_NAME] for s in hdr_states]
            new_states = set(new_states)
            cur_states = set(self.names(properties=False))
            if new_states != cur_states:
                extra = new_states - cur_states
                missing = cur_states - new_states
                if extra and missing:
                    namelist = (
                        '(' + ','.join(extra) + ')',
                        'instead of',
                        '(' + ','.join(missing) + ')',
                    )
                    raise AddedCSVColumnError(
                        namelist, 'has different', column_type='state'
                    )
                elif extra:
                    raise AddedCSVColumnError(extra, 'has extra', column_type='state')
                elif missing:
                    raise AddedCSVColumnError(
                        missing, 'is missing', column_type='state'
                    )
                else:
                    raise RuntimeError('unexpected branch')
            # check that at least one property is in common
            new_prop = new_keys - new_states
            if not new_prop:
                return 0  # no data
            cur_prop = set(self.names(states=False))

            # Add columns for all properties only found on the input,
            # and initialize values to a list of NaN's as long as the
            # current table, so data in all fields will be the same length.

            # Initialize rowmap with mapping for state columns
            rowmap = [-1] * len(hdr_names)
            idx = 0
            for i, s in enumerate(hdr_data):
                if s[Fields.COLTYPE] == Fields.C_PROP:
                    continue
                rowmap[i] = idx
                idx += 1

            nan_list = [float('nan')] * self.num_rows
            idx = 0
            for i, value in enumerate(hdr_data):
                if value[Fields.COLTYPE] == Fields.C_STATE:
                    continue
                name = value[Fields.DATA_NAME]
                if name not in cur_prop:
                    value[Fields.DATA_NAME] = name
                    value[Fields.DATA_VALUES] = nan_list[:]
                    value[Fields.DATA_ERRORS] = nan_list[:]
                    value[Fields.COLTYPE] = Fields.C_PROP
                    self._data.append(value)
                    rowmap[i] = len(self.properties) - 1
                else:
                    rowmap[i] = idx + self._nstates
                idx += 1

        # print("@@ rowmap = {}".format(rowmap))

        # Parse the new data
        num_added = 0
        new_rowlen = 1 + 2 * len(self.names())
        for row in csv_file:
            if rowmap:
                # Re-order according to the rowmap.
                # By initializing with NaN, any columns not in the
                # input, but in the current data, will be replaced with NaN
                # values.
                row2 = [float('nan')] * new_rowlen
                # print('@@ row={} row2-init={}'.format(row, row2))
                for i, j in enumerate(rowmap):
                    row2[j * 2 + 1] = row[i * 2 + 1]  # value
                    row2[j * 2 + 2] = row[i * 2 + 2]  # error
                row = row2
            self._parse_csv_row(self._data, row, error_column=True)
            num_added += 1
            self._nrows += 1

        return num_added

    @classmethod
    def _prop_parse_csv_headers(cls, nstates, headers):
        """Parse a row of CSV headers which are pairs
        of columns like "<name> [(units)], <error-type> Error".

        Returns:
             (names, data). Names is a list of all the column names.
             Data is a dict with two keys, "properties" and "states".
             Each value will be a list of property/state objects.
        """
        names, data = cls._parse_csv_headers(headers, error_column=True)
        for i in range(0, nstates):
            data[i][Fields.COLTYPE] = Fields.C_STATE
        for i in range(nstates, len(data)):
            data[i][Fields.COLTYPE] = Fields.C_PROP
        return names, data


class PropertyMetadata(tabular.Metadata):
    """Class to import property metadata.
    """

    pass


class PropertyColumn(tabular.Column):
    """Data column for a property.
    """

    type_name = 'Property'

    def __init__(self, name, data):
        tabular.Column.__init__(self, name, data)
        self.errors = data[Fields.DATA_ERRORS]
        self.error_type = data[Fields.DATA_ERRTYPE]

    def data(self):
        return {
            Fields.DATA_UNITS: self.units,
            Fields.DATA_VALUES: self.values,
            Fields.DATA_ERRORS: self.errors,
            Fields.DATA_ERRTYPE: self.error_type,
        }


class StateColumn(tabular.Column):
    """Data column for a state.
    """

    type_name = 'State'

    def __init__(self, name, data):
        tabular.Column.__init__(self, name, data)
        self.errors = [0.0] * len(self)
        self.error_type = 'none'

    def data(self):
        return {Fields.DATA_UNITS: self.units, Fields.DATA_VALUES: self.values}


def convert_csv(meta_csv, datatype, data_csv, nstates, output):
    meta = PropertyMetadata.from_csv(meta_csv)
    meta.datatype = datatype
    data = PropertyData.from_csv(data_csv, nstates)
    obj = PropertyTable(data=data, metadata=meta)
    ofile = get_file(output, mode='w')
    obj.dump(ofile)
