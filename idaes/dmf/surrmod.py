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
Surrogate modeling helper classes and functions.
This is used to run ALAMO on property data.
"""
# stdlib
import logging
# third-party
import numpy as np
from pandas import DataFrame
# package
from idaes.dmf import resource, propdata
from idaes.dmf.resource import Predicates
from idaes.dmf import errors
# alamo
from idaes.surrogate import alamopy

__author__ = 'Dan Gunter <dkgunter@lbl.gov>'

_log = logging.getLogger(__name__)

alamo_enabled = alamopy.multos.has_alamo()


class SurrogateModel(object):
    """Run ALAMO to generate surrogate models.

    Automatically track the objects in the DMF.

    Example::

        model = SurrogateModel(dmf, simulator='linsim.py')
        rsrc = dmf.fetch_one(1) # get resource ID 1
        data = rsrc.property_table.data
        model.set_input_data(data, ['temp'], 'density')
        results = model.run()
    """

    PARAM_DATA_KEY = 'parameters'  #: Key in resource 'data' for params

    def __init__(self, experiment, **kwargs):
        """Create surrogate model generator.

        This invokes `alamopy.doalamo()` when the `model` attribute
        is retrieved. Subsequent calls return the model that was generated
        by the first call.

        Args:
            experiment (Experiment): Associated parent resource.
            **kwargs: Keyword arguments passed to doalamo()
        Raises:
            errors.AlamoDisabledError: if alamopy cannot be imported or used.
        """
        if not alamo_enabled:
            raise errors.AlamoDisabledError()
        self._x, self._z, self._xv, self._zv = None, None, None, None
        self._kwargs = kwargs
        self._exp = experiment
        self._rsrc = self._create_resource()
        self._exp.add(self._rsrc)

    def set_input_data(self, data, x_colnames, z_colname):
        """Set input from provided dataframe or property data.

        Args:
            data (PropertyData|pandas.DataFrame): Input data
            x_colnames (List[str]|str): One or more column names for parameters
            z_colname (str): Column for response variable
        Returns:
            None
        Raises:
            KeyError: if columns are not found in data
        """
        df = self._get_df(data)
        self._x = np.array(df[x_colnames])
        self._z = np.array(df[z_colname])
        self._add_data(data, dtype='input', x_colnames=x_colnames,
                       z_colname=z_colname)
        if 'xlabels' not in self._kwargs:
            self._kwargs['xlabels'] = x_colnames

    def set_input_data_np(self, x, z, xlabels=None, zlabel='z'):
        """Set input data from numpy arrays.

        Args:
            x (arr): Numpy array with parameters
            xlabels (List[str]): List of labels for x
            zlabel (str): Label for z
            z (arr): Numpy array with response variables

        Returns:
            None
        """
        self._x, self._z = x, z
        xlabels = self._get_xlabels(x, xlabels)
        self._add_data(self._make_dataframe(x, z, xlabels, zlabel),
                       dtype='input', x_colnames=xlabels, z_colname=zlabel)

    def set_validation_data(self, data, x_colnames, z_colname):
        """Set validation data from provided data.

        Args:
            data (PropertyData|pandas.DataFrame): Input data
            x_colnames (List[str]|str): One or more column names for parameters
            z_colname (str): Column for response variable
        Returns:
            None
        Raises:
            KeyError: if columns are not found in data
        """
        df = self._get_df(data)
        self._xv = np.array(df[x_colnames])
        self._zv = np.array(df[z_colname])
        self._add_data(data, dtype='validation', x_colnames=x_colnames,
                       z_colname=z_colname)

    def set_validation_data_np(self, x, z, xlabels=None, zlabel='z'):
        """Set input data from numpy arrays.

        Args:
            x (arr): Numpy array with parameters
            xlabels (List[str]): List of labels for x
            zlabel (str): Label for z
            z (arr): Numpy array with response variables

        Returns:
            None
        """
        self._xv, self._zv = x, z
        xlabels = self._get_xlabels(x, xlabels)
        df = self._make_dataframe(x, z, xlabels, zlabel)
        self._add_data(df, dtype='validation', x_colnames=xlabels,
                       z_colname=zlabel)

    @staticmethod
    def _get_xlabels(x, xlabels):
        if xlabels is None:
            # generate labels x1, x2, .., xN for each column in x
            xlabels = ['x{:d}'.format(i) for i in range(x.shape[1])]
        return xlabels

    @staticmethod
    def _make_dataframe(x, z, xlabels, zlabel=None):
        """Make a dataframe from some numpy inputs.

        Args:
            x (np.Array): Table of x (variables)
            z (np.Array): Vector of z (responses)
            xlabels (List[str]): Names for each x
            zlabel (str): Name for response

        Returns:
            DataFrame: DataFrame with appropriate header
        """
        df = DataFrame(x, columns=xlabels)
        df[zlabel] = z

        return df

    def run(self, **kwargs):
        """Run ALAMO.

        Args:
            **kwargs: Additional arguments merged with those passed to
                      the class constructor. Any duplicate values will
                      override the earlier ones.

        Returns:
            dict: The dictionary returned from :meth:`alamopy.doalamo`
        """
        if kwargs:
            self._kwargs.update(kwargs)
            # update stored resource
            self._rsrc.data = {self.PARAM_DATA_KEY: self._kwargs}
            self._exp.dmf.update(self._rsrc)
        results = self._run_alamo()
        # TODO: add result object to DMF
        return results

    def _run_alamo(self):
        kwargs = self._kwargs.copy()
        if self._xv is not None and self._zv is not None:
            kwargs['xval'], kwargs['zval'] = self._xv, self._zv
        # run alamo to generate model
        try:
            results = alamopy.doalamo(self._x, self._z, **kwargs)
        except alamopy.AlamoError as err:
            raise errors.AlamoError(str(err))
        # return the generated model
        return results

    @staticmethod
    def _get_df(data):
        """Return input, if already a dataframe, or extract one."""
        if isinstance(data, propdata.PropertyData):
            df = data.values_dataframe()
        elif isinstance(data, dict):
            df = DataFrame(data)
        else:
            df = data  # assume dataframe-like object
        return df

    def _add_data(self, data, dtype='', x_colnames=None, z_colname=None,
                  metadata=None):
        """Add the data as a DMF resource.

        Args:
            data (propdata.PropertyData|DataFrame): Data values
            metadata (propdata.Metadata): If present, metadata
            dtype (str): Either 'input' or 'validation', the type of data
            x_colnames (List[str]): List of x (variable) column names
            z_colname (str): z (response) column name
        """
        if isinstance(data, propdata.PropertyData):
            pdata = data.as_arr()
        else:
            # create property data from dataframe
            assert x_colnames is not None
            assert z_colname is not None
            pdata = []
            for col in x_colnames:
                values = list(data[col])
                pdata.append({"name": "{} x:{}".format(dtype, col),
                              "units": "",
                              "values": values,
                              "errors": [0] * len(values),
                              "error_type": "absolute",
                              "type": "property"})
            values = list(data[z_colname])
            pdata.append({'name': '{} z:{}'.format(dtype, z_colname),
                          'units': '',
                          'values': values,
                          'errors': [0] * len(values),
                          'error_type': 'absolute',
                          'type': 'property'})
        # create resource
        r = resource.Resource(type_=resource.ResourceTypes.property)
        r.data = {'data': pdata,
                  'meta': metadata.as_dict() if metadata else {}}
        r.v['aliases'].append(dtype)
        _log.debug('adding resource dtype={}'.format(dtype))
        # add resource to experiment
        self._exp.add(r)
        # add link between resource and the surrogate-model resource
        self._exp.link(r, Predicates.uses, self._rsrc)

    def _create_resource(self):
        r = resource.Resource(type_=resource.ResourceTypes.surrogate_model)
        r.v['desc'] = 'SurrogateModel'
        r.data = {self.PARAM_DATA_KEY: self._kwargs}
        return r
