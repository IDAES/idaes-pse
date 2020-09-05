# -*- coding: UTF-8 -*-
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
A module of functions and classes for configuring NMPC/MHE problems
"""
import enum
from pyomo.environ import SolverFactory
from pyomo.core.base.var import _GeneralVarData
from pyomo.common.config import ConfigEnum


class ControlInitOption(ConfigEnum):
    FROM_PREVIOUS = 11
    BY_TIME_ELEMENT = 12
    FROM_INITIAL_CONDITIONS = 13


class ElementInitializationInputOption(ConfigEnum):
    SET_POINT = 21
    INITIAL = 22
    CURRENT_VALUES = 23


class TimeResolutionOption(ConfigEnum):
    COLLOCATION_POINTS = 31
    FINITE_ELEMENTS = 32
    SAMPLE_POINTS = 33
    INITIAL_POINT = 34


class ControlPenaltyType(ConfigEnum):
    ERROR = 41
    ACTION = 42
    NONE = 43


class VariableCategory(ConfigEnum):
    DIFFERENTIAL = 51
    ALGEBRAIC = 52
    DERIVATIVE = 53
    INPUT = 54
    FIXED = 55
    SCALAR = 56

# This function is used as the domain for the user-provided
# list of inputs at time.first().
def validate_list_of_vardata(varlist):
    if not isinstance(varlist, list):
        raise TypeError('Not a list of VarData')
    for var in varlist:
        if not isinstance(var, _GeneralVarData):
            raise TypeError('Not a list of VarData')
    return varlist


def validate_list_of_vardata_value_tuples(varvaluelist):
    if not isinstance(varvaluelist, list):
        raise TypeError('Not a list')
    for item in varvaluelist:
        if not isinstance(item, tuple):
            raise TypeError('Item in list is not a tuple')
        if not len(item) == 2:
            raise ValueError('Tuple in list does not have correct length')
        if not isinstance(item[0], _GeneralVarData):
            raise TypeError('First entry is not a VarData')
        item = (item[0], float(item[1]))
    return varvaluelist


def validate_solver(solver):
    if type(solver) is str:
        solver = SolverFactory(solver)
    if not hasattr(solver, 'solve'):
        raise TypeError(
            'Solver does not implement solve method')
    return solver


def get_ncp(continuous_set):
    scheme = continuous_set.get_discretization_info()['scheme']
    if scheme == 'LAGRANGE-RADAU':
        return continuous_set.get_discretization_info()['ncp']
    elif scheme == 'BACKWARD Difference':
        return 1
    else:
        raise NotImplementedError(
            '%s discretization scheme is not supported' % scheme)

