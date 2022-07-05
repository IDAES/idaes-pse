# -*- coding: UTF-8 -*-
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
A module of functions and classes for configuring NMPC/MHE problems
"""
import enum
from pyomo.environ import SolverFactory
from pyomo.core.base.var import _GeneralVarData


class ControlInitOption(enum.Enum):
    FROM_PREVIOUS = 11
    BY_TIME_ELEMENT = 12
    FROM_INITIAL_CONDITIONS = 13
    SETPOINT = 14


class ElementInitializationInputOption(enum.Enum):
    SETPOINT = 21
    INITIAL = 22
    CURRENT_VALUES = 23


class InputOption(enum.Enum):
    CURRENT = 20
    INITIAL = 21
    SETPOINT = 22
    # TODO: PREVIOUS = 23


class TimeResolutionOption(enum.Enum):
    COLLOCATION_POINTS = 31
    FINITE_ELEMENTS = 32
    SAMPLE_POINTS = 33
    INITIAL_POINT = 34


class ControlPenaltyType(enum.Enum):
    ERROR = 41
    ACTION = 42
    NONE = 43


class VariableCategory(enum.Enum):
    DIFFERENTIAL = 1
    ALGEBRAIC = 2
    DERIVATIVE = 3
    INPUT = 4
    FIXED = 5
    SCALAR = 6
    UNUSED = 7
    DISTURBANCE = 8
    MEASUREMENT = 9


class ConstraintCategory(enum.Enum):
    DIFFERENTIAL = 1
    ALGEBRAIC = 2
    DISCRETIZATION = 3
    INPUT = 4
    INITIAL = 5
    TERMINAL = 6
    SENSITIVITY = 7
    SCALAR = 8
    UNUSED = 9
    INEQUALITY = 10


class NoiseBoundOption(enum.Enum):
    FAIL = 60
    DISCARD = 61
    PUSH = 62


class PlantHorizonType(enum.Enum):
    FULL = 71
    ROLLING = 72


# This function is used as the domain for the user-provided
# list of inputs at time.first().
def validate_list_of_vardata(varlist):
    if not isinstance(varlist, list):
        raise TypeError("Not a list of VarData")
    for var in varlist:
        if not isinstance(var, _GeneralVarData):
            raise TypeError("Not a list of VarData")
    return varlist


def validate_list_of_vardata_value_tuples(varvaluelist):
    if not isinstance(varvaluelist, list):
        raise TypeError("Not a list")
    for item in varvaluelist:
        if not isinstance(item, tuple):
            raise TypeError("Item in list is not a tuple")
        if not len(item) == 2:
            raise ValueError("Tuple in list does not have correct length")
        if not isinstance(item[0], _GeneralVarData):
            raise TypeError("First entry is not a VarData")
        item = (item[0], float(item[1]))
    return varvaluelist


def validate_solver(solver):
    if type(solver) is str:
        solver = SolverFactory(solver)
    if not hasattr(solver, "solve"):
        raise TypeError("Solver does not implement solve method")
    return solver


def get_ncp(continuous_set):
    scheme = continuous_set.get_discretization_info()["scheme"]
    if scheme == "LAGRANGE-RADAU":
        return continuous_set.get_discretization_info()["ncp"]
    elif scheme == "BACKWARD Difference":
        return 1
    else:
        raise NotImplementedError("%s discretization scheme is not supported" % scheme)
