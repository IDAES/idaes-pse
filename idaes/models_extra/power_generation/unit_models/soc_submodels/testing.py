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

__author__ = "Douglas Allan"

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock


def _cell_flowsheet_model(dynamic, time_set, zfaces):
    # function that creates a unit model with cell-level variables for testing
    # subcomponents that require them
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=time_set, time_units=pyo.units.s)
    tset = m.fs.config.time
    znodes = m.fs.znodes = pyo.Set(
        initialize=[(zfaces[i] + zfaces[i + 1]) / 2.0 for i in range(len(zfaces) - 1)]
    )
    iznodes = m.fs.iznodes = pyo.Set(initialize=range(1, len(znodes) + 1))

    m.fs.length_z = pyo.Var(initialize=0.25, units=pyo.units.m)
    m.fs.length_y = pyo.Var(initialize=0.25, units=pyo.units.m)

    m.fs.current_density = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.A / pyo.units.m**2
    )

    m.fs.temperature_z = pyo.Var(tset, iznodes, initialize=1000, units=pyo.units.K)

    m.fs.length_y.fix(0.08)
    m.fs.length_z.fix(0.08)
    m.fs.temperature_z.fix(1000)
    m.fs.current_density.fix(0)

    return m


def _build_test_utility(block, comp_dict, references=None):
    # Takes a unit model and four dictionaries: references (variables),
    # not_references (variables), constraints, and expressions. They should
    # have the attribute name as a key and the length of the attribute as
    # the value. This function goes through and ensures that all these
    # components exist and are the right length, and furthermore they are no
    # unlisted components of these types

    if references is not None:
        for attr in references:
            try:
                comp = getattr(block, attr)
            except AttributeError:
                raise AttributeError(
                    f"Reference {attr} missing from block {block.name}."
                )
            if not comp.is_reference():
                raise AssertionError(
                    f"Attribute {attr} found on block {block.name}, but "
                    "was not Reference."
                )
        for comp in block.component_data_objects(descend_into=False):
            if comp.is_reference():
                if not comp in references:
                    raise AssertionError(
                        f"Unexpected Reference {comp.name} encountered "
                        f"in block {block.name}."
                    )
    for ctype, sub_dict in comp_dict.items():
        for attr, length in sub_dict.items():
            try:
                comp = getattr(block, attr)
            except AttributeError:
                raise AttributeError(f"{ctype} {attr} missing from block {block.name}.")

            if not len(comp) == length:
                raise AssertionError(
                    f"{ctype} {comp.name} was not expected length in block "
                    f"{block.name}."
                )
        for comp in block.component_data_objects(ctype=ctype, descend_into=False):
            short_name = comp.local_name.split("[")[0]
            if not short_name in sub_dict.keys():
                raise AssertionError(
                    f"Unexpected {ctype} {comp.name} encountered in block "
                    f"{block.name}."
                )
