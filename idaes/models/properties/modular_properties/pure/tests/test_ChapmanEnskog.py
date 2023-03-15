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
Tests for Lennard-Jones implementation of Chapman Enskog viscosity from
Properties of Gases and Liquids, 5th Ed., 9-4-1

Authors: Douglas Allan
"""

import pytest
import types

import pyomo.environ as pyo
from pyomo.environ import ConcreteModel, Block, value, units as pyunits, Var
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.ChapmanEnskog import (
    ChapmanEnskogLennardJones,
    collision_integral_neufeld_callback,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata


def construct_dummy_model(param_dict):
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    # Parameters for sulfur dioxide, from Properties of Gases and Liquids 5th ed. Appendix B.
    m.params.config.parameter_data = {
        "lennard_jones_sigma": (param_dict["sigma"], pyunits.angstrom),
        "lennard_jones_epsilon_reduced": (
            param_dict["epsilon"],
            pyunits.K,
        ),
    }

    m.meta_object = PropertyClassMetadata()
    m.meta_object._default_units.set_units(
        temperature=pyunits.K,
        mass=pyunits.kg,
        length=pyunits.m,
        time=pyunits.s,
        amount=pyunits.mol,
    )

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Create variables that should exist on param block
    m.params.mw = Var(initialize=param_dict["mw"], units=pyunits.kg / pyunits.mol)
    m.params.omega = Var(initialize=param_dict["omega"], units=pyunits.dimensionless)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=273.16, units=pyunits.K)
    return m


@pytest.mark.unit
def test_visc_vap_comp_sulfur_dioxide():
    m = construct_dummy_model(
        {
            "sigma": 4.112,
            "epsilon": 335.4,
            "mw": 0.064065,
            "omega": 0.257,
        }
    )

    ChapmanEnskogLennardJones.visc_d_phase_comp.build_parameters(m.params, "Vap")

    assert isinstance(m.params.lennard_jones_sigma, Var)
    assert value(m.params.lennard_jones_sigma) == pytest.approx(4.112e-10, rel=1e-12)
    assert isinstance(m.params.lennard_jones_epsilon_reduced, Var)
    assert value(m.params.lennard_jones_epsilon_reduced) == 335.4
    assert (
        m.params.viscosity_collision_integral_callback
        is collision_integral_neufeld_callback
    )

    expr = ChapmanEnskogLennardJones.visc_d_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)

    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    T_list = [10, 100, 300, 700]  # Gas temperature in C
    visc_list = [120, 163, 246, 376]  # Experimental viscosities in millipoise
    for i in range(4):
        m.props[1].temperature.value = T_list[i] + 273.15
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_visc_vap_comp_methanol():
    m = construct_dummy_model(
        {
            "sigma": 3.626,
            "epsilon": 481.8,
            "mw": 32.042e-3,
            "omega": 0.565,
        }
    )
    ChapmanEnskogLennardJones.visc_d_phase_comp.build_parameters(m.params, "Vap")

    expr = ChapmanEnskogLennardJones.visc_d_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids
    T_list = [67, 127, 277]  # Gas temperature in C
    visc_list = [112, 132, 181]  # Experimental viscosities in millipoise
    for i in range(3):
        m.props[1].temperature.value = T_list[i] + 273.15
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_visc_vap_comp_ethane():
    m = construct_dummy_model(
        {
            "sigma": 4.443,
            "epsilon": 215.7,
            "mw": 30.070e-3,
            "omega": 0.099,
        }
    )
    ChapmanEnskogLennardJones.visc_d_phase_comp.build_parameters(m.params, "Vap")

    expr = ChapmanEnskogLennardJones.visc_d_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids
    T_list = [47, 117, 247]  # Gas temperature in C
    visc_list = [100, 120, 156]  # Experimental viscosities in millipoise

    for i in range(3):
        m.props[1].temperature.value = T_list[i] + 273.15
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)
