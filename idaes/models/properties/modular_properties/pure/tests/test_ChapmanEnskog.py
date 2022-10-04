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
Tests for methods from Perry's

All methods and parameters from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

All parameter indicies based on conventions used by the source

Authors: Andrew Lee
"""

import pytest
import types

from pyomo.environ import ConcreteModel, Block, value, Var, Param, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.pure.ChapmanEnskog import *
from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

@pytest.mark.unit
def test_visc_vap_comp_sulfur_dioxide():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    # Parameters for sulfur dioxide, from Properties of Gases and Liquids 5th ed. Appendix B.
    m.params.config.parameter_data = {
        "lennard_jones_sigma": 4.112 * pyunits.angstrom,
        "lennard_jones_epsilon_reduced": 335.4 * pyunits.K,
    }

    # Also need to dummy configblock on the model for the test
    # m.config = ConfigBlock(implicit=True)
    # m.config.include_enthalpy_of_formation = True

    m.meta_object = PropertyClassMetadata()
    m.meta_object.default_units["temperature"] = pyunits.K
    m.meta_object.default_units["mass"] = pyunits.kg
    m.meta_object.default_units["length"] = pyunits.m
    m.meta_object.default_units["time"] = pyunits.s
    m.meta_object.default_units["amount"] = pyunits.mol
    m.meta_object.default_units["current"] = pyunits.ampere

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Create variables that should exist on param block
    m.params.mw = Var(initialize=0.064065, units=pyunits.kg / pyunits.mol)
    m.params.omega = Var(initialize=0.257, units=pyunits.dimensionless)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=273.16, units=pyunits.K)

    ChapmanEnskogLennardJones.visc_d_phase_comp.build_parameters(m.params, "Vap")

    assert isinstance(m.params.lennard_jones_sigma, Var)
    assert value(m.params.lennard_jones_sigma) == pytest.approx(4.112e-10, rel=1e-12)
    assert isinstance(m.params.lennard_jones_epsilon_reduced, Var)
    assert value(m.params.lennard_jones_epsilon_reduced) == 335.4
    assert m.params.viscosity_collision_integral_callback is collision_integral_neufeld_callback

    expr = ChapmanEnskogLennardJones.visc_d_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)

    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    T_list = [10, 100, 300, 700]  # Gas temperature in C
    visc_list = [120, 163, 246, 376]  # Experimental viscosities in millipoise
    for i in range(4):
        m.props[1].temperature.value = (T_list[i] + 273.15)
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_visc_vap_comp_methanol():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    # Properties of Methanol from Properties of Gases and Liquids 5th Ed. Appendix B.
    m.params.config.parameter_data = {
        "lennard_jones_sigma": 3.626 * pyunits.angstrom,
        "lennard_jones_epsilon_reduced": 481.8 * pyunits.K,
    }

    m.meta_object = PropertyClassMetadata()
    m.meta_object.default_units["temperature"] = pyunits.K
    m.meta_object.default_units["mass"] = pyunits.kg
    m.meta_object.default_units["length"] = pyunits.m
    m.meta_object.default_units["time"] = pyunits.s
    m.meta_object.default_units["amount"] = pyunits.mol
    m.meta_object.default_units["current"] = pyunits.ampere

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Create variables that should exist on param block
    m.params.mw = Var(initialize=32.042e-3, units=pyunits.kg / pyunits.mol)
    m.params.omega = Var(initialize=0.565, units=pyunits.dimensionless)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    # Poling et al. like to leave out the .15 from Kelvin conversion
    # and then give more sig figs than they ought
    m.props[1].temperature = Var(initialize=273 + 300, units=pyunits.K)

    ChapmanEnskogLennardJones.visc_d_phase_comp.build_parameters(m.params, "Vap")

    expr = ChapmanEnskogLennardJones.visc_d_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids
    T_list = [67, 127, 277]  # Gas temperature in C
    visc_list = [112, 132, 181]  # Experimental viscosities in millipoise
    for i in range(3):
        m.props[1].temperature.value = (T_list[i] + 273.15)
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)

@pytest.mark.unit
def test_visc_vap_comp_ethane():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = Block()

    m.params.config = ConfigBlock(implicit=True)
    # Properties of Ethane from Properties of Gases and Liquids 5th Ed. Appendix B.
    m.params.config.parameter_data = {
        "lennard_jones_sigma": 4.443 * pyunits.angstrom,
        "lennard_jones_epsilon_reduced": 215.7 * pyunits.K,
    }

    m.meta_object = PropertyClassMetadata()
    m.meta_object.default_units["temperature"] = pyunits.K
    m.meta_object.default_units["mass"] = pyunits.kg
    m.meta_object.default_units["length"] = pyunits.m
    m.meta_object.default_units["time"] = pyunits.s
    m.meta_object.default_units["amount"] = pyunits.mol
    m.meta_object.default_units["current"] = pyunits.ampere

    def get_metadata(self):
        return m.meta_object

    m.get_metadata = types.MethodType(get_metadata, m)
    m.params.get_metadata = types.MethodType(get_metadata, m.params)

    # Create variables that should exist on param block
    m.params.mw = Var(initialize=30.070e-3, units=pyunits.kg / pyunits.mol)
    m.params.omega = Var(initialize=0.099, units=pyunits.dimensionless)

    # Create a dummy state block
    m.props = Block([1])
    add_object_reference(m.props[1], "params", m.params)

    m.props[1].temperature = Var(initialize=273 + 300, units=pyunits.K)

    ChapmanEnskogLennardJones.visc_d_phase_comp.build_parameters(m.params, "Vap")

    expr = ChapmanEnskogLennardJones.visc_d_phase_comp.return_expression(
        m.props[1], m.params, "Vap", m.props[1].temperature
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids
    T_list = [47, 117, 247]  # Gas temperature in C
    visc_list = [100, 120, 156]  # Experimental viscosities in millipoise

    for i in range(3):
        m.props[1].temperature.value = (T_list[i] + 273.15)
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)
