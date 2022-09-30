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
from sys import modules

from pyomo.environ import ConcreteModel, Block, value, Var, Param, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent
import pyomo.environ as pyo

from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import (
    ViscosityWilkePhase, wilke_phi_ij_callback
)

from idaes.core.util.misc import add_object_reference
from idaes.core.base.property_meta import PropertyClassMetadata
from idaes.core import declare_process_block_class, LiquidPhase, VaporPhase

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import ConstantProperties
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True

    def parameters(self):
        self.parameters_set = True


def define_state(b):
    b.state_defined = True

# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass

@pytest.mark.unit
def test_visc_phase_N2_CO2():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(
        components={
            "N2": {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (28.014, pyunits.g/pyunits.mol),
                    "visc_d_Vap_comp_coeff": (175.8, pyunits.micropoise)
                }
            },
            "CO2": {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (44.009, pyunits.g/pyunits.mol),
                    "visc_d_Vap_comp_coeff": (146.6, pyunits.micropoise)
                }
            },
        },
        phases={
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Ideal,
                "visc_d_phase": ViscosityWilkePhase
            },
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
            "current": pyunits.ampere,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.state_block_class([1], defined_state=False, parameters=m.params)
    ViscosityWilkePhase.visc_d.build_parameters(m.params.Vap)

    # Add common variables
    m.props[1].temperature = Var(initialize=300, units=pyunits.K)
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )
    m.props[1].visc_d_phase_comp = Var(m.props[1].phase_list, m.props[1].component_list, initialize=0)
    m.props[1].visc_d_phase_comp["Vap", "N2"].value = pyunits.convert(
        175.8 * pyunits.micropoise,
        to_units=pyunits.Pa*pyunits.s
    )
    m.props[1].visc_d_phase_comp["Vap", "CO2"].value = pyunits.convert(
        146.6 * pyunits.micropoise,
        to_units=pyunits.Pa*pyunits.s
    )

    assert m.params.Vap.viscosity_phi_ij_callback is wilke_phi_ij_callback

    expr = ViscosityWilkePhase.visc_d.return_expression(
        m.props[1], m.params.Vap
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)

    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.213, 0.495, 0.767, 1.0]  # Mole fraction of N2
    visc_list = [146.6, 153.5, 161.8, 175.8]   # Experimental viscosities in millipoise
    err_list = [0.0, -1.3, -1.8, -2.8]  # Percent error in Wilke's method
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "N2"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "CO2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(err_list[i], abs=0.3)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_visc_mol_comp_methanol():
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

    ChapmanEnskogLennardJones.viscosity_dynamic_vap_comp.build_parameters(m.params)

    expr = ChapmanEnskogLennardJones.viscosity_dynamic_vap_comp.return_expression(
        m.props[1], m.params
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
def test_visc_mol_comp_ethane():
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

    ChapmanEnskogLennardJones.viscosity_dynamic_vap_comp.build_parameters(m.params)

    expr = ChapmanEnskogLennardJones.viscosity_dynamic_vap_comp.return_expression(
        m.props[1], m.params
    )
    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids
    T_list = [47, 117, 247]  # Gas temperature in C
    visc_list = [100, 120, 156]  # Experimental viscosities in millipoise

    for i in range(3):
        m.props[1].temperature.value = (T_list[i] + 273.15)
        assert pyo.value(expr_micropoise) == pytest.approx(visc_list[i], rel=0.03)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)
