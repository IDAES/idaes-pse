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
Tests for Wassiljew-Mason-Saxena mixing rules for thermal conductivity.
"""

import pytest
import types
from sys import modules

from pyomo.environ import ConcreteModel, Block, value, Var, Param, units as pyunits
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent
import pyomo.environ as pyo

from idaes.models.properties.modular_properties.transport_properties.thermal_conductivity_wms import (
    ThermalConductivityWMSPhase
)
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import wilke_phi_ij_callback
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

def construct_dummy_model(component_dict):
    m = ConcreteModel()
    components = [comp for comp in component_dict.keys()]
    comp1 = components[0]
    comp2 = components[1]

    m.params = DummyParameterBlock(
        components={
            comp1: {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "therm_cond_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (component_dict[comp1]["mw"], pyunits.g / pyunits.mol),
                }
            },
            comp2: {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "therm_cond_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (component_dict[comp2]["mw"], pyunits.g / pyunits.mol),
                }
            },
        },
        phases={
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Ideal,
                "therm_cond_phase": ThermalConductivityWMSPhase
            },
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.state_block_class([1], defined_state=False, parameters=m.params)
    ThermalConductivityWMSPhase.therm_cond_phase.build_parameters(m.params.Vap)

    # Add common variables
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )
    m.props[1].visc_d_phase_comp = Var(
        m.props[1].phase_list,
        m.props[1].component_list,
        initialize=0,
        units=pyunits.Pa * pyunits.s,
    )
    m.props[1].visc_d_phase_comp["Vap", comp1].value = pyunits.convert(
        component_dict[comp1]["visc_d_Vap"] * pyunits.micropoise,
        to_units=pyunits.Pa * pyunits.s,
    )
    m.props[1].visc_d_phase_comp["Vap", comp2].value = pyunits.convert(
        component_dict[comp2]["visc_d_Vap"] * pyunits.micropoise,
        to_units=pyunits.Pa * pyunits.s,
    )
    m.props[1].therm_cond_phase_comp = Var(
        m.props[1].phase_list,
        m.props[1].component_list,
        initialize=0,
        units=pyunits.W /pyunits.m /pyunits.K,
    )
    m.props[1].therm_cond_phase_comp["Vap", comp1].value = component_dict[comp1]["therm_cond_Vap"]
    m.props[1].therm_cond_phase_comp["Vap", comp2].value = component_dict[comp2]["therm_cond_Vap"]

    return m

# TODO: come back to this and compare calculated values to those contained in Figure 10-6
# Will need to calculate pure component viscosities through some correlation
# @pytest.mark.unit
# def test_wms_therm_cond_phase_():
#     # Taken from Example 10-5 of the properties of gases and liquids
#     component_dict = {
#         "benzene": {"mw": 78.114, "visc_d_Vap": 92.5,"therm_cond_Vap": 1.66e-2},
#         "Ar": {"mw": 39.948, "visc_d_Vap": 271, "therm_cond_Vap": 2.14e-2}
#     }
#     m = construct_dummy_model(component_dict)
#
#     assert m.params.Vap.viscosity_phi_ij_callback is wilke_phi_ij_callback
#
#     expr = ThermalConductivityWMSPhase.therm_cond_phase.return_expression(
#         m.props[1], m.params.Vap
#     )
#     m.props[1].mole_frac_phase_comp["Vap", "benzene"].value = 0.25
#     m.props[1].mole_frac_phase_comp["Vap", "Ar"].value = 0.75
#
#     phi_bb = m.props[1].visc_d_phi_ij["benzene", "benzene"]
#     assert_units_equivalent(phi_bb, pyunits.dimensionless)
#     assert pyo.value(phi_bb) == pytest.approx(1, rel=1e-12)
#
#     phi_aa = m.props[1].visc_d_phi_ij["Ar", "Ar"]
#     assert_units_equivalent(phi_aa, pyunits.dimensionless)
#     assert pyo.value(phi_aa) == pytest.approx(1, rel=1e-12)
#
#     phi_ba = m.props[1].visc_d_phi_ij["benzene", "Ar"]
#     assert_units_equivalent(phi_ba, pyunits.dimensionless)
#     assert pyo.value(phi_ba) == pytest.approx(0.459, rel=2e-3)
#
#     phi_ab = m.props[1].visc_d_phi_ij["Ar", "benzene"]
#     assert_units_equivalent(phi_ab, pyunits.dimensionless)
#     # Eq. 9-5.15, Properties of Gases and Liquids 5th Ed.
#     assert pyo.value(phi_ab) == pytest.approx(2.630, rel=2e-3)
#
#     assert pyo.value(expr) == pytest.approx(1.84e-2, rel=2e-3)
#
#     assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)
#
#
# @pytest.mark.unit
# def test_wilke_visc_d_phase_ammonia_hydrogen():
#     component_dict = {"NH3": {"mw": 17.031, "visc_d_Vap": 105.9}, "H2": {"mw": 2.016, "visc_d_Vap": 90.6}}
#     m = construct_dummy_model(component_dict)
#
#     assert m.params.Vap.viscosity_phi_ij_callback is wilke_phi_ij_callback
#
#     expr = ViscosityWilkePhase.visc_d.return_expression(
#         m.props[1], m.params.Vap
#     )
#
#     expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
#     # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
#     mole_frac_list = [0.0, 0.195, 0.399, 0.536, 0.677, 1.0]  # Mole fraction of NH3
#     visc_list = [90.6, 118.4, 123.8, 122.4, 120.0, 105.9]  # Experimental viscosities in millipoise
#     err_list = [0.0, -11.0, -12.0, -11.0, -9.7, 0.0]  # Percent error in Wilke's method
#     for i in range(6):
#         m.props[1].mole_frac_phase_comp["Vap", "NH3"] = mole_frac_list[i]
#         m.props[1].mole_frac_phase_comp["Vap", "H2"] = 1 - mole_frac_list[i]
#         err = 100 * value(expr_micropoise / visc_list[i] - 1)
#         assert err == pytest.approx(err_list[i], abs=0.3)
#
#     assert_units_equivalent(expr, pyunits.Pa * pyunits.s)
