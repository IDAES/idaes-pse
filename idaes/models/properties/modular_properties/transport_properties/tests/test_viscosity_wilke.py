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
Test Wilke mixing rules for viscosity from Properties of Gases and Liquids, 5th Ed.
Section 9-5-2

Authors: Douglas Allan
"""

import pytest
from sys import modules

from pyomo.environ import ConcreteModel, value, Var, units as pyunits
from pyomo.util.check_units import assert_units_equivalent
import pyomo.environ as pyo

from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import (
    ViscosityWilke,
    wilke_phi_ij_callback,
    herring_zimmer_phi_ij_callback,
)

from idaes.core import declare_process_block_class, LiquidPhase, VaporPhase

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import ConstantProperties


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


def construct_dummy_model(component_dict, callback):
    m = ConcreteModel()
    components = [comp for comp in component_dict.keys()]
    comp1 = components[0]
    comp2 = components[1]

    m.params = DummyParameterBlock(
        components={
            comp1: {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (component_dict[comp1]["mw"], pyunits.g / pyunits.mol),
                },
            },
            comp2: {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (component_dict[comp2]["mw"], pyunits.g / pyunits.mol),
                },
            },
        },
        phases={
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Ideal,
                "visc_d_phase": ViscosityWilke,
                "transport_property_options": {"viscosity_phi_ij_callback": callback},
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
    ViscosityWilke.visc_d_phase.build_parameters(m.params.Vap)

    # Add common variables
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )
    m.props[1]._visc_d_phase_comp = Var(
        m.props[1].phase_list,
        m.props[1].component_list,
        initialize=0,
        units=pyunits.Pa * pyunits.s,
    )
    m.props[1]._visc_d_phase_comp["Vap", comp1].value = pyunits.convert(
        component_dict[comp1]["visc_d_Vap"] * pyunits.micropoise,
        to_units=pyunits.Pa * pyunits.s,
    )
    m.props[1]._visc_d_phase_comp["Vap", comp2].value = pyunits.convert(
        component_dict[comp2]["visc_d_Vap"] * pyunits.micropoise,
        to_units=pyunits.Pa * pyunits.s,
    )
    return m


@pytest.mark.unit
def test_wilke_visc_d_phase_N2_CO2():
    component_dict = {
        "N2": {"mw": 28.014, "visc_d_Vap": 175.8},
        "CO2": {"mw": 44.009, "visc_d_Vap": 146.6},
    }
    m = construct_dummy_model(component_dict, wilke_phi_ij_callback)

    expr = ViscosityWilke.visc_d_phase.return_expression(m.props[1], "Vap")

    phi_N2_N2 = m.props[1].visc_d_phi_ij["N2", "N2"]
    assert_units_equivalent(phi_N2_N2, pyunits.dimensionless)
    assert pyo.value(phi_N2_N2) == pytest.approx(1, rel=1e-12)

    phi_CO2_CO2 = m.props[1].visc_d_phi_ij["CO2", "CO2"]
    assert_units_equivalent(phi_CO2_CO2, pyunits.dimensionless)
    assert pyo.value(phi_CO2_CO2) == pytest.approx(1, rel=1e-12)

    phi_N2_CO2 = m.props[1].visc_d_phi_ij["N2", "CO2"]
    assert_units_equivalent(phi_N2_CO2, pyunits.dimensionless)
    assert pyo.value(phi_N2_CO2) == pytest.approx(1.369409683998896, rel=1e-12)

    phi_CO2_N2 = m.props[1].visc_d_phi_ij["CO2", "N2"]
    assert_units_equivalent(phi_CO2_N2, pyunits.dimensionless)
    # Eq. 9-5.15, Properties of Gases and Liquids 5th Ed.
    assert pyo.value(phi_CO2_N2) == pytest.approx(
        pyo.value(phi_N2_CO2 * 146.6 / 175.8 * 28.014 / 44.009), rel=1e-12
    )

    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.213, 0.495, 0.767, 1.0]  # Mole fraction of N2
    visc_list = [
        146.6,
        153.5,
        161.8,
        172.1,
        175.8,
    ]  # Experimental viscosities in millipoise
    err_list = [0.0, -1.3, -1.8, -2.8, 0.0]  # Percent error in Wilke's method
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "N2"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "CO2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(err_list[i], abs=0.3)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_wilke_visc_d_phase_ammonia_hydrogen():
    component_dict = {
        "NH3": {"mw": 17.031, "visc_d_Vap": 105.9},
        "H2": {"mw": 2.016, "visc_d_Vap": 90.6},
    }
    m = construct_dummy_model(component_dict, wilke_phi_ij_callback)

    expr = ViscosityWilke.visc_d_phase.return_expression(m.props[1], "Vap")

    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.195, 0.399, 0.536, 0.677, 1.0]  # Mole fraction of NH3
    visc_list = [
        90.6,
        118.4,
        123.8,
        122.4,
        120.0,
        105.9,
    ]  # Experimental viscosities in millipoise
    err_list = [0.0, -11.0, -12.0, -11.0, -9.7, 0.0]  # Percent error in Wilke's method
    for i in range(6):
        m.props[1].mole_frac_phase_comp["Vap", "NH3"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "H2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(err_list[i], abs=0.3)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
@pytest.mark.unit
def test_wilke_visc_d_phase_N2O_SO2():
    component_dict = {
        "N2O": {"mw": 44.013, "visc_d_Vap": 173.0},
        "SO2": {"mw": 64.066, "visc_d_Vap": 152.3},
    }
    m = construct_dummy_model(component_dict, wilke_phi_ij_callback)

    expr = ViscosityWilke.visc_d_phase.return_expression(m.props[1], "Vap")

    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.325, 0.625, 0.817, 1.0]  # Mole fraction of NH3
    visc_list = [
        152.3,
        161.7,
        167.8,
        170.7,
        173.0,
    ]  # Experimental viscosities in millipoise
    err_list = [0.0, -2.2, -2.2, -1.3, 0.0]  # Percent error in Wilke's method
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "N2O"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "SO2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(err_list[i], abs=0.3)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_herning_zipperer_visc_d_phase_N2_CO2():
    component_dict = {
        "N2": {"mw": 28.014, "visc_d_Vap": 175.8},
        "CO2": {"mw": 44.009, "visc_d_Vap": 146.6},
    }
    m = construct_dummy_model(component_dict, herring_zimmer_phi_ij_callback)

    expr = ViscosityWilke.visc_d_phase.return_expression(m.props[1], "Vap")

    phi_N2_N2 = m.props[1].visc_d_phi_ij["N2", "N2"]
    assert_units_equivalent(phi_N2_N2, pyunits.dimensionless)
    assert pyo.value(phi_N2_N2) == pytest.approx(1, rel=1e-12)

    phi_CO2_CO2 = m.props[1].visc_d_phi_ij["CO2", "CO2"]
    assert_units_equivalent(phi_CO2_CO2, pyunits.dimensionless)
    assert pyo.value(phi_CO2_CO2) == pytest.approx(1, rel=1e-12)

    phi_N2_CO2 = m.props[1].visc_d_phi_ij["N2", "CO2"]
    assert_units_equivalent(phi_N2_CO2, pyunits.dimensionless)
    assert pyo.value(phi_N2_CO2) == pytest.approx(1.253381233999109, rel=1e-12)

    phi_CO2_N2 = m.props[1].visc_d_phi_ij["CO2", "N2"]
    assert_units_equivalent(phi_CO2_N2, pyunits.dimensionless)
    # Eq. 9-5.15, Properties of Gases and Liquids 5th Ed.
    assert pyo.value(phi_CO2_N2) == pytest.approx(1 / pyo.value(phi_N2_CO2), rel=1e-12)

    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.213, 0.495, 0.767, 1.0]  # Mole fraction of N2
    visc_list = [
        146.6,
        153.5,
        161.8,
        172.1,
        175.8,
    ]  # Experimental viscosities in millipoise
    err_list = [0.0, -1, -1.5, -2.5, 0.0]  # Percent error in H&Z's method
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "N2"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "CO2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(err_list[i], abs=0.3)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
def test_herning_zipperer_visc_d_phase_ammonia_hydrogen():
    component_dict = {
        "NH3": {"mw": 17.031, "visc_d_Vap": 105.9},
        "H2": {"mw": 2.016, "visc_d_Vap": 90.6},
    }
    m = construct_dummy_model(component_dict, herring_zimmer_phi_ij_callback)

    expr = ViscosityWilke.visc_d_phase.return_expression(m.props[1], "Vap")

    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.195, 0.399, 0.536, 0.677, 1.0]  # Mole fraction of NH3
    visc_list = [
        90.6,
        118.4,
        123.8,
        122.4,
        120.0,
        105.9,
    ]  # Experimental viscosities in millipoise
    err_list = [0.0, -18.0, -19.0, -16.0, -14.0, 0.0]  # Percent error in H&Z's method
    for i in range(6):
        m.props[1].mole_frac_phase_comp["Vap", "NH3"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "H2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(
            err_list[i], abs=1
        )  # Not as many sig figs in the table

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)


@pytest.mark.unit
@pytest.mark.unit
def test_herning_zipperer_visc_d_phase_N2O_SO2():
    component_dict = {
        "N2O": {"mw": 44.013, "visc_d_Vap": 173.0},
        "SO2": {"mw": 64.066, "visc_d_Vap": 152.3},
    }
    m = construct_dummy_model(component_dict, herring_zimmer_phi_ij_callback)

    expr = ViscosityWilke.visc_d_phase.return_expression(m.props[1], "Vap")

    expr_micropoise = pyunits.convert(expr, pyunits.micropoise)
    # Pulled from Table 9.2, Properties of Gases and Liquids 5th Ed.
    mole_frac_list = [0.0, 0.325, 0.625, 0.817, 1.0]  # Mole fraction of NH3
    visc_list = [
        152.3,
        161.7,
        167.8,
        170.7,
        173.0,
    ]  # Experimental viscosities in millipoise
    err_list = [0.0, -2.2, -2.1, -1.2, 0.0]  # Percent error in H&Z's method
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "N2O"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "SO2"] = 1 - mole_frac_list[i]
        err = 100 * value(expr_micropoise / visc_list[i] - 1)
        assert err == pytest.approx(err_list[i], abs=0.3)

    assert_units_equivalent(expr, pyunits.Pa * pyunits.s)
