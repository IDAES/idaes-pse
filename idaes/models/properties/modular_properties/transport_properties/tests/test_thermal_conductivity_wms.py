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
from sys import modules

from pyomo.environ import ConcreteModel, value, Var, units as pyunits
from pyomo.util.check_units import assert_units_equivalent
import pyomo.environ as pyo

from idaes.models.properties.modular_properties.transport_properties import (
    ThermalConductivityWMS,
)
from idaes.models.properties.modular_properties.transport_properties.viscosity_wilke import (
    wilke_phi_ij_callback,
)
from idaes.core import declare_process_block_class, LiquidPhase, VaporPhase

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import (
    ConstantProperties,
    ChapmanEnskog,
)


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


def construct_dummy_model(component_dict, chapman_enskog=False):
    m = ConcreteModel()
    components = [comp for comp in component_dict.keys()]
    comp1 = components[0]
    comp2 = components[1]

    config_dict = {
        comp1: {
            "therm_cond_phase_comp": {"Vap": ConstantProperties.Constant},
            "parameter_data": {
                "mw": (component_dict[comp1]["mw"], pyunits.g / pyunits.mol),
                "therm_cond_Vap_comp_coeff": (
                    component_dict[comp1]["therm_cond_Vap"],
                    pyunits.W / pyunits.m / pyunits.K,
                ),
            },
        },
        comp2: {
            "therm_cond_phase_comp": {"Vap": ConstantProperties.Constant},
            "parameter_data": {
                "mw": (component_dict[comp2]["mw"], pyunits.g / pyunits.mol),
                "therm_cond_Vap_comp_coeff": (
                    component_dict[comp2]["therm_cond_Vap"],
                    pyunits.W / pyunits.m / pyunits.K,
                ),
            },
        },
    }
    if chapman_enskog:
        config_dict[comp1]["visc_d_phase_comp"] = {
            "Vap": ChapmanEnskog.ChapmanEnskogLennardJones
        }
        config_dict[comp2]["visc_d_phase_comp"] = {
            "Vap": ChapmanEnskog.ChapmanEnskogLennardJones
        }
        config_dict[comp1]["parameter_data"]["lennard_jones_sigma"] = (
            component_dict[comp1]["sigma"],
            pyunits.angstrom,
        )
        config_dict[comp2]["parameter_data"]["lennard_jones_sigma"] = (
            component_dict[comp2]["sigma"],
            pyunits.angstrom,
        )
        config_dict[comp1]["parameter_data"]["lennard_jones_epsilon_reduced"] = (
            component_dict[comp1]["epsilon"],
            pyunits.K,
        )
        config_dict[comp2]["parameter_data"]["lennard_jones_epsilon_reduced"] = (
            component_dict[comp2]["epsilon"],
            pyunits.K,
        )
    else:
        config_dict[comp1]["visc_d_phase_comp"] = {"Vap": ConstantProperties.Constant}
        config_dict[comp2]["visc_d_phase_comp"] = {"Vap": ConstantProperties.Constant}
        config_dict[comp1]["parameter_data"]["visc_d_Vap_comp_coeff"] = (
            component_dict[comp1]["visc_d_Vap"],
            pyunits.micropoise,
        )
        config_dict[comp2]["parameter_data"]["visc_d_Vap_comp_coeff"] = (
            component_dict[comp2]["visc_d_Vap"],
            pyunits.micropoise,
        )

    m.params = DummyParameterBlock(
        components=config_dict,
        phases={
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Ideal,
                "therm_cond_phase": ThermalConductivityWMS,
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
    ThermalConductivityWMS.therm_cond_phase.build_parameters(m.params.Vap)
    # Add common variables
    m.props[1].temperature = Var(initialize=300, units=pyunits.K)
    m.props[1].mole_frac_phase_comp = Var(
        m.params.phase_list, m.params.component_list, initialize=0.5
    )

    def rule_visc_d_phase_comp(b, p, j):
        cobj = b.params.get_component(j)
        if (
            cobj.config.visc_d_phase_comp is not None
            and p in cobj.config.visc_d_phase_comp
        ):
            return cobj.config.visc_d_phase_comp[p].visc_d_phase_comp.return_expression(
                b, cobj, p, b.temperature
            )
        else:
            return pyo.Expression.Skip

    m.props[1]._visc_d_phase_comp = pyo.Expression(
        m.props[1].phase_list,
        m.props[1].component_list,
        doc="Viscosity for each phase-component pair",
        rule=rule_visc_d_phase_comp,
    )

    m.props[1]._therm_cond_phase_comp = Var(
        m.props[1].phase_list,
        m.props[1].component_list,
        initialize=0,
        units=pyunits.W / pyunits.m / pyunits.K,
    )
    m.props[1]._therm_cond_phase_comp["Vap", comp1].value = component_dict[comp1][
        "therm_cond_Vap"
    ]
    m.props[1]._therm_cond_phase_comp["Vap", comp2].value = component_dict[comp2][
        "therm_cond_Vap"
    ]

    return m


@pytest.mark.unit
def test_wms_therm_cond_phase_():
    # Taken from Example 10-5 of the properties of gases and liquids
    component_dict = {
        "benzene": {"mw": 78.114, "visc_d_Vap": 92.5, "therm_cond_Vap": 1.66e-2},
        "Ar": {"mw": 39.948, "visc_d_Vap": 271, "therm_cond_Vap": 2.14e-2},
    }
    m = construct_dummy_model(component_dict, chapman_enskog=False)

    assert (
        m.params.Vap.config.transport_property_options["viscosity_phi_ij_callback"]
        is wilke_phi_ij_callback
    )

    expr = ThermalConductivityWMS.therm_cond_phase.return_expression(m.props[1], "Vap")
    m.props[1].mole_frac_phase_comp["Vap", "benzene"].value = 0.25
    m.props[1].mole_frac_phase_comp["Vap", "Ar"].value = 0.75

    phi_bb = m.props[1].visc_d_phi_ij["benzene", "benzene"]
    assert_units_equivalent(phi_bb, pyunits.dimensionless)
    assert pyo.value(phi_bb) == pytest.approx(1, rel=1e-12)

    phi_aa = m.props[1].visc_d_phi_ij["Ar", "Ar"]
    assert_units_equivalent(phi_aa, pyunits.dimensionless)
    assert pyo.value(phi_aa) == pytest.approx(1, rel=1e-12)

    phi_ba = m.props[1].visc_d_phi_ij["benzene", "Ar"]
    assert_units_equivalent(phi_ba, pyunits.dimensionless)
    assert pyo.value(phi_ba) == pytest.approx(0.459, rel=2e-3)

    phi_ab = m.props[1].visc_d_phi_ij["Ar", "benzene"]
    assert_units_equivalent(phi_ab, pyunits.dimensionless)
    # Eq. 9-5.15, Properties of Gases and Liquids 5th Ed.
    assert pyo.value(phi_ab) == pytest.approx(2.630, rel=2e-3)

    assert pyo.value(expr) == pytest.approx(1.84e-2, rel=2e-3)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)


@pytest.mark.unit
def test_wms_therm_cond_phase_benzene_n_hexane():
    # Unfortunately, in order to test the thermal conductivity mixing rule, we need to know component
    # viscosities. In order to estimate the visocisity at the particular experimental conditions
    # we use Chapman Enskog with viscosity Lennard Jones parameters from
    # Properties of Gases and Liquids, 5th Ed., Appendix B
    component_dict = {
        "benzene": {
            "mw": 78.114,
            "sigma": 5.349,
            "epsilon": 412.3,
            "therm_cond_Vap": 1.90093e-2,
        },
        "n-hexane": {
            "mw": 86.178,
            "sigma": 5.949,
            "epsilon": 399.3,
            "therm_cond_Vap": 2.29209e-2,
        },
    }
    m = construct_dummy_model(component_dict, chapman_enskog=True)

    assert (
        m.params.Vap.config.transport_property_options["viscosity_phi_ij_callback"]
        is wilke_phi_ij_callback
    )

    expr = ThermalConductivityWMS.therm_cond_phase.return_expression(m.props[1], "Vap")
    # Pulled off Figure 10-6 from Properties of Gases and Liquids, 5th Ed.
    # That table cites: Bennett, L. A., and R. G. Vines: J. Chem. Phys., 23: 1587 (1955).
    mole_frac_list = [0.0, 0.25188, 0.50090, 0.74994, 1.0]  # Mole fraction of benzene
    therm_cond_list = [
        2.29209e-2,
        2.186596e-2,
        2.087267e-2,
        1.997177e-2,
        1.900933e-2,
    ]  # Experimental tcond in W/(m K)
    # The temperature is mislabeled in Figure 10-6. This is the actual temperature from Bennett and Vines (1955)
    m.props[1].temperature.value = 125.0 + 273.15
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "benzene"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "n-hexane"] = 1 - mole_frac_list[i]
        assert value(expr) == pytest.approx(therm_cond_list[i], rel=0.01)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)


@pytest.mark.unit
def test_wms_therm_cond_phase_methanol_n_hexane():
    # Unfortunately, in order to test the thermal conductivity mixing rule, we need to know component
    # viscosities. In order to estimate the visocisity at the particular experimental conditions
    # we use Chapman Enskog with viscosity Lennard Jones parameters from
    # Properties of Gases and Liquids, 5th Ed., Appendix B
    component_dict = {
        "methanol": {
            "mw": 32.04,
            "sigma": 3.626,
            "epsilon": 481.8,
            "therm_cond_Vap": 2.178133e-2,
        },
        "n-hexane": {
            "mw": 86.178,
            "sigma": 5.949,
            "epsilon": 399.3,
            "therm_cond_Vap": 1.98716e-2,
        },
    }
    m = construct_dummy_model(component_dict, chapman_enskog=True)

    assert (
        m.params.Vap.config.transport_property_options["viscosity_phi_ij_callback"]
        is wilke_phi_ij_callback
    )

    expr = ThermalConductivityWMS.therm_cond_phase.return_expression(m.props[1], "Vap")
    # Pulled off Figure 10-6 from Properties of Gases and Liquids, 5th Ed.
    # That table cites: Bennett, L. A., and R. G. Vines:J. Chem. Phys., 23: 1587 (1955).
    mole_frac_list = [
        0.0,
        0.25021,
        0.501179,
        0.748935,
        1.0,
    ]  # Mole fraction of methanol
    therm_cond_list = [
        1.98716e-2,
        2.12654e-2,
        2.22279e-2,
        2.22279e-2,
        2.17813e-2,
    ]  # Experimental tcond in W/(m K)
    m.props[1].temperature.value = 98.4 + 273.15
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "methanol"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "n-hexane"] = 1 - mole_frac_list[i]
        # This mixing rule doesn't appear to be particularly reliable for polar-nonpolar solutions
        assert value(expr) == pytest.approx(therm_cond_list[i], rel=0.08)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)


@pytest.mark.unit
def test_wms_therm_cond_phase_benzene_argon():
    # Unfortunately, in order to test the thermal conductivity mixing rule, we need to know component
    # viscosities. In order to estimate the visocisity at the particular experimental conditions
    # we use Chapman Enskog with viscosity Lennard Jones parameters from
    # Properties of Gases and Liquids, 5th Ed., Appendix B
    component_dict = {
        "benzene": {
            "mw": 78.114,
            "sigma": 5.349,
            "epsilon": 412.3,
            "therm_cond_Vap": 1.663763e-2,
        },
        "Ar": {
            "mw": 39.95,
            "sigma": 3.542,
            "epsilon": 93.3,
            "therm_cond_Vap": 2.13809e-2,
        },
    }
    m = construct_dummy_model(component_dict, chapman_enskog=True)

    assert (
        m.params.Vap.config.transport_property_options["viscosity_phi_ij_callback"]
        is wilke_phi_ij_callback
    )

    expr = ThermalConductivityWMS.therm_cond_phase.return_expression(m.props[1], "Vap")
    # Pulled off Figure 10-6 from Properties of Gases and Liquids, 5th Ed.
    # That table cites: Bennett, L. A., and R. G. Vines:J. Chem. Phys., 23: 1587 (1955).
    mole_frac_list = [0.0, 0.25130, 0.50026, 0.74780, 1.0]  # Mole fraction of benzene
    therm_cond_list = [
        2.13809e-2,
        1.907851e-2,
        1.776181e-2,
        1.709196e-2,
        1.663763e-2,
    ]  # Experimental tcond in W/(m K)
    m.props[1].temperature.value = 100.6 + 273.15
    for i in range(5):
        m.props[1].mole_frac_phase_comp["Vap", "benzene"] = mole_frac_list[i]
        m.props[1].mole_frac_phase_comp["Vap", "Ar"] = 1 - mole_frac_list[i]
        assert value(expr) == pytest.approx(therm_cond_list[i], rel=0.05)

    assert_units_equivalent(expr, pyunits.W / pyunits.m / pyunits.K)
