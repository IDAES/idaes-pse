#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for CubicComplementarityVLE formulation

Authors: Andrew Lee
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Param,
    value,
    Var,
    units as pyunits,
)

from idaes.core.base.phases import PhaseType
from idaes.core.util.exceptions import ConfigurationError

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.phase_equil.smooth_VLE_2 import (
    CubicComplementarityVLE,
    _calculate_temperature_slacks,
    _calculate_ceos_derivative_slacks,
    EPS_INIT,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity

from idaes.models.properties.modular_properties.pure.NIST import NIST


@pytest.mark.unit
def test_different_cubics():
    m = ConcreteModel()

    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {
                    "pressure_crit": (220.6e5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    "omega": 0.344,
                },
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        phases={
            "Liq": {
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
            "Vap": {
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.SRK},
            },
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        phases_in_equilibrium=[("Vap", "Liq")],
        phase_equilibrium_state={("Vap", "Liq"): CubicComplementarityVLE},
        parameter_data={
            "PR_kappa": {
                ("H2O", "H2O"): 0.000,
            },
            "SRK_kappa": {
                ("H2O", "H2O"): 0.000,
            },
        },
    )

    with pytest.raises(
        ConfigurationError,
        match="params - CubicComplementarityVLE formulation requires that both phases use the same "
        "type of cubic equation of state.",
    ):
        m.props = m.params.state_block_class([1], parameters=m.params)


@pytest.mark.unit
def test_non_cubic():
    m = ConcreteModel()

    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {
                    "pressure_crit": (220.6e5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    "omega": 0.344,
                },
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        phases={
            "Liq": {
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
            "Vap": {
                "equation_of_state": Ideal,
            },
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        phases_in_equilibrium=[("Vap", "Liq")],
        phase_equilibrium_state={("Vap", "Liq"): CubicComplementarityVLE},
        parameter_data={
            "PR_kappa": {
                ("H2O", "H2O"): 0.000,
            }
        },
    )

    with pytest.raises(
        ConfigurationError,
        match="params - CubicComplementarityVLE formulation only supports cubic equations of state.",
    ):
        m.props = m.params.state_block_class([1], parameters=m.params)


@pytest.fixture()
def frame():
    m = ConcreteModel()

    # Create a dummy parameter block
    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {
                    "pressure_crit": (220.6e5, pyunits.Pa),
                    "temperature_crit": (647, pyunits.K),
                    "omega": 0.344,
                },
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        phases={
            "Liq": {
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
            "Vap": {
                "equation_of_state": Cubic,
                "equation_of_state_options": {"type": CubicType.PR},
            },
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        phases_in_equilibrium=[("Vap", "Liq")],
        phase_equilibrium_state={("Vap", "Liq"): CubicComplementarityVLE},
        parameter_data={
            "PR_kappa": {
                ("H2O", "H2O"): 0.000,
            }
        },
    )

    # Create a dummy state block
    m.props = m.params.state_block_class([1], parameters=m.params)

    return m


@pytest.mark.unit
def test_build(frame):
    assert isinstance(frame.props[1].s_Vap_Liq, Var)
    assert isinstance(frame.props[1].gp_Vap_Liq, Var)
    assert isinstance(frame.props[1].gn_Vap_Liq, Var)

    assert isinstance(frame.props[1].cubic_second_derivative_Vap_Liq, Expression)

    assert isinstance(frame.props[1].eps_t_Vap_Liq, Param)
    assert value(frame.props[1].eps_t_Vap_Liq) == pytest.approx(1, rel=1e-8)
    assert isinstance(frame.props[1].eps_z_Vap_Liq, Param)
    assert value(frame.props[1].eps_z_Vap_Liq) == pytest.approx(1, rel=1e-8)

    assert isinstance(frame.props[1]._teq_constraint_Vap_Liq, Constraint)
    assert isinstance(
        frame.props[1].temperature_slack_complementarity_Vap_Liq, Constraint
    )
    assert isinstance(frame.props[1].cubic_slack_complementarity_Vap_Liq, Constraint)

    for k in ["Liq", "Vap"]:
        assert k in frame.props[1].s_Vap_Liq
        assert k in frame.props[1].gp_Vap_Liq
        assert k in frame.props[1].gn_Vap_Liq
        assert k in frame.props[1].cubic_second_derivative_Vap_Liq
        assert k in frame.props[1].temperature_slack_complementarity_Vap_Liq
        assert k in frame.props[1].cubic_slack_complementarity_Vap_Liq


# TODO: Need tests for formulation


@pytest.mark.unit
def test_cubic_second_derivative(frame):
    Z = frame.props[1].compress_fact_phase

    # For pure water
    R = 8.314
    b = 0.07780 * R * 647 / 220.6e5

    for P in range(1, 11):
        for T in range(300, 500, 10):
            frame.props[1].pressure.set_value(P * 1e5)
            frame.props[1].temperature.set_value(T)
            frame.props[1]._teq[("Vap", "Liq")].set_value(T)

            B = b * P * 1e5 / R / T

            for p in ["Vap", "Liq"]:
                der = value(6 * Z[p] + 2 * -(1 + B - 2 * B))
                assert value(
                    frame.props[1].cubic_second_derivative_Vap_Liq[p]
                ) == pytest.approx(der, rel=1e-8)


@pytest.mark.unit
def test_calculate_temperature_slacks(frame):
    s = frame.props[1].s_Vap_Liq
    # Teq > T
    frame.props[1].temperature.set_value(300)
    frame.props[1]._teq[("Vap", "Liq")].set_value(400)
    _calculate_temperature_slacks(frame.props[1], ("Vap", "Liq"), "Liq", "Vap")

    assert value(s["Vap"]) == pytest.approx(100, rel=1e-8)
    assert value(s["Liq"]) == pytest.approx(EPS_INIT, abs=1e-8)

    # Teq < T
    frame.props[1]._teq[("Vap", "Liq")].set_value(200)
    _calculate_temperature_slacks(frame.props[1], ("Vap", "Liq"), "Liq", "Vap")

    assert value(s["Liq"]) == pytest.approx(100, rel=1e-8)
    assert value(s["Vap"]) == pytest.approx(EPS_INIT, abs=1e-8)

    # Teq == T
    frame.props[1]._teq[("Vap", "Liq")].set_value(300)
    _calculate_temperature_slacks(frame.props[1], ("Vap", "Liq"), "Liq", "Vap")

    assert value(s["Liq"]) == pytest.approx(EPS_INIT, abs=1e-8)
    assert value(s["Vap"]) == pytest.approx(EPS_INIT, abs=1e-8)


@pytest.mark.unit
def test_calculate_ceos_derivative_slacks(frame):
    gp = frame.props[1].gp_Vap_Liq
    gn = frame.props[1].gn_Vap_Liq
    der = frame.props[1].cubic_second_derivative_Vap_Liq
    Z = frame.props[1].compress_fact_phase

    # For pure water
    R = 8.314
    b = 0.07780 * R * 647 / 220.6e5

    # Atmospheric conditions, vapor derivative +ve, liquid -ve
    frame.props[1].pressure.set_value(1e5)
    frame.props[1].temperature.set_value(300)
    frame.props[1]._teq[("Vap", "Liq")].set_value(300)

    _calculate_ceos_derivative_slacks(frame.props[1], ("Vap", "Liq"), "Liq", "Vap")

    B = value(b * frame.props[1].pressure / R / frame.props[1].temperature)
    der_l = value(6 * Z["Liq"] + 2 * -(1 + B - 2 * B))
    der_v = value(6 * Z["Vap"] + 2 * -(1 + B - 2 * B))

    assert value(gp["Liq"]) == pytest.approx(EPS_INIT, abs=1e-8)
    assert value(gn["Liq"]) == pytest.approx(-der_l, rel=1e-8)

    assert value(gp["Vap"]) == pytest.approx(der_v, rel=1e-8)
    assert value(gn["Vap"]) == pytest.approx(EPS_INIT, abs=1e-8)

    # Supercritical conditions, vapor derivative +ve, liquid +ve
    frame.props[1].pressure.set_value(250e5)
    frame.props[1].temperature.set_value(700)
    frame.props[1]._teq[("Vap", "Liq")].set_value(700)

    _calculate_ceos_derivative_slacks(frame.props[1], ("Vap", "Liq"), "Liq", "Vap")

    B = value(b * frame.props[1].pressure / R / frame.props[1].temperature)
    der_v = value(6 * Z["Vap"] + 2 * -(1 + B - 2 * B))

    assert value(gp["Liq"]) == pytest.approx(EPS_INIT, abs=1e-8)
    assert value(gn["Liq"]) == pytest.approx(EPS_INIT, abs=1e-8)

    assert value(gp["Vap"]) == pytest.approx(der_v, rel=1e-8)
    assert value(gn["Vap"]) == pytest.approx(EPS_INIT, abs=1e-8)

    # Supercritical conditions, vapor derivative -ve, liquid -ve
    frame.props[1].pressure.set_value(250e5)
    frame.props[1].temperature.set_value(300)
    frame.props[1]._teq[("Vap", "Liq")].set_value(300)

    _calculate_ceos_derivative_slacks(frame.props[1], ("Vap", "Liq"), "Liq", "Vap")

    B = value(b * frame.props[1].pressure / R / frame.props[1].temperature)
    der_l = value(6 * Z["Liq"] + 2 * -(1 + B - 2 * B))

    assert value(gp["Liq"]) == pytest.approx(EPS_INIT, abs=1e-8)
    assert value(gn["Liq"]) == pytest.approx(-der_l, rel=1e-8)

    assert value(gp["Vap"]) == pytest.approx(EPS_INIT, abs=1e-8)
    assert value(gn["Vap"]) == pytest.approx(EPS_INIT, abs=1e-8)


class TestWithNonCondensable:
    # TODO These tests could be fleshed out
    @pytest.fixture()
    def H2O_H2_model(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = GenericParameterBlock(
            components={
                "H2O": {
                    "parameter_data": {
                        "pressure_crit": (220.6e5, pyunits.Pa),
                        "temperature_crit": (647, pyunits.K),
                        "omega": 0.344,
                        "pressure_sat_comp_coeff": {  # NIST <- Stull 1947
                            "A": 4.6543,
                            "B": 1435.264,
                            "C": -64.848,
                        },
                    },
                    "pressure_sat_comp": NIST,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                },
                "H2": {
                    "valid_phase_types": [PhaseType.vaporPhase],
                    "parameter_data": {
                        "pressure_crit": (13e5, pyunits.Pa),
                        "temperature_crit": (33.2, pyunits.K),
                        "omega": -0.218,
                    },
                },
            },
            phases={
                "Liq": {
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
                "Vap": {
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
            },
            state_definition=FTPx,
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            phases_in_equilibrium=[("Vap", "Liq")],
            phase_equilibrium_state={("Vap", "Liq"): CubicComplementarityVLE},
            parameter_data={
                "PR_kappa": {
                    ("H2O", "H2O"): 0.000,
                    ("H2", "H2O"): 0.000,
                    ("H2O", "H2"): 0.000,
                    ("H2", "H2"): 0.000,
                }
            },
        )

        # Create a dummy state block
        m.props = m.params.state_block_class([1], parameters=m.params)

        return m

    @pytest.mark.unit
    def test_calculate_teq_permanent_gas(self, H2O_H2_model):
        m = H2O_H2_model
        m.props[1].pressure.set_value(1e5)
        m.props[1].temperature.set_value(300)
        CubicComplementarityVLE.calculate_teq(m.props[1], ("Vap", "Liq"))
        assert not m.props[1].is_property_constructed("temperature_dew")
        assert not m.props[1].is_property_constructed("temperature_bubble")


class TestWithNonVolatile:
    # TODO These tests could be fleshed out
    @pytest.fixture()
    def H2O_TEG_model(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = GenericParameterBlock(
            components={
                "H2O": {
                    "parameter_data": {
                        "pressure_crit": (220.6e5, pyunits.Pa),
                        "temperature_crit": (647, pyunits.K),
                        "omega": 0.344,
                        "pressure_sat_comp_coeff": {  # NIST <- Stull 1947
                            "A": 4.6543,
                            "B": 1435.264,
                            "C": -64.848,
                        },
                    },
                    "pressure_sat_comp": NIST,
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                },
                "TEG": {  # Triethylene Glycol
                    "valid_phase_types": [PhaseType.liquidPhase],
                    # Values from WolframAlpha
                    "parameter_data": {
                        "pressure_crit": (3.3e6, pyunits.Pa),
                        "temperature_crit": (797, pyunits.K),
                        "omega": 0.51,
                    },
                },
            },
            phases={
                "Liq": {
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
                "Vap": {
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
            },
            state_definition=FTPx,
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            phases_in_equilibrium=[("Vap", "Liq")],
            phase_equilibrium_state={("Vap", "Liq"): CubicComplementarityVLE},
            parameter_data={
                "PR_kappa": {
                    ("H2O", "H2O"): 0.000,
                    ("TEG", "H2O"): 0.000,
                    ("H2O", "TEG"): 0.000,
                    ("TEG", "TEG"): 0.000,
                }
            },
        )

        # Create a dummy state block
        m.props = m.params.state_block_class([1], parameters=m.params)

        return m

    @pytest.mark.unit
    def test_calculate_teq_permanent_gas(self, H2O_TEG_model):
        m = H2O_TEG_model
        m.props[1].pressure.set_value(1e5)
        m.props[1].temperature.set_value(300)
        CubicComplementarityVLE.calculate_teq(m.props[1], ("Vap", "Liq"))
        assert not m.props[1].is_property_constructed("temperature_dew")
        assert not m.props[1].is_property_constructed("temperature_bubble")
