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
Tests for Cubic equation of state methods

Author: Radhakrishna Tumbalam Gooty
"""

import pytest
from sys import modules
import os

from pyomo.environ import (check_optimal_termination,
                           ConcreteModel,
                           Constraint,
                           Expression,
                           Param,
                           Constraint,
                           value,
                           Var,
                           units as pyunits)

from idaes.core import (declare_process_block_class,
                        LiquidPhase, VaporPhase, SolidPhase)
from idaes.generic_models.properties.core.eos.ceos import Cubic, CubicType
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.core.util.constants import Constants as const
from idaes import bin_directory
from idaes.power_generation.properties.natural_gas_PR import get_prop
from idaes.core.util.misc import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom


# Reference: CoolProp: http://www.coolprop.org/index.html
data = {
            "phase": "Vap",
            "pressure": 101325,
            "temperature": 311.8730783132979,
            "enth_mol_phase": 24188.517506218643,
            "heat_capacity_ratio_phase": 1.290590129342867,
            "cp_mol_phase": 37.726122878482165,
            "cv_mol_phase": 29.2316840341025,
            "isentropic_speed_sound_phase": 279.2959607849875,
            "mole_frac_phase_comp": {
                "CO2": 0.94,
                "H2O": 0.06
            }
        }


@pytest.fixture()
def build_model():
    m = ConcreteModel()

    # Properties
    comp_props = get_prop(components=["CO2", "H2O"], phases=["Vap", "Liq"])

    # Parameters block
    m.params = GenericParameterBlock(default=comp_props)

    m.props = m.params.build_state_block(
        default={
            "defined_state": True,
            "parameters": m.params,
            "has_phase_equilibrium": True
        }
    )

    m.props.flow_mol.fix(100)
    m.props.pressure.fix(101325)
    m.props.mole_frac_comp["CO2"].fix(0.94)
    m.props.mole_frac_comp["H2O"].fix(0.06)

    m.props.temperature_constraint = Constraint(
        expr=m.props.temperature == m.props.temperature_dew["Vap", "Liq"]
    )

    assert degrees_of_freedom(m) == 0

    m.props.initialize(state_vars_fixed=True)
    results = get_solver(options={"bound_push": 1e-8}).solve(m)

    assert check_optimal_termination(results)

    return m


@pytest.mark.unit
def test_cp_mol_phase(build_model):
    m = build_model
    assert (
        str(pyunits.get_units(Cubic.cp_mol_phase(m.props, "Vap"))) ==
        'kg*m**2/K/mol/s**2'
    )

    assert (
        pytest.approx(value(Cubic.cp_mol_phase(m.props, "Vap")), rel=0.1) 
        == data["cp_mol_phase"] )


@pytest.mark.unit
def test_cv_mol_phase(build_model):
    m = build_model
    assert (
        str(pyunits.get_units(Cubic.cv_mol_phase(m.props, "Vap"))) ==
        'kg*m**2/K/mol/s**2'
    )

    assert (
        pytest.approx(value(Cubic.cv_mol_phase(m.props, "Vap")), rel=0.1) 
        == data["cv_mol_phase"] )


@pytest.mark.unit
def test_heat_capacity_ratio_phase(build_model):
    m = build_model
    assert (
        str(pyunits.get_units(Cubic.heat_capacity_ratio_phase(m.props, "Vap"))) ==
        'None'
    )

    assert (
        pytest.approx(value(Cubic.heat_capacity_ratio_phase(m.props, "Vap")), rel=0.1) 
        == data["heat_capacity_ratio_phase"] )


@pytest.mark.unit
def test_isentropic_speed_sound_phase(build_model):
    m = build_model
    assert (
        str(pyunits.get_units(Cubic.isentropic_speed_sound_phase(m.props, "Vap"))) ==
        'm/s'
    )

    assert (
        pytest.approx(value(Cubic.isentropic_speed_sound_phase(m.props, "Vap")), rel=0.1)
        == data["isentropic_speed_sound_phase"] )

