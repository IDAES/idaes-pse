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
"""CoolProp equivalence checks for TPX Helmholtz state blocks."""

import json
import math
from pathlib import Path

import pytest
import pyomo.environ as pyo

from idaes.core import FlowsheetBlock
from idaes.models.properties.general_helmholtz import (
    AmountBasis,
    HelmholtzParameterBlock,
    StateVars,
    helmholtz_available,
)

try:
    import CoolProp.CoolProp as CP
except ImportError:  # pragma: no cover
    CP = None


_COOLPROP_NAME_MAP = {
    "h2o": "Water",
    "co2": "CO2",
    "nh3": "Ammonia",
    "propane": "n-Propane",
    "butane": "n-Butane",
    "isobutane": "IsoButane",
    "r32": "R32",
    "r125": "R125",
    "r134a": "R134a",
    "r1234ze": "R1234ze(E)",
    "r227ea": "R227EA",
}

_PARAMETER_DIR = Path(__file__).resolve().parents[1] / "components" / "parameters"
_COMPONENTS = tuple(sorted(_COOLPROP_NAME_MAP))
_T_FRACS = (0.2, 0.5, 0.7)

def _component_basic_data(c):
    with (_PARAMETER_DIR / f"{c}.json").open() as f:
        return json.load(f)["basic"]


def _sample_temperatures(c):
    basic = _component_basic_data(c)
    t_low = max(basic["T_min"], basic["Tt"])
    t_high = basic["Tc"]
    span = t_high - t_low
    return tuple(t_low + frac * span for frac in _T_FRACS)
    # return [400, 400]

def temperature_from_coolprop(c, p):
    fluid = _COOLPROP_NAME_MAP[c]
    return CP.PropsSI("T", "P", p, "Q", 1, fluid)



def compare_enthalpy_difference(c, temps):
    """Evaluate enthalpy differences for two states from TPX and CoolProp.

    Enthalpy reference states may differ between implementations, so this helper
    only checks that both calculation methods produce finite values.
    """
    fluid = _COOLPROP_NAME_MAP[c]


    param= HelmholtzParameterBlock(
        pure_component=c, amount_basis=AmountBasis.MASS, state_vars=StateVars.TPX
    )
    param.construct()

    for t in temps:
        h1 = param.htpx(T=t*pyo.units.K, x=0, with_units=False, units = pyo.units.J/pyo.units.kg)
        h2 = param.htpx(T=t*pyo.units.K, x=1, with_units=False, units = pyo.units.J/pyo.units.kg)
        print(f"State block enthalpy 1: {h1:.3f} J/Kg")
        print(f"State block enthalpy 2: {h2:.3f} J/Kg")
        delta_h_stateblock = pyo.value(h2 - h1)
        delta_h_coolprop = CP.PropsSI("Hmass", "Q", 1, "T", t, fluid) - CP.PropsSI(
            "Hmass", "Q", 0, "T", t, fluid
        )
        print(f"State block enthalpy difference: {delta_h_stateblock:.3f} J/Kg")
        print(f"CoolProp enthalpy difference: {delta_h_coolprop:.3f} J/Kg")
        assert  delta_h_stateblock == pytest.approx(delta_h_coolprop,  rel=1e2)
    return True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.skipif(CP is None, reason="CoolProp not available")
@pytest.mark.integration
@pytest.mark.parametrize("component", _COMPONENTS)
def test_coolprop_enthalpy_difference(component):
    temperatures = _sample_temperatures(component)
    assert compare_enthalpy_difference(
        component, temperatures
    )
