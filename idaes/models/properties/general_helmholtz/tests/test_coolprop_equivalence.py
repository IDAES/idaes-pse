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
from idaes.models.properties.general_helmholtz.helmholtz_functions import (
    HelmholtzThermoExpressions,
)

try:
    import CoolProp.CoolProp as CP
except ImportError:  # pragma: no cover
    CP = None


_COOLPROP_NAME_MAP = {
    "1-butene": "1-Butene",
    "acetone": "Acetone",
    "nh3": "Ammonia",
    "argon": "Argon",
    "benzene": "Benzene",
    "co2": "CO2",
    "carbonmonoxide": "CarbonMonoxide",
    "carbonylsulfide": "CarbonylSulfide",
    "cyclohexane": "CycloHexane",
    "cyclopentane": "Cyclopentane",
    "d4": "D4",
    "d5": "D5",
    "deuterium": "Deuterium",
    "dichloroethane": "Dichloroethane",
    "diethylether": "DiethylEther",
    "dimethylcarbonate": "DimethylCarbonate",
    "dimethylether": "DimethylEther",
    "ethane": "Ethane",
    "ethanol": "Ethanol",
    "ethylbenzene": "EthylBenzene",
    "ethylene": "Ethylene",
    "ethyleneoxide": "EthyleneOxide",
    "h2o": "Water",
    "heavywater": "HeavyWater",
    "helium": "Helium",
    "hydrogen": "Hydrogen",
    "hydrogenchloride": "HydrogenChloride",
    "hydrogensulfide": "HydrogenSulfide",
    "isobutane": "IsoButane",
    "isobutene": "IsoButene",
    "isohexane": "Isohexane",
    "isopentane": "Isopentane",
    "krypton": "Krypton",
    "md2m": "MD2M",
    "md3m": "MD3M",
    "md4m": "MD4M",
    "mdm": "MDM",
    "mm": "MM",
    "methane": "Methane",
    "neon": "Neon",
    "neopentane": "Neopentane",
    "nitrogen": "Nitrogen",
    "nitrousoxide": "NitrousOxide",
    "novec649": "Novec649",
    "orthodeuterium": "OrthoDeuterium",
    "orthohydrogen": "OrthoHydrogen",
    "oxygen": "Oxygen",
    "paradeuterium": "ParaDeuterium",
    "parahydrogen": "ParaHydrogen",
    "propylene": "Propylene",
    "r113": "R113",
    "r115": "R115",
    "r116": "R116",
    "r12": "R12",
    "r1233zd(e)": "R1233zd(E)",
    "r1234yf": "R1234yf",
    "r1234ze": "R1234ze(E)",
    "r1234ze(z)": "R1234ze(Z)",
    "r124": "R124",
    "r1243zf": "R1243zf",
    "r125": "R125",
    "r1336mzz(e)": "R1336mzz(E)",
    "r134a": "R134a",
    "r13i1": "R13I1",
    "r141b": "R141b",
    "r142b": "R142b",
    "r152a": "R152A",
    "r161": "R161",
    "r218": "R218",
    "r227ea": "R227EA",
    "r23": "R23",
    "r236ea": "R236EA",
    "r236fa": "R236FA",
    "r245ca": "R245ca",
    "r245fa": "R245fa",
    "r32": "R32",
    "r365mfc": "R365MFC",
    "r40": "R40",
    "r404a": "R404A",
    "r41": "R41",
    "r410a": "R410A",
    "r507a": "R507A",
    "sulfurdioxide": "SulfurDioxide",
    "sulfurhexafluoride": "SulfurHexafluoride",
    "toluene": "Toluene",
    "xenon": "Xenon",
    "cis-2-butene": "cis-2-Butene",
    "m-xylene": "m-Xylene",
    "butane": "n-Butane",
    "n-decane": "n-Decane",
    "n-dodecane": "n-Dodecane",
    "n-hexane": "n-Hexane",
    "n-nonane": "n-Nonane",
    "n-octane": "n-Octane",
    "n-pentane": "n-Pentane",
    "propane": "n-Propane",
    "o-xylene": "o-Xylene",
    "p-xylene": "p-Xylene",
    "trans-2-butene": "trans-2-Butene",
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


def temperature_from_coolprop(c, p):
    fluid = _COOLPROP_NAME_MAP[c]
    return CP.PropsSI("T", "P", p, "Q", 1, fluid)


def compare_enthalpy_difference(c, temps):
    """Evaluate enthalpy differences for two states from TPX and CoolProp.

    Enthalpy reference states may differ between implementations, so this test
    evaluates differences irrespective of reference offset.
    """
    fluid = _COOLPROP_NAME_MAP[c]

    param = HelmholtzParameterBlock(
        pure_component=c, amount_basis=AmountBasis.MASS, state_vars=StateVars.TPX
    )
    param.construct()

    for t in temps:
        h1 = param.htpx(
            T=t * pyo.units.K, x=0, with_units=False, units=pyo.units.J / pyo.units.kg
        )
        h2 = param.htpx(
            T=t * pyo.units.K, x=1, with_units=False, units=pyo.units.J / pyo.units.kg
        )
        print(f"State block enthalpy 1: {h1:.3f} J/Kg")
        print(f"State block enthalpy 2: {h2:.3f} J/Kg")
        delta_h_stateblock = pyo.value(h2 - h1)
        delta_h_coolprop = CP.PropsSI("Hmass", "Q", 1, "T", t, fluid) - CP.PropsSI(
            "Hmass", "Q", 0, "T", t, fluid
        )
        print(f"State block enthalpy difference: {delta_h_stateblock:.3f} J/Kg")
        print(f"CoolProp enthalpy difference: {delta_h_coolprop:.3f} J/Kg")
        assert delta_h_stateblock == pytest.approx(delta_h_coolprop, rel=1e2)
    return True


def compare_entropy_difference(c, temps):
    """Evaluate entropy differences for two states from TPX and CoolProp.

    Entropy reference states may differ between implementations, so this test
    evaluates differences irrespective of reference offset.
    """
    fluid = _COOLPROP_NAME_MAP[c]

    param = HelmholtzParameterBlock(
        pure_component=c, amount_basis=AmountBasis.MASS, state_vars=StateVars.TPX
    )
    param.construct()

    for t in temps:
        s1 = param.stpx(
            T=t * pyo.units.K,
            x=0,
            with_units=False,
            units=pyo.units.J / pyo.units.K / pyo.units.kg,
        )
        s2 = param.stpx(
            T=t * pyo.units.K,
            x=1,
            with_units=False,
            units=pyo.units.J / pyo.units.K / pyo.units.kg,
        )
        print(f"State block entropy 1: {s1:.3f} J/k/Kg")
        print(f"State block entropy 2: {s2:.3f} J/k/Kg")
        delta_s_stateblock = pyo.value(s2 - s1)
        delta_s_coolprop = CP.PropsSI("Smass", "Q", 1, "T", t, fluid) - CP.PropsSI(
            "Smass", "Q", 0, "T", t, fluid
        )
        print(f"State block entropy difference: {delta_s_stateblock:.3f} J/k/Kg")
        print(f"CoolProp entropy difference: {delta_s_coolprop:.3f} J/k/Kg")
        assert delta_s_stateblock == pytest.approx(delta_s_coolprop, rel=1e2)
    return True


def compare_mass_density(c, temps):
    """Evaluate mass density from TPX and CoolProp.
    Finds saturation pressure for each temperature, then evaluates density at shifted pressures for the vapor and liquid phases.
    """
    fluid = _COOLPROP_NAME_MAP[c]
    m = pyo.ConcreteModel()

    m.param = HelmholtzParameterBlock(
        pure_component=c, amount_basis=AmountBasis.MASS, state_vars=StateVars.TPX
    )
    te = HelmholtzThermoExpressions(m, m.param)

    for t in temps:
        p_sat_l = pyo.value(te.p(T=t * pyo.units.K, x=0) * pyo.units.Pa)
        p_sat_v = pyo.value(te.p(T=t * pyo.units.K, x=1) * pyo.units.Pa)
        p_val = p_sat_l * 1.1
        rho_l = pyo.value(te.rho_liq(T=t * pyo.units.K, p=p_val * pyo.units.Pa))
        rho_l_cp = CP.PropsSI("Dmass", "T", t, "P", p_val, fluid)
        print(f"IDAES rho_liq for T {t} P {p_val}: {rho_l}")
        print(f"CoolProp rho_liq for T {t} P {p_val}: {rho_l_cp}")
        assert rho_l == pytest.approx(rho_l_cp, rel=1e-1)
        for frac in [0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]:
            p_val = p_sat_v * frac
            rho_v = pyo.value(
                te.rho_vap(T=t * pyo.units.kelvin, p=p_val * pyo.units.Pa)
            )
            rho_v_cp = CP.PropsSI("Dmass", "T", t, "P", p_val, fluid)
            print(f"IDAES rho_vap for T {t} P {p_val}: {rho_v}")
            print(f"CoolProp rho_vap for T {t} P {p_val}: {rho_v_cp}")
            assert rho_v == pytest.approx(rho_v_cp, rel=1e-1)
    return True


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.skipif(CP is None, reason="CoolProp not available")
@pytest.mark.integration
@pytest.mark.parametrize("component", _COMPONENTS)
def test_coolprop_comparisons(component):
    temperatures = _sample_temperatures(component)
    assert compare_enthalpy_difference(component, temperatures)
    assert compare_entropy_difference(component, temperatures)
    assert compare_mass_density(component, temperatures)
