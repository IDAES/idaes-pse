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
"""Regression tests for the supercritical boiler subflowsheet."""

import pytest

from pyomo.environ import value
from pyomo.network import Arc
from pyomo.opt import check_optimal_termination

from idaes.core import FlowsheetBlock
from idaes.core.scaling.util import jacobian_cond
from idaes.core.util.diagnostics_tools import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties import iapws95
from idaes.models.unit_models import Heater, Mixer
from idaes.models.unit_models.separator import Separator
from idaes.models.properties.general_helmholtz import helmholtz_available
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
)
from idaes.models_extra.power_generation.flowsheets.supercritical_power_plant import (
    boiler_subflowsheet_build as boiler,
)


pytestmark = [
    pytest.mark.skipif(
        not helmholtz_available(), reason="General Helmholtz not available"
    ),
    pytest.mark.integration,
    pytest.mark.solver,
]


@pytest.fixture(scope="module")
def solved_boiler():
    """Build, initialize, scale, and solve the flowsheet once."""
    model, solver = boiler.main()
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 100,
        "halt_on_ampl_error": "yes",
    }

    boiler.initialize(model)
    results = boiler.scale_solve(model)

    assert check_optimal_termination(results)
    assert degrees_of_freedom(model.fs) == 0
    return model


def _assert_baseline(expected):
    """Compare named model quantities with their accepted baseline values."""
    for name, (expression, baseline) in expected.items():
        assert value(expression) == pytest.approx(baseline, rel=1e-5), name


def test_structural_issues(solved_boiler):
    """Check for structural warnings in the fully specified flowsheet."""
    diagnostics = DiagnosticsToolbox(solved_boiler)
    diagnostics.report_structural_issues()
    diagnostics.assert_no_structural_warnings(
        ignore_unit_consistency=True,
        ignore_evaluation_errors=True,
    )


def test_build_flowsheet(solved_boiler):
    """Confirm the expected property packages, unit models, and arcs exist."""
    fs = solved_boiler.fs

    assert isinstance(fs, FlowsheetBlock)

    assert isinstance(fs.prop_water, iapws95.Iapws95ParameterBlock)
    assert isinstance(fs.prop_fluegas, FlueGasParameterBlock)

    for unit in (fs.ECON, fs.PrSH, fs.FSH, fs.RH):
        assert isinstance(unit, BoilerHeatExchanger)
    for unit in (fs.Water_wall, fs.PlSH):
        assert isinstance(unit, Heater)
    assert isinstance(fs.Spl1, Separator)
    for unit in (fs.mix1, fs.ATMP1):
        assert isinstance(unit, Mixer)

    for arc in (
        fs.econ2ww,
        fs.ww2prsh,
        fs.prsh2plsh,
        fs.plsh2fsh,
        fs.FSHtoATMP1,
        fs.fg_fsh2_separator,
        fs.fg_fsh2rh,
        fs.fg_fsh2PrSH,
        fs.fg_rhtomix,
        fs.fg_prsh2mix,
        fs.fg_mix2econ,
    ):
        assert isinstance(arc, Arc)


def test_primary_superheater_baseline(solved_boiler):
    """Check the accepted primary-superheater performance baseline."""
    unit = solved_boiler.fs.PrSH

    _assert_baseline(
        {
            "gas viscosity": (
                unit.hot_side.properties_in[0].visc_d,
                4.664336525463096e-5,
            ),
            "gas conductivity": (
                unit.hot_side.properties_in[0].therm_cond,
                0.07707457207753378,
            ),
            "tube velocity": (unit.v_tube[0], 3.59101431405311),
            "shell velocity": (unit.v_shell[0], 3.249018360523842),
            "tube Reynolds number": (unit.N_Re_tube[0], 788085.9423771807),
            "shell Reynolds number": (unit.N_Re_shell[0], 1350.9147570939767),
            "tube convection coefficient": (
                unit.hconv_tube[0],
                2816.879673687501,
            ),
            "shell radiation coefficient": (
                unit.hconv_shell_rad[0],
                12.458852026141935,
            ),
            "shell convection coefficient": (
                unit.hconv_shell_conv[0],
                21.875711082480247,
            ),
            "total shell coefficient": (
                unit.hconv_shell_total[0],
                34.334563108622184,
            ),
            "temperature driving force": (
                unit.delta_temperature[0],
                128.93723599571848,
            ),
            "inlet temperature difference": (
                unit.deltaT_1[0],
                403.6854949232379,
            ),
            "outlet temperature difference": (
                unit.deltaT_2[0],
                20.704375050837836,
            ),
            "tube pressure drop": (unit.deltaP_tube[0], -77110.48800420559),
            "shell pressure drop": (unit.deltaP_shell[0], -6.277316512429605),
            "mean beam length": (unit.mbl, 0.15544209377028503),
            "gas emissivity": (unit.gas_emissivity[0], 0.11925661199084645),
            "gas emissivity div2": (
                unit.gas_emissivity_div2[0],
                0.10043355541088558,
            ),
            "gas emissivity mul2": (
                unit.gas_emissivity_mul2[0],
                0.13934151172588213,
            ),
            "gas gray fraction": (
                unit.gas_gray_fraction[0],
                0.1639463768040806,
            ),
            "inlet liquid density": (
                unit.cold_side.properties_in[0].dens_mass_phase["Liq"],
                117.75352461273377,
            ),
            "outlet liquid density": (
                unit.cold_side.properties_out[0].dens_mass_phase["Liq"],
                98.60023285853784,
            ),
            "heat-transfer area": (unit.area, 13956.037360130911),
            "overall heat-transfer coefficient": (
                unit.overall_heat_transfer_coefficient[0],
                45.52948107169046,
            ),
        }
    )


def test_boiler_temperature_baseline(solved_boiler):
    """Check the accepted temperature profile through the boiler."""
    fs = solved_boiler.fs

    _assert_baseline(
        {
            "economizer water inlet": (
                fs.ECON.cold_side.properties_in[0].temperature,
                465.9419072298742,
            ),
            "economizer water outlet": (
                fs.ECON.cold_side.properties_out[0].temperature,
                546.880832860516,
            ),
            "economizer gas inlet": (
                fs.ECON.hot_side.properties_in[0].temperature,
                723.288504741046,
            ),
            "economizer gas outlet": (
                fs.ECON.hot_side.properties_out[0].temperature,
                497.9670014334912,
            ),
            "water-wall inlet": (
                fs.Water_wall.control_volume.properties_in[0].temperature,
                546.880832860516,
            ),
            "water-wall outlet": (
                fs.Water_wall.control_volume.properties_out[0].temperature,
                723.9095101569864,
            ),
            "primary superheater steam inlet": (
                fs.PrSH.cold_side.properties_in[0].temperature,
                723.9095101569864,
            ),
            "primary superheater steam outlet": (
                fs.PrSH.cold_side.properties_out[0].temperature,
                764.1677527834821,
            ),
            "primary superheater gas inlet": (
                fs.PrSH.hot_side.properties_in[0].temperature,
                1167.8532477067201,
            ),
            "primary superheater gas outlet": (
                fs.PrSH.hot_side.properties_out[0].temperature,
                744.6138852078241,
            ),
            "platen superheater steam inlet": (
                fs.PlSH.control_volume.properties_in[0].temperature,
                764.1677527834821,
            ),
            "platen superheater steam outlet": (
                fs.PlSH.control_volume.properties_out[0].temperature,
                797.4895466622573,
            ),
            "finishing superheater steam inlet": (
                fs.FSH.cold_side.properties_in[0].temperature,
                797.4895466622573,
            ),
            "finishing superheater steam outlet": (
                fs.FSH.cold_side.properties_out[0].temperature,
                873.3699075510588,
            ),
            "finishing superheater gas inlet": (
                fs.FSH.hot_side.properties_in[0].temperature,
                1300.335,
            ),
            "finishing superheater gas outlet": (
                fs.FSH.hot_side.properties_out[0].temperature,
                1167.8532477067201,
            ),
            "attemperator outlet": (
                fs.ATMP1.mixed_state[0].temperature,
                873.3698772408023,
            ),
            "reheater steam inlet": (
                fs.RH.cold_side.properties_in[0].temperature,
                581.224646872069,
            ),
            "reheater steam outlet": (
                fs.RH.cold_side.properties_out[0].temperature,
                871.6876469877258,
            ),
            "reheater gas inlet": (
                fs.RH.hot_side.properties_in[0].temperature,
                1167.8532477067201,
            ),
            "reheater gas outlet": (
                fs.RH.hot_side.properties_out[0].temperature,
                716.1530530112955,
            ),
        }
    )


def test_attemperator_enthalpy_baseline(solved_boiler):
    """Check the accepted attemperator inlet and outlet enthalpies."""
    fs = solved_boiler.fs

    assert value(fs.ATMP1.Steam.enth_mol[0]) == pytest.approx(
        62772.18366237012, rel=1e-5
    )
    assert value(fs.ATMP1.outlet.enth_mol[0]) == pytest.approx(
        62772.18202098774, rel=1e-5
    )


def test_flowsheet_jacobian_condition_number(solved_boiler):
    """Confirm the accepted scaled and unscaled Jacobian condition numbers."""
    unscaled = jacobian_cond(solved_boiler.fs, scaled=False)
    scaled = jacobian_cond(solved_boiler.fs, scaled=True)

    assert unscaled == pytest.approx(5.140885438597462e19, rel=1e-3)
    assert scaled == pytest.approx(169081.2224457273, rel=1e-3)
