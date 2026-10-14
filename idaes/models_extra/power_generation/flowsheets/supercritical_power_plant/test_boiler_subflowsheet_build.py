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
from idaes.core.util.testing import assert_solution_equivalent
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


expected_results = {
    # External water/steam boundary conditions
    "ECON.cold_side.properties_in[0].flow_mol": {
        None: (24194.177, 1e-8, None),
    },
    "ECON.cold_side.properties_in[0].enth_mol": {
        None: (14990.97, 1e-8, None),
    },
    "ECON.cold_side.properties_in[0].pressure": {
        None: (26922222.222, 1e-8, None),
    },
    "RH.cold_side.properties_in[0].flow_mol": {
        None: (21235.27, 1e-8, None),
    },
    "RH.cold_side.properties_in[0].enth_mol": {
        None: (53942.7569, 1e-8, None),
    },
    "RH.cold_side.properties_in[0].pressure": {
        None: (3677172.33638, 1e-8, None),
    },
    # External flue-gas boundary conditions
    "FSH.hot_side.properties_in[0].flow_mol_comp": {
        "H2O": (21290.6999 * 0.0869, 1e-8, None),
        "CO2": (21290.6999 * 0.1449, 1e-8, None),
        "N2": (21290.6999 * 0.7434, 1e-8, None),
        "O2": (21290.6999 * 0.0247, 1e-8, None),
        "NO": (21290.6999 * 0.0006, 1e-8, None),
        "SO2": (21290.6999 * 0.002, 1e-8, None),
    },
    "FSH.hot_side.properties_in[0].pressure": {
        None: (100145, 1e-8, None),
    },
    # Fixed duties, spray-water flow, and gas split
    "Water_wall.heat_duty": {0: (7.51e8, 1e-8, None)},
    "PlSH.heat_duty": {0: (5.5e7, 1e-8, None)},
    "ATMP1.SprayWater_state[0].flow_mol": {
        None: (0.001, 1e-8, None),
    },
    "ATMP1.SprayWater_state[0].pressure": {
        None: (1.22e8, 1e-8, None),
    },
    "Spl1.split_fraction": {(0, "outlet_1"): (0.75, 1e-8, None)},
    "PrSH.hot_side.properties_in[0].visc_d": {
        None: (4.664336525463096e-5, 1e-5, None),
    },
    "PrSH.hot_side.properties_in[0].therm_cond": {
        None: (0.07707457207753378, 1e-5, None),
    },
    "PrSH.v_tube": {0: (3.59101431405311, 1e-5, None)},
    "PrSH.v_shell": {0: (3.249018360523842, 1e-5, None)},
    "PrSH.N_Re_tube": {0: (788085.9423771807, 1e-5, None)},
    "PrSH.N_Re_shell": {0: (1350.9147570939767, 1e-5, None)},
    "PrSH.hconv_tube": {0: (2816.879673687501, 1e-5, None)},
    "PrSH.hconv_shell_rad": {0: (12.458852026141935, 1e-5, None)},
    "PrSH.hconv_shell_conv": {0: (21.875711082480247, 1e-5, None)},
    "PrSH.hconv_shell_total": {0: (34.334563108622184, 1e-5, None)},
    "PrSH.delta_temperature": {0: (128.93723599571848, 1e-5, None)},
    "PrSH.deltaT_1": {0: (403.6854949232379, 1e-5, None)},
    "PrSH.deltaT_2": {0: (20.704375050837836, 1e-5, None)},
    "PrSH.deltaP_tube": {0: (-77110.48800420559, 1e-5, None)},
    "PrSH.deltaP_shell": {0: (-6.277316512429605, 1e-5, None)},
    "PrSH.mbl": {None: (0.15544209377028503, 1e-5, None)},
    "PrSH.gas_emissivity": {0: (0.11925661199084645, 1e-5, None)},
    "PrSH.gas_emissivity_div2": {0: (0.10043355541088558, 1e-5, None)},
    "PrSH.gas_emissivity_mul2": {0: (0.13934151172588213, 1e-5, None)},
    "PrSH.gas_gray_fraction": {0: (0.1639463768040806, 1e-5, None)},
    "PrSH.cold_side.properties_in[0].dens_mass_phase": {
        "Liq": (117.75352461273377, 1e-5, None),
    },
    "PrSH.cold_side.properties_out[0].dens_mass_phase": {
        "Liq": (98.60023285853784, 1e-5, None),
    },
    "PrSH.area": {None: (13956.037360130911, 1e-5, None)},
    "PrSH.overall_heat_transfer_coefficient": {
        0: (45.52948107169046, 1e-5, None),
    },
    "ECON.cold_side.properties_in[0].temperature": {
        None: (465.9419072298742, 1e-5, None),
    },
    "ECON.cold_side.properties_out[0].temperature": {
        None: (546.880832860516, 1e-5, None),
    },
    "ECON.hot_side.properties_in[0].temperature": {
        None: (723.288504741046, 1e-5, None),
    },
    "ECON.hot_side.properties_out[0].temperature": {
        None: (497.9670014334912, 1e-5, None),
    },
    "Water_wall.control_volume.properties_in[0].temperature": {
        None: (546.880832860516, 1e-5, None),
    },
    "Water_wall.control_volume.properties_out[0].temperature": {
        None: (723.9095101569864, 1e-5, None),
    },
    "PrSH.cold_side.properties_in[0].temperature": {
        None: (723.9095101569864, 1e-5, None),
    },
    "PrSH.cold_side.properties_out[0].temperature": {
        None: (764.1677527834821, 1e-5, None),
    },
    "PrSH.hot_side.properties_in[0].temperature": {
        None: (1167.8532477067201, 1e-5, None),
    },
    "PrSH.hot_side.properties_out[0].temperature": {
        None: (744.6138852078241, 1e-5, None),
    },
    "PlSH.control_volume.properties_in[0].temperature": {
        None: (764.1677527834821, 1e-5, None),
    },
    "PlSH.control_volume.properties_out[0].temperature": {
        None: (797.4895466622573, 1e-5, None),
    },
    "FSH.cold_side.properties_in[0].temperature": {
        None: (797.4895466622573, 1e-5, None),
    },
    "FSH.cold_side.properties_out[0].temperature": {
        None: (873.3699075510588, 1e-5, None),
    },
    "FSH.hot_side.properties_in[0].temperature": {
        None: (1300.335, 1e-5, None),
    },
    "FSH.hot_side.properties_out[0].temperature": {
        None: (1167.8532477067201, 1e-5, None),
    },
    "ATMP1.mixed_state[0].temperature": {
        None: (873.3698772408023, 1e-5, None),
    },
    "ATMP1.Steam.enth_mol": {0: (62772.18366237012, 1e-5, None)},
    "ATMP1.outlet.enth_mol": {0: (62772.18202098774, 1e-5, None)},
    "RH.cold_side.properties_in[0].temperature": {
        None: (581.224646872069, 1e-5, None),
    },
    "RH.cold_side.properties_out[0].temperature": {
        None: (871.6876469877258, 1e-5, None),
    },
    "RH.hot_side.properties_in[0].temperature": {
        None: (1167.8532477067201, 1e-5, None),
    },
    "RH.hot_side.properties_out[0].temperature": {
        None: (716.1530530112955, 1e-5, None),
    },
}


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


def test_structural_issues(solved_boiler):
    """Check for structural warnings in the fully specified flowsheet."""
    diagnostics = DiagnosticsToolbox(solved_boiler)
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


def test_solution(solved_boiler):
    """Check the KPI values produced by the original boiler flowsheet."""
    assert_solution_equivalent(solved_boiler.fs, expected_results)


def test_water_steam_flow_conservation(solved_boiler):
    """Check molar-flow conservation along both water/steam paths."""
    fs = solved_boiler.fs
    tol = 1e-7

    for unit in (fs.ECON, fs.PrSH, fs.FSH, fs.RH):
        assert value(unit.cold_side.properties_out[0].flow_mol) == pytest.approx(
            value(unit.cold_side.properties_in[0].flow_mol), rel=tol
        )

    for unit in (fs.Water_wall, fs.PlSH):
        assert value(unit.control_volume.properties_out[0].flow_mol) == pytest.approx(
            value(unit.control_volume.properties_in[0].flow_mol), rel=tol
        )

    assert value(fs.ATMP1.outlet.flow_mol[0]) == pytest.approx(
        value(fs.ATMP1.Steam.flow_mol[0]) + value(fs.ATMP1.SprayWater.flow_mol[0]),
        rel=tol,
    )


def test_flue_gas_component_conservation(solved_boiler):
    """Check component-flow conservation through the complete flue-gas path."""
    fs = solved_boiler.fs
    components = tuple(fs.prop_fluegas.component_list)
    tol = 1e-7

    for component in components:
        # Non-reacting heat exchangers conserve each gas component.
        for unit in (fs.FSH, fs.RH, fs.PrSH, fs.ECON):
            assert value(
                unit.hot_side.properties_out[0].flow_mol_comp[component]
            ) == pytest.approx(
                value(unit.hot_side.properties_in[0].flow_mol_comp[component]),
                rel=tol,
                abs=1e-10,
            )

        # The splitter branches sum to its inlet.
        assert value(fs.Spl1.inlet.flow_mol_comp[0, component]) == pytest.approx(
            value(fs.Spl1.outlet_1.flow_mol_comp[0, component])
            + value(fs.Spl1.outlet_2.flow_mol_comp[0, component]),
            rel=tol,
            abs=1e-10,
        )

        # The mixer outlet equals the sum of its two branch inlets.
        assert value(fs.mix1.outlet.flow_mol_comp[0, component]) == pytest.approx(
            value(fs.mix1.Reheat_out.flow_mol_comp[0, component])
            + value(fs.mix1.PrSH_out.flow_mol_comp[0, component]),
            rel=tol,
            abs=1e-10,
        )


def test_heat_exchanger_energy_conservation(solved_boiler):
    """Check equal and opposite heat transfer on both sides of each exchanger."""
    fs = solved_boiler.fs

    for unit in (fs.ECON, fs.PrSH, fs.FSH, fs.RH):
        cold_side_heat = value(unit.cold_side.heat[0])
        hot_side_heat = value(unit.hot_side.heat[0])

        assert cold_side_heat + hot_side_heat == pytest.approx(
            0, abs=max(abs(cold_side_heat), abs(hot_side_heat)) * 1e-7
        )


def test_pressure_and_splitter_kpis(solved_boiler):
    """Check pressure-loss direction and the solved flue-gas split."""
    fs = solved_boiler.fs

    for unit in (fs.ECON, fs.PrSH, fs.FSH, fs.RH):
        assert value(unit.cold_side.properties_out[0].pressure) < value(
            unit.cold_side.properties_in[0].pressure
        )
        assert value(unit.hot_side.properties_out[0].pressure) < value(
            unit.hot_side.properties_in[0].pressure
        )

    inlet_flow = sum(
        value(fs.Spl1.inlet.flow_mol_comp[0, component])
        for component in fs.prop_fluegas.component_list
    )
    outlet_1_flow = sum(
        value(fs.Spl1.outlet_1.flow_mol_comp[0, component])
        for component in fs.prop_fluegas.component_list
    )
    outlet_2_flow = sum(
        value(fs.Spl1.outlet_2.flow_mol_comp[0, component])
        for component in fs.prop_fluegas.component_list
    )

    assert outlet_1_flow / inlet_flow == pytest.approx(0.75, rel=1e-7)
    assert outlet_2_flow / inlet_flow == pytest.approx(0.25, rel=1e-7)


def test_flowsheet_jacobian_condition_number(solved_boiler):
    """Confirm the accepted scaled and unscaled Jacobian condition numbers."""
    unscaled = jacobian_cond(solved_boiler.fs, scaled=False)
    scaled = jacobian_cond(solved_boiler.fs, scaled=True)

    assert unscaled == pytest.approx(5.140885438597462e19, rel=1e-3)
    assert scaled == pytest.approx(169081.2224457273, rel=1e-3)
