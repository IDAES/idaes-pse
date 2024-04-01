#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pytest
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
import idaes.core.util.scaling as iscale
from idaes.models.unit_models import HeatExchangerFlowPattern
from idaes.models.properties.modular_properties import GenericParameterBlock
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)
from idaes.models_extra.power_generation.unit_models import CrossFlowHeatExchanger1D
import idaes.core.util.model_statistics as mstat
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

# Set up solver
optarg = {
    # 'bound_push' : 1e-6,
    "constr_viol_tol": 1e-8,
    "nlp_scaling_method": "user-scaling",
    "linear_solver": "ma57",
    "OF_ma57_automatic_scaling": "yes",
    "max_iter": 350,
    "tol": 1e-8,
    "halt_on_ampl_error": "no",
}
solver = get_solver("ipopt", options=optarg)


def _create_model(pressure_drop):
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.h2_side_prop_params = GenericParameterBlock(
        **get_prop(["H2", "H2O", "Ar", "N2"], {"Vap"}, eos=EosType.IDEAL),
        doc="H2O + H2 gas property parameters",
    )
    m.fs.heat_exchanger = CrossFlowHeatExchanger1D(
        has_holdup=True,
        dynamic=False,
        cold_side={
            "property_package": m.fs.h2_side_prop_params,
            "has_holdup": False,
            "dynamic": False,
            "has_pressure_change": pressure_drop,
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "FORWARD",
        },
        hot_side={
            "property_package": m.fs.h2_side_prop_params,
            "has_holdup": False,
            "dynamic": False,
            "has_pressure_change": pressure_drop,
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "BACKWARD",
        },
        shell_is_hot=True,
        flow_type=HeatExchangerFlowPattern.countercurrent,
        finite_elements=12,
        tube_arrangement="staggered",
    )

    hx = m.fs.heat_exchanger

    hx.hot_side_inlet.flow_mol.fix(2619.7)
    hx.hot_side_inlet.temperature.fix(971.6)
    hx.hot_side_inlet.pressure.fix(1.2e5)
    hx.hot_side_inlet.mole_frac_comp[0, "H2"].fix(0.79715)
    hx.hot_side_inlet.mole_frac_comp[0, "H2O"].fix(0.20177)
    hx.hot_side_inlet.mole_frac_comp[0, "Ar"].fix(0.00086358)
    hx.hot_side_inlet.mole_frac_comp[0, "N2"].fix(0.00021589)

    hx.cold_side_inlet.flow_mol.fix(2619.7)
    hx.cold_side_inlet.temperature.fix(446.21)
    hx.cold_side_inlet.pressure.fix(1.2e5)
    hx.cold_side_inlet.mole_frac_comp[0, "H2"].fix(0.36203)
    hx.cold_side_inlet.mole_frac_comp[0, "H2O"].fix(0.63689)
    hx.cold_side_inlet.mole_frac_comp[0, "Ar"].fix(0.00086358)
    hx.cold_side_inlet.mole_frac_comp[0, "N2"].fix(0.00021589)

    hx.di_tube.fix(0.0525018)
    hx.thickness_tube.fix(0.0039116)
    hx.length_tube_seg.fix(4.3)
    hx.nseg_tube.fix(12)
    hx.ncol_tube.fix(50)
    hx.nrow_inlet.fix(25)

    hx.pitch_x.fix(0.1)
    hx.pitch_y.fix(0.1)
    hx.therm_cond_wall = 43.0
    hx.rfouling_tube = 0.0001
    hx.rfouling_shell = 0.0001
    hx.fcorrection_htc_tube.fix(1)
    hx.fcorrection_htc_shell.fix(1)
    if pressure_drop:
        hx.fcorrection_dp_tube.fix(1)
        hx.fcorrection_dp_shell.fix(1)

    hx.cp_wall.value = 502.4

    pp = m.fs.h2_side_prop_params
    pp.set_default_scaling("enth_mol_phase", 1e-3)
    pp.set_default_scaling("pressure", 1e-5)
    pp.set_default_scaling("temperature", 1e-2)
    pp.set_default_scaling("flow_mol", 1e-3)

    _mf_scale = {
        "H2": 1,
        "H2O": 1,
        "N2": 10,
        "Ar": 10,
    }
    for comp, s in _mf_scale.items():
        pp.set_default_scaling("mole_frac_comp", s, index=comp)
        pp.set_default_scaling("mole_frac_phase_comp", s, index=("Vap", comp))
        pp.set_default_scaling("flow_mol_phase_comp", s * 1e-3, index=("Vap", comp))

    shell = hx.hot_side
    tube = hx.cold_side
    iscale.set_scaling_factor(shell.area, 1e-1)
    # ssf(hx.shell.heat, 1e-6)
    iscale.set_scaling_factor(tube.area, 1)
    # ssf(hx.tube.heat, 1e-6)
    iscale.set_scaling_factor(shell._enthalpy_flow, 1e-8)  # pylint: disable=W0212
    iscale.set_scaling_factor(tube._enthalpy_flow, 1e-8)  # pylint: disable=W0212
    iscale.set_scaling_factor(shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(hx.heat_holdup, 1e-8)

    iscale.calculate_scaling_factors(m)

    return m


def _check_model_statistics(m, deltaP):
    fixed_unused_var_set = {
        "fs.h2_side_prop_params.H2.omega",
        "fs.h2_side_prop_params.H2.pressure_crit",
        "fs.h2_side_prop_params.H2.temperature_crit",
        "fs.h2_side_prop_params.H2.cp_mol_ig_comp_coeff_G",
        "fs.h2_side_prop_params.H2.cp_mol_ig_comp_coeff_H",
        "fs.h2_side_prop_params.H2O.omega",
        "fs.h2_side_prop_params.H2O.pressure_crit",
        "fs.h2_side_prop_params.H2O.temperature_crit",
        "fs.h2_side_prop_params.H2O.cp_mol_ig_comp_coeff_G",
        "fs.h2_side_prop_params.H2O.cp_mol_ig_comp_coeff_H",
        "fs.h2_side_prop_params.H2O.pressure_sat_comp_coeff_A",
        "fs.h2_side_prop_params.H2O.pressure_sat_comp_coeff_B",
        "fs.h2_side_prop_params.H2O.pressure_sat_comp_coeff_C",
        "fs.h2_side_prop_params.Ar.omega",
        "fs.h2_side_prop_params.Ar.pressure_crit",
        "fs.h2_side_prop_params.Ar.temperature_crit",
        "fs.h2_side_prop_params.Ar.cp_mol_ig_comp_coeff_G",
        "fs.h2_side_prop_params.Ar.cp_mol_ig_comp_coeff_H",
        "fs.h2_side_prop_params.N2.omega",
        "fs.h2_side_prop_params.N2.pressure_crit",
        "fs.h2_side_prop_params.N2.temperature_crit",
        "fs.h2_side_prop_params.N2.cp_mol_ig_comp_coeff_G",
        "fs.h2_side_prop_params.N2.cp_mol_ig_comp_coeff_H",
    }
    if not deltaP:
        fixed_unused_var_set.add("fs.heat_exchanger.delta_elevation")

    for var in mstat.fixed_unused_variables_set(m):
        assert var.name in fixed_unused_var_set

    unfixed_unused_var_set = {
        "fs.heat_exchanger.hot_side.material_flow_dx[0.0,0.0,Vap,H2]",
        "fs.heat_exchanger.hot_side.material_flow_dx[0.0,0.0,Vap,H2O]",
        "fs.heat_exchanger.hot_side.material_flow_dx[0.0,0.0,Vap,Ar]",
        "fs.heat_exchanger.hot_side.material_flow_dx[0.0,0.0,Vap,N2]",
        "fs.heat_exchanger.hot_side.enthalpy_flow_dx[0.0,0.0,Vap]",
        "fs.heat_exchanger.hot_side.pressure_dx[0.0,0.0]",
        "fs.heat_exchanger.cold_side.material_flow_dx[0.0,1.0,Vap,H2]",
        "fs.heat_exchanger.cold_side.material_flow_dx[0.0,1.0,Vap,H2O]",
        "fs.heat_exchanger.cold_side.material_flow_dx[0.0,1.0,Vap,Ar]",
        "fs.heat_exchanger.cold_side.material_flow_dx[0.0,1.0,Vap,N2]",
        "fs.heat_exchanger.cold_side.enthalpy_flow_dx[0.0,1.0,Vap]",
        "fs.heat_exchanger.cold_side.pressure_dx[0.0,1.0]",
    }

    for var in mstat.unused_variables_set(m) - mstat.fixed_unused_variables_set(m):
        assert var.name in unfixed_unused_var_set

    assert len(mstat.deactivated_constraints_set(m)) == 0


@pytest.fixture
def model_no_dP():
    m = _create_model(pressure_drop=False)
    return m


@pytest.mark.component
def test_initialization(model_no_dP):
    m = model_no_dP

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=False)

    m.fs.heat_exchanger.initialize_build(optarg=optarg)

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=False)
    assert pyo.value(
        m.fs.heat_exchanger.hot_side_outlet.temperature[0]
    ) == pytest.approx(485.34, abs=1e-1)
    assert pyo.value(
        m.fs.heat_exchanger.cold_side_outlet.temperature[0]
    ) == pytest.approx(911.47, abs=1e-1)


@pytest.mark.integration
def test_units(model_no_dP):
    assert_units_consistent(model_no_dP.fs.heat_exchanger)


@pytest.fixture
def model_dP():
    m = _create_model(pressure_drop=True)
    return m


@pytest.mark.component
def test_initialization_dP(model_dP):
    m = model_dP

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=True)

    m.fs.heat_exchanger.initialize_build(optarg=optarg)

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=True)

    assert pyo.value(
        m.fs.heat_exchanger.hot_side_outlet.temperature[0]
    ) == pytest.approx(485.34, abs=1e-1)
    assert pyo.value(
        m.fs.heat_exchanger.cold_side_outlet.temperature[0]
    ) == pytest.approx(911.47, abs=1e-1)
    assert pyo.value(m.fs.heat_exchanger.hot_side_outlet.pressure[0]) == pytest.approx(
        118870.08569, abs=1
    )
    assert pyo.value(m.fs.heat_exchanger.cold_side_outlet.pressure[0]) == pytest.approx(
        111418.71399, abs=1
    )


@pytest.mark.integration
def test_units_dP(model_dP):
    assert_units_consistent(model_dP.fs.heat_exchanger)
