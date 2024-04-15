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
from idaes.models.properties.modular_properties import GenericParameterBlock
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)
from idaes.models_extra.power_generation.unit_models import Heater1D, Heater1DInitializer
import idaes.core.util.model_statistics as mstat
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

# Set up solver
optarg = {
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
    m.fs.heater = Heater1D(
        property_package=m.fs.h2_side_prop_params,
        has_holdup=True,
        dynamic=False,
        has_fluid_holdup=False,
        has_pressure_change=pressure_drop,
        finite_elements=4,
        tube_arrangement="in-line",
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    heater = m.fs.heater

    heater.inlet.flow_mol.fix(5102.5)
    heater.inlet.temperature.fix(938.83)
    heater.inlet.pressure.fix(1.2e5)
    heater.inlet.mole_frac_comp[0, "H2"].fix(0.57375)
    heater.inlet.mole_frac_comp[0, "H2O"].fix(0.42517)
    heater.inlet.mole_frac_comp[0, "Ar"].fix(0.00086358)
    heater.inlet.mole_frac_comp[0, "N2"].fix(0.00021589)

    heater.di_tube.fix(0.0525018)
    heater.thickness_tube.fix(0.0039116)
    heater.pitch_x.fix(0.1)
    heater.pitch_y.fix(0.1)
    heater.length_tube_seg.fix(10)
    heater.nseg_tube.fix(1)
    heater.rfouling = 0.0001
    heater.fcorrection_htc_shell.fix(1)
    heater.cp_wall = 502.4
    if pressure_drop:
        heater.fcorrection_dp_shell.fix(1)

    heater.ncol_tube.fix(40)
    heater.nrow_inlet.fix(40)
    heater.electric_heat_duty.fix(3.6504e06)

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

    shell = heater.control_volume
    iscale.set_scaling_factor(shell.area, 1e-1)
    iscale.set_scaling_factor(shell.heat, 1e-6)
    iscale.set_scaling_factor(shell._enthalpy_flow, 1e-8)  # pylint: disable=W0212
    iscale.set_scaling_factor(shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(heater.heat_holdup, 1e-8)

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

    for var in mstat.fixed_unused_variables_set(m):
        assert var.name in fixed_unused_var_set

    unfixed_unused_var_set = {
        "fs.heater.control_volume.material_flow_dx[0.0,0.0,Vap,H2]",
        "fs.heater.control_volume.material_flow_dx[0.0,0.0,Vap,H2O]",
        "fs.heater.control_volume.material_flow_dx[0.0,0.0,Vap,Ar]",
        "fs.heater.control_volume.material_flow_dx[0.0,0.0,Vap,N2]",
        "fs.heater.control_volume.enthalpy_flow_dx[0.0,0.0,Vap]",
        "fs.heater.control_volume.pressure_dx[0.0,0.0]",
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

    initializer = Heater1DInitializer(
        solver="ipopt",
        solver_options=optarg
    )
    initializer.initialize(model=m.fs.heater)

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=False)
    assert pyo.value(m.fs.heater.outlet.temperature[0]) == pytest.approx(
        959.55, abs=1e-1
    )


@pytest.mark.integration
def test_units(model_no_dP):
    assert_units_consistent(model_no_dP.fs.heater)


@pytest.fixture
def model_dP():
    m = _create_model(pressure_drop=True)
    return m


@pytest.mark.component
def test_initialization_dP(model_dP):
    m = model_dP

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=True)

    initializer = Heater1DInitializer(
        solver="ipopt",
        solver_options=optarg
    )
    initializer.initialize(model=m.fs.heater)

    assert degrees_of_freedom(m) == 0
    _check_model_statistics(m, deltaP=True)

    assert pyo.value(m.fs.heater.outlet.temperature[0]) == pytest.approx(
        959.55, abs=1e-1
    )
    assert pyo.value(m.fs.heater.outlet.pressure[0]) == pytest.approx(119762.3, abs=1)


@pytest.mark.integration
def test_units_dP(model_dP):
    assert_units_consistent(model_dP.fs.heater)
