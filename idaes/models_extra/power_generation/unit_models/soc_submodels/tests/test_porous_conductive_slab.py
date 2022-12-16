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

__author__ = "Douglas Allan"

import pytest
import numpy as np

import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.models_extra.power_generation.unit_models.soc_submodels as soc
import idaes.models_extra.power_generation.unit_models.soc_submodels.testing as soc_testing


def common_components(nt, nz, nx, ncomp):
    return {
        pyo.Var: {
            "temperature_z": nz * nt,
            "current_density": nz * nt,
            "porosity": 1,
            "tortuosity": 1,
            "length_x": 1,
            "length_z": 1,
            "length_y": 1,
            "conc_mol_comp_ref": nz * nt * ncomp,
            "conc_mol_comp_deviation_x0": nz * nt * ncomp,
            "conc_mol_comp_deviation_x": nz * nx * nt * ncomp,
            "conc_mol_comp_deviation_x1": nz * nt * ncomp,
            "material_flux_x0": nz * nt * ncomp,
            "material_flux_x1": nz * nt * ncomp,
            "temperature_deviation_x0": nz * nt,
            "heat_flux_x0": nz * nt,
            "temperature_deviation_x1": nz * nt,
            "heat_flux_x1": nz * nt,
            "temperature_deviation_x": nx * nz * nt,
            "enth_mol": nx * nz * nt,
            "pressure": nx * nz * nt,
            "mole_frac_comp": nx * nz * nt * ncomp,
            "diff_eff_coeff": nt * nx * nz * ncomp,
            "resistivity_log_preexponential_factor": 1,
            "resistivity_thermal_exponent_dividend": 1,
            "solid_heat_capacity": 1,
            "solid_density": 1,
            "solid_thermal_conductivity": 1,
        },
        pyo.Constraint: {
            "conc_mol_comp_eqn": nt * nz * nx * ncomp,
            "mole_frac_comp_eqn": nt * nz * nx,
            "enth_mol_eqn": nt * nz * nx,
            "material_flux_x0_eqn": nz * nt * ncomp,
            "material_flux_x1_eqn": nz * nt * ncomp,
            "heat_flux_x0_eqn": nz * nt,
            "heat_flux_x1_eqn": nz * nt,
            "diff_eff_coeff_eqn": nt * nx * nz * ncomp,
            "material_balance_eqn": nx * nz * nt * ncomp,
            "energy_balance_solid_eqn": nx * nz * nt,
        },
        pyo.Expression: {
            "temperature_x0": nz * nt,
            "temperature": nx * nz * nt,
            "temperature_x1": nz * nt,
            "conc_mol_comp_x0": nz * nt * ncomp,
            "conc_mol_comp": nz * nx * nt * ncomp,
            "conc_mol_comp_x1": nz * nt * ncomp,
            "diff_eff_coeff_xfaces": nt * (nx + 1) * nz * ncomp,
            "diff_eff_coeff_zfaces": nt * nx * (nz + 1) * ncomp,
            "temperature_xfaces": nt * (nx + 1) * nz,
            "temperature_zfaces": nt * nx * (nz + 1),
            "dconc_mol_compdt": nz * nx * nt * ncomp,
            "vol_mol": nt * nx * nz,
            "dz": nz,
            "dx": nx,
            "node_volume": nx * nz,
            "zface_area": nx,
            "xface_area": nz,
            "dcdx": ncomp * (nx + 1) * nz * nt,
            "dcdz": ncomp * (nz + 1) * nx * nt,
            "material_flux_x": ncomp * (nx + 1) * nz * nt,
            "material_flux_z": ncomp * (nz + 1) * nx * nt,
            "dTdx": (nx + 1) * nz * nt,
            "dTdz": (nz + 1) * nx * nt,
            "heat_flux_x": (nx + 1) * nz * nt,
            "heat_flux_z": (nz + 1) * nx * nt,
            "resistivity": nx * nz * nt,
            "resistance": nx * nz * nt,
            "current": nt * nz,
            "voltage_drop": nx * nz * nt,
            "resistance_total": nz * nt,
            "voltage_drop_total": nz * nt,
            "joule_heating": nx * nz * nt,
        },
    }


def fix_boundary_conditions(electrode):
    electrode.temperature_deviation_x0.fix()
    electrode.conc_mol_comp_ref.fix()
    electrode.conc_mol_comp_deviation_x0.fix()
    electrode.material_flux_x0.fix()
    electrode.heat_flux_x0.fix()


@pytest.fixture
def modelNoHoldup():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0, 1], time_units=pyo.units.s)
    m.fs.oxygen_electrode = soc.PorousConductiveSlab(
        has_holdup=False,
        control_volume_zfaces=np.linspace(0, 1, 4).tolist(),
        control_volume_xfaces=np.linspace(0, 1, 5).tolist(),
        component_list=["O2", "N2"],
    )
    electrode = m.fs.oxygen_electrode
    electrode.length_x.fix(1e-3)
    electrode.porosity.fix(0.48)
    electrode.tortuosity.fix(5.4)
    electrode.solid_heat_capacity.fix(450)
    electrode.solid_density.fix(3210.0)
    electrode.solid_thermal_conductivity.fix(1.86)
    electrode.resistivity_log_preexponential_factor.fix(pyo.log(2.98e-5))
    electrode.resistivity_thermal_exponent_dividend.fix(-1392.0)

    fix_boundary_conditions(electrode)

    electrode.length_y.fix(0.08)
    electrode.length_z.fix(0.08)
    electrode.temperature_z.fix(1000)
    electrode.current_density.fix(0)
    return m


@pytest.fixture
def modelHoldupNotDynamic():
    time_set = [0]
    zfaces = np.linspace(0, 1, 8).tolist()
    xfaces_electrode = np.linspace(0, 1, 12).tolist()
    m = soc_testing._cell_flowsheet_model(
        dynamic=False, time_set=time_set, zfaces=zfaces
    )
    iznodes = m.fs.iznodes
    tset = m.fs.config.time
    comps = m.fs.comps = pyo.Set(initialize=["H2", "H2O", "N2"])
    m.fs.temperature_deviation_x0 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.K
    )
    m.fs.conc_mol_comp_ref = pyo.Var(
        tset, iznodes, comps, initialize=1, units=pyo.units.mol / pyo.units.m**3
    )
    m.fs.conc_mol_comp_deviation_x0 = pyo.Var(
        tset, iznodes, comps, initialize=0, units=pyo.units.mol / pyo.units.m**3
    )
    m.fs.dconc_mol_comp_refdt = pyo.Param(
        tset,
        iznodes,
        comps,
        initialize=0,
        units=pyo.units.mol / (pyo.units.m**3 * pyo.units.s),
    )
    m.fs.material_flux_x0 = pyo.Var(
        tset,
        iznodes,
        comps,
        initialize=0,
        units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
    )
    m.fs.heat_flux_x0 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.W / pyo.units.m**2
    )
    m.fs.fuel_electrode = soc.PorousConductiveSlab(
        has_holdup=True,
        has_gas_holdup=True,
        control_volume_zfaces=zfaces,
        control_volume_xfaces=xfaces_electrode,
        component_list=["H2", "H2O", "N2"],
        length_z=m.fs.length_z,
        length_y=m.fs.length_y,
        conc_mol_comp_ref=m.fs.conc_mol_comp_ref,
        dconc_mol_comp_refdt=m.fs.dconc_mol_comp_refdt,
        conc_mol_comp_deviation_x0=m.fs.conc_mol_comp_deviation_x0,
        material_flux_x0=m.fs.material_flux_x0,
        heat_flux_x0=m.fs.heat_flux_x0,
        temperature_z=m.fs.temperature_z,
        temperature_deviation_x0=m.fs.temperature_deviation_x0,
        current_density=m.fs.current_density,
    )
    electrode = m.fs.fuel_electrode
    electrode.length_x.fix(1e-3)
    electrode.porosity.fix(0.48)
    electrode.tortuosity.fix(5.4)
    electrode.solid_heat_capacity.fix(450)
    electrode.solid_density.fix(3210.0)
    electrode.solid_thermal_conductivity.fix(1.86)
    electrode.resistivity_log_preexponential_factor.fix(pyo.log(2.98e-5))
    electrode.resistivity_thermal_exponent_dividend.fix(-1392.0)

    fix_boundary_conditions(electrode)
    return m


@pytest.mark.build
@pytest.mark.unit
def test_build_modelNoHoldup(modelNoHoldup):
    electrode = modelNoHoldup.fs.oxygen_electrode
    nt = len(electrode.flowsheet().time)
    nz = len(electrode.iznodes)
    nx = len(electrode.xfaces) - 1
    ncomp = len(electrode.component_list)

    comp_dict = common_components(nt, nz, nx, ncomp)
    comp_dict[pyo.Param] = {}
    comp_dict[pyo.Param]["dconc_mol_comp_refdt"] = nt * nz * ncomp
    comp_dict[pyo.Param]["dconc_mol_comp_deviation_xdt"] = nt * nz * nx * ncomp
    comp_dict[pyo.Param]["dcedt"] = nt * nx * nz
    comp_dict[pyo.Param]["dcedt_solid"] = nt * nx * nz

    soc_testing._build_test_utility(
        electrode,
        comp_dict=comp_dict,
    )

    assert degrees_of_freedom(electrode) == 0


@pytest.mark.build
@pytest.mark.unit
def test_build_modelHoldupNotDynamic(modelHoldupNotDynamic):
    electrode = modelHoldupNotDynamic.fs.fuel_electrode
    nt = len(electrode.flowsheet().time)
    nz = len(electrode.iznodes)
    nx = len(electrode.xfaces) - 1
    ncomp = len(electrode.component_list)

    comp_dict = common_components(nt, nz, nx, ncomp)
    comp_dict[pyo.Var]["int_energy_mol"] = nx * nz * nt
    comp_dict[pyo.Var]["int_energy_density"] = nx * nz * nt
    comp_dict[pyo.Var]["int_energy_density_solid"] = nx * nz * nt
    comp_dict[pyo.Param] = {}
    comp_dict[pyo.Param]["dconc_mol_comp_refdt"] = nt * nz * ncomp
    comp_dict[pyo.Param]["dconc_mol_comp_deviation_xdt"] = nt * nz * nx * ncomp
    comp_dict[pyo.Param]["dcedt"] = nt * nx * nz
    comp_dict[pyo.Param]["dcedt_solid"] = nt * nx * nz
    comp_dict[pyo.Constraint]["int_energy_mol_eqn"] = nt * nx * nz
    comp_dict[pyo.Constraint]["int_energy_density_eqn"] = nt * nx * nz
    comp_dict[pyo.Constraint]["int_energy_density_solid_eqn"] = nt * nx * nz

    soc_testing._build_test_utility(
        electrode,
        comp_dict=comp_dict,
        references=[
            "temperature_z",
            "current_density",
            "length_z",
            "length_y",
            "conc_mol_comp_ref",
            "conc_mol_comp_deviation_x0",
            "dconc_mol_comp_refdt",
            "material_flux_x0",
            "temperature_deviation_x0",
            "heat_flux_x0",
        ],
    )

    assert degrees_of_freedom(electrode) == 0
