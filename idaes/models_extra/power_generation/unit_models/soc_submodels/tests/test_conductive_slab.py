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


# Class for the electrolyte and interconnect
def common_components(nt, nz, nx):
    return {
        pyo.Var: {
            "temperature_z": nz * nt,
            "current_density": nz * nt,
            "length_z": 1,
            "length_y": 1,
            "temperature_deviation_x0": nz * nt,
            "heat_flux_x0": nz * nt,
            "temperature_deviation_x1": nz * nt,
            "heat_flux_x1": nz * nt,
            "temperature_deviation_x": nx * nz * nt,
            "length_x": 1,
            "resistivity_log_preexponential_factor": 1,
            "resistivity_thermal_exponent_dividend": 1,
            "heat_capacity": 1,
            "density": 1,
            "thermal_conductivity": 1,
        },
        pyo.Constraint: {
            "heat_flux_x0_eqn": nz * nt,
            "heat_flux_x1_eqn": nz * nt,
            "energy_balance_solid_eqn": nx * nz * nt,
        },
        pyo.Expression: {
            "temperature_x0": nz * nt,
            "temperature": nx * nz * nt,
            "temperature_x1": nz * nt,
            "dz": nz,
            "dx": nx,
            "node_volume": nx * nz,
            "zface_area": nx,
            "xface_area": nz,
            "current": nz * nt,
            "dTdx": (nx + 1) * nz * nt,
            "dTdz": (nz + 1) * nx * nt,
            "heat_flux_x": (nx + 1) * nz * nt,
            "heat_flux_z": (nz + 1) * nx * nt,
            "resistivity": nx * nz * nt,
            "resistance": nx * nz * nt,
            "voltage_drop": nx * nz * nt,
            "resistance_total": nz * nt,
            "voltage_drop_total": nz * nt,
            "joule_heating": nx * nz * nt,
        },
    }


@pytest.fixture
def modelNoHoldup():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=[0, 1], time_units=pyo.units.s)
    m.fs.slab = soc.SocConductiveSlab(
        has_holdup=False,
        control_volume_zfaces=np.linspace(0, 1, 4).tolist(),
        control_volume_xfaces=np.linspace(0, 1, 5).tolist(),
    )
    slab = m.fs.slab

    m.fs.slab.length_x.fix(1.4e-4)
    m.fs.slab.heat_capacity.fix(470)
    m.fs.slab.density.fix(5160)
    m.fs.slab.thermal_conductivity.fix(2.16)
    m.fs.slab.resistivity_log_preexponential_factor.fix(pyo.log(1.07e-4))
    m.fs.slab.resistivity_thermal_exponent_dividend.fix(7237)

    slab.temperature_deviation_x0.fix(0)
    slab.heat_flux_x0.fix(0)
    slab.length_y.fix(0.08)
    slab.length_z.fix(0.08)
    slab.temperature_z.fix(1000)
    slab.current_density.fix(0)

    return m


@pytest.fixture
def modelHoldupNotDynamic():
    time_set = [0]
    zfaces = np.linspace(0, 1, 8).tolist()
    xfaces = np.linspace(0, 1, 12).tolist()
    m = soc_testing._cell_flowsheet_model(
        dynamic=False,
        time_set=time_set,
        zfaces=zfaces,
    )
    iznodes = m.fs.iznodes
    tset = m.fs.config.time
    m.fs.temperature_deviation_x0 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.K
    )
    m.fs.heat_flux_x0 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.W / pyo.units.m**2
    )
    m.fs.temperature_deviation_x1 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.K
    )
    m.fs.heat_flux_x1 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.W / pyo.units.m**2
    )

    m.fs.slab = soc.SocConductiveSlab(
        has_holdup=True,
        control_volume_zfaces=zfaces,
        control_volume_xfaces=xfaces,
        length_z=m.fs.length_z,
        length_y=m.fs.length_y,
        current_density=m.fs.current_density,
        temperature_z=m.fs.temperature_z,
        temperature_deviation_x0=m.fs.temperature_deviation_x0,
        heat_flux_x0=m.fs.heat_flux_x0,
        temperature_deviation_x1=m.fs.temperature_deviation_x1,
        heat_flux_x1=m.fs.heat_flux_x1,
    )
    m.fs.temperature_deviation_x0.fix(0)
    m.fs.heat_flux_x0.fix(0)

    m.fs.slab.length_x.fix(1.4e-4)
    m.fs.slab.heat_capacity.fix(470)
    m.fs.slab.density.fix(5160)
    m.fs.slab.thermal_conductivity.fix(2.16)
    m.fs.slab.resistivity_log_preexponential_factor.fix(pyo.log(1.07e-4))
    m.fs.slab.resistivity_thermal_exponent_dividend.fix(7237)

    return m


@pytest.mark.build
@pytest.mark.unit
def test_build_model_no_holdup(modelNoHoldup):
    slab = modelNoHoldup.fs.slab
    nx = len(slab.ixnodes)
    nz = len(slab.iznodes)
    nt = len(slab.flowsheet().time)
    comp_dict = common_components(nt, nz, nx)
    comp_dict[pyo.Param] = {}
    comp_dict[pyo.Param]["dcedt_solid"] = nt * nx * nz
    soc_testing._build_test_utility(
        slab,
        comp_dict=comp_dict,
    )

    assert degrees_of_freedom(slab) == 0


@pytest.mark.build
@pytest.mark.unit
def test_build_model_holdup_not_dynamic(modelHoldupNotDynamic):
    slab = modelHoldupNotDynamic.fs.slab
    nx = len(slab.ixnodes)
    nz = len(slab.iznodes)
    nt = len(slab.flowsheet().time)
    comp_dict = common_components(nt, nz, nx)
    comp_dict[pyo.Var]["int_energy_density_solid"] = nx * nz * nt
    comp_dict[pyo.Param] = {}
    comp_dict[pyo.Param]["dcedt_solid"] = nt * nx * nz
    comp_dict[pyo.Constraint]["int_energy_density_solid_eqn"] = nt * nx * nz
    soc_testing._build_test_utility(
        slab,
        comp_dict=comp_dict,
        references=[
            "temperature_z",
            "current_density",
            "length_z",
            "length_y",
            "temperature_deviation_x0",
            "heat_flux_x0",
            "temperature_deviation_x1",
            "heat_flux_x1",
        ],
    )

    assert degrees_of_freedom(slab) == 0


@pytest.mark.component
def test_initialization_exception(modelNoHoldup):
    with pytest.raises(
        NotImplementedError,
        match="Initialization for the conductive_slab unit model is not implemented because there is no obvious "
        "set of boundary conditions to fix during solid_oxide_cell initialization, and it is not meant to be "
        "initialized in isolation.",
    ):
        modelNoHoldup.fs.slab.initialize_build()


# @pytest.mark.component
# def test_units(modelHoldupNotDynamic):
#     # TODO come back to this later
#     return
