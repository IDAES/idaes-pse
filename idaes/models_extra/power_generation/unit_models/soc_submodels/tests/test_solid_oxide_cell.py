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
from idaes.core import FlowsheetBlock, UnitModelBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern
import idaes.models_extra.power_generation.unit_models.soc_submodels as soc
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
import idaes.models_extra.power_generation.unit_models.soc_submodels.testing as soc_testing

solver = pyo.SolverFactory("ipopt")

def build_tester(cell, nt, nz):
    soc_testing._build_test_utility(
        cell,
        comp_dict={
            pyo.Var: {
                "current_density": nz * nt,
                "potential": nt,
                "temperature_z": nz * nt,
                "length_z": 1,
                "length_y": 1,
            },
            pyo.Constraint: {
                "mean_temperature_eqn": nz * nt,
                "potential_eqn": nz * nt,
                "no_qflux_fuel_interconnect_eqn": nz * nt,
                "no_qflux_oxygen_interconnect_eqn": nz * nt,
            },
            pyo.Expression: {
                "eta_contact": nz * nt,
                "eta_ohm": nz * nt,
                "electrical_work": 1,
            },
        },
    )

@pytest.fixture
def model():
    time_set = [0]
    zfaces = np.linspace(0, 1, 4).tolist()
    xfaces_electrode = np.linspace(0, 1, 6).tolist()
    xfaces_electrolyte = np.linspace(0, 1, 8).tolist()

    fuel_comps = ["H2", "H2O", "N2"]
    fuel_tpb_stoich_dict = {"H2": -0.5, "H2O": 0.5, "N2": 0, "Vac": 0.5, "O^2-": -0.5}
    oxygen_comps = ["O2", "N2"]
    oxygen_tpb_stoich_dict = {"O2": -0.25, "N2": 0, "Vac": -0.5, "O^2-": 0.5}
    
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": False,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    m.fs.cell = soc.SolidOxideCell(
        default={
            "has_holdup": True,
            "control_volume_zfaces": zfaces,
            "control_volume_xfaces_fuel_electrode": xfaces_electrode,
            "control_volume_xfaces_oxygen_electrode": xfaces_electrode,
            "control_volume_xfaces_electrolyte": xfaces_electrolyte,
            "fuel_component_list": fuel_comps,
            "fuel_tpb_stoich_dict": fuel_tpb_stoich_dict,
            "oxygen_component_list": oxygen_comps,
            "oxygen_tpb_stoich_dict": oxygen_tpb_stoich_dict,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            "include_contact_resistance": True
        }
    )
    return m

@pytest.fixture
def model_no_contact_resistance():
    time_set = [0]
    zfaces = np.linspace(0, 1, 6).tolist()
    xfaces_electrode = np.linspace(0, 1, 4).tolist()
    xfaces_electrolyte = np.linspace(0, 1, 12).tolist()

    fuel_comps = ["H2", "H2O"]
    fuel_tpb_stoich_dict = {"H2": -0.5, "H2O": 0.5, "Vac": 0.5, "O^2-": -0.5}
    oxygen_comps = ["O2", "N2", "H2O"]
    oxygen_tpb_stoich_dict = {"O2": -0.25, "N2": 0, "H2O":0, "Vac": -0.5, "O^2-": 0.5}
    
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": False,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    m.fs.cell = soc.SolidOxideCell(
        default={
            "has_holdup": False,
            "control_volume_zfaces": zfaces,
            "control_volume_xfaces_fuel_electrode": xfaces_electrode,
            "control_volume_xfaces_oxygen_electrode": xfaces_electrode,
            "control_volume_xfaces_electrolyte": xfaces_electrolyte,
            "fuel_component_list": fuel_comps,
            "fuel_tpb_stoich_dict": fuel_tpb_stoich_dict,
            "oxygen_component_list": oxygen_comps,
            "oxygen_tpb_stoich_dict": oxygen_tpb_stoich_dict,
            "flow_pattern": HeatExchangerFlowPattern.cocurrent,
            "include_contact_resistance": False
        }
    )
    return m

@pytest.mark.build
@pytest.mark.unit
def test_build(model):
    cell = model.fs.cell
    nt = len(model.fs.time)
    nz = len(cell.zfaces) - 1
    build_tester(cell, nt, nz)

    channels = [cell.fuel_chan, cell.oxygen_chan]

    for chan in channels:
        assert cell.temperature_z is chan.temperature_z.referent
        assert cell.length_y is chan.length_y.referent
        assert cell.length_z is chan.length_z.referent

    contact_resistors = [
        cell.contact_interconnect_fuel_flow_mesh,
        cell.contact_interconnect_oxygen_flow_mesh,
        cell.contact_flow_mesh_fuel_electrode,
        cell.contact_flow_mesh_oxygen_electrode,
    ]

    for unit in contact_resistors:
        assert cell.temperature_z is unit.temperature_z.referent
        assert cell.current_density is unit.current_density.referent
        assert cell.length_y is unit.length_y.referent
        assert cell.length_z is unit.length_z.referent

    assert (
        cell.fuel_chan.qflux_x0
        is cell.contact_interconnect_fuel_flow_mesh.qflux_x1.referent
    )
    assert (
        cell.fuel_chan.Dtemp_x0
        is cell.contact_interconnect_fuel_flow_mesh.Dtemp.referent
    )

    assert (
        cell.fuel_chan.qflux_x1
        is cell.contact_flow_mesh_fuel_electrode.qflux_x0.referent
    )
    assert (
        cell.fuel_chan.Dtemp_x1
        is cell.contact_flow_mesh_fuel_electrode.Dtemp.referent
    )

    assert (
        cell.oxygen_chan.qflux_x1
        is cell.contact_interconnect_oxygen_flow_mesh.qflux_x0.referent
    )
    assert (
        cell.oxygen_chan.Dtemp_x1
        is cell.contact_interconnect_oxygen_flow_mesh.Dtemp.referent
    )

    assert (
        cell.oxygen_chan.qflux_x0
        is cell.contact_flow_mesh_oxygen_electrode.qflux_x1.referent
    )
    assert (
        cell.oxygen_chan.Dtemp_x0
        is cell.contact_flow_mesh_oxygen_electrode.Dtemp.referent
    )

    electrodes = [cell.fuel_electrode, cell.oxygen_electrode]

    for trode in electrodes:
        assert cell.temperature_z is trode.temperature_z.referent
        assert cell.current_density is trode.current_density.referent
        assert cell.length_y is trode.length_y.referent
        assert cell.length_z is trode.length_z.referent

    assert cell.fuel_chan.Dtemp_x1 is cell.fuel_electrode.Dtemp_x0.referent
    assert (
        cell.contact_flow_mesh_fuel_electrode.qflux_x1
        is cell.fuel_electrode.qflux_x0.referent
    )
    assert cell.fuel_chan.conc is cell.fuel_electrode.conc_ref.referent
    assert cell.fuel_chan.Dconc_x1 is cell.fuel_electrode.Dconc_x0.referent
    assert cell.fuel_chan.dcdt is cell.fuel_electrode.dconc_refdt.referent
    assert cell.fuel_chan.xflux_x1 is cell.fuel_electrode.xflux_x0.referent

    assert cell.oxygen_chan.Dtemp_x0 is cell.oxygen_electrode.Dtemp_x1.referent
    assert (
        cell.contact_flow_mesh_oxygen_electrode.qflux_x0
        is cell.oxygen_electrode.qflux_x1.referent
    )
    assert cell.oxygen_chan.conc is cell.oxygen_electrode.conc_ref.referent
    assert cell.oxygen_chan.Dconc_x0 is cell.oxygen_electrode.Dconc_x1.referent
    assert cell.oxygen_chan.dcdt is cell.oxygen_electrode.dconc_refdt.referent
    assert cell.oxygen_chan.xflux_x0 is cell.oxygen_electrode.xflux_x1.referent

    tpb_list = [cell.fuel_tpb, cell.oxygen_tpb]

    for tpb in tpb_list:
        assert cell.temperature_z is tpb.temperature_z.referent
        assert cell.current_density is tpb.current_density.referent
        assert cell.length_y is tpb.length_y.referent
        assert cell.length_z is tpb.length_z.referent

    assert cell.fuel_tpb.Dtemp.referent is cell.fuel_electrode.Dtemp_x1
    assert cell.fuel_tpb.qflux_x0.referent is cell.fuel_electrode.qflux_x1
    assert cell.fuel_tpb.conc_ref.referent is cell.fuel_chan.conc
    assert cell.fuel_tpb.Dconc.referent is cell.fuel_electrode.Dconc_x1

    assert cell.oxygen_tpb.Dtemp.referent is cell.oxygen_electrode.Dtemp_x0
    assert cell.oxygen_tpb.qflux_x1.referent is cell.oxygen_electrode.qflux_x0
    assert cell.oxygen_tpb.conc_ref.referent is cell.oxygen_chan.conc
    assert cell.oxygen_tpb.Dconc.referent is cell.oxygen_electrode.Dconc_x0

    assert cell.temperature_z is cell.electrolyte.temperature_z.referent
    assert cell.current_density is cell.electrolyte.current_density.referent
    assert cell.length_y is cell.electrolyte.length_y.referent
    assert cell.length_z is cell.electrolyte.length_z.referent
    assert cell.fuel_electrode.Dtemp_x1 is cell.electrolyte.Dtemp_x0.referent
    assert cell.oxygen_electrode.Dtemp_x0 is cell.electrolyte.Dtemp_x1.referent
    assert cell.fuel_tpb.qflux_x1 is cell.electrolyte.qflux_x0.referent
    assert cell.oxygen_tpb.qflux_x0 is cell.electrolyte.qflux_x1.referent

@pytest.mark.build
@pytest.mark.unit
def test_build_no_contact_resistance(model_no_contact_resistance):
    cell = model_no_contact_resistance.fs.cell
    nt = len(model_no_contact_resistance.fs.time)
    nz = len(cell.zfaces) - 1
    build_tester(cell, nt, nz)

    channels = [cell.fuel_chan, cell.oxygen_chan]

    for chan in channels:
        assert cell.temperature_z is chan.temperature_z.referent
        assert cell.length_y is chan.length_y.referent
        assert cell.length_z is chan.length_z.referent

    electrodes = [cell.fuel_electrode, cell.oxygen_electrode]

    for trode in electrodes:
        assert cell.temperature_z is trode.temperature_z.referent
        assert cell.current_density is trode.current_density.referent
        assert cell.length_y is trode.length_y.referent
        assert cell.length_z is trode.length_z.referent

    assert cell.fuel_chan.Dtemp_x1 is cell.fuel_electrode.Dtemp_x0.referent
    assert (
        cell.fuel_chan.qflux_x1
        is cell.fuel_electrode.qflux_x0.referent
    )
    assert cell.fuel_chan.conc is cell.fuel_electrode.conc_ref.referent
    assert cell.fuel_chan.Dconc_x1 is cell.fuel_electrode.Dconc_x0.referent
    assert cell.fuel_chan.dcdt is cell.fuel_electrode.dconc_refdt.referent
    assert cell.fuel_chan.xflux_x1 is cell.fuel_electrode.xflux_x0.referent

    assert cell.oxygen_chan.Dtemp_x0 is cell.oxygen_electrode.Dtemp_x1.referent
    assert (
        cell.oxygen_chan.qflux_x0
        is cell.oxygen_electrode.qflux_x1.referent
    )
    assert cell.oxygen_chan.conc is cell.oxygen_electrode.conc_ref.referent
    assert cell.oxygen_chan.Dconc_x0 is cell.oxygen_electrode.Dconc_x1.referent
    assert cell.oxygen_chan.dcdt is cell.oxygen_electrode.dconc_refdt.referent
    assert cell.oxygen_chan.xflux_x0 is cell.oxygen_electrode.xflux_x1.referent

    tpb_list = [cell.fuel_tpb, cell.oxygen_tpb]

    for tpb in tpb_list:
        assert cell.temperature_z is tpb.temperature_z.referent
        assert cell.current_density is tpb.current_density.referent
        assert cell.length_y is tpb.length_y.referent
        assert cell.length_z is tpb.length_z.referent

    assert cell.fuel_tpb.Dtemp.referent is cell.fuel_electrode.Dtemp_x1
    assert cell.fuel_tpb.qflux_x0.referent is cell.fuel_electrode.qflux_x1
    assert cell.fuel_tpb.conc_ref.referent is cell.fuel_chan.conc
    assert cell.fuel_tpb.Dconc.referent is cell.fuel_electrode.Dconc_x1

    assert cell.oxygen_tpb.Dtemp.referent is cell.oxygen_electrode.Dtemp_x0
    assert cell.oxygen_tpb.qflux_x1.referent is cell.oxygen_electrode.qflux_x0
    assert cell.oxygen_tpb.conc_ref.referent is cell.oxygen_chan.conc
    assert cell.oxygen_tpb.Dconc.referent is cell.oxygen_electrode.Dconc_x0

    assert cell.temperature_z is cell.electrolyte.temperature_z.referent
    assert cell.current_density is cell.electrolyte.current_density.referent
    assert cell.length_y is cell.electrolyte.length_y.referent
    assert cell.length_z is cell.electrolyte.length_z.referent
    assert cell.fuel_electrode.Dtemp_x1 is cell.electrolyte.Dtemp_x0.referent
    assert cell.oxygen_electrode.Dtemp_x0 is cell.electrolyte.Dtemp_x1.referent
    assert cell.fuel_tpb.qflux_x1 is cell.electrolyte.qflux_x0.referent
    assert cell.oxygen_tpb.qflux_x0 is cell.electrolyte.qflux_x1.referent