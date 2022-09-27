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
from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern
import idaes.models_extra.power_generation.unit_models.soc_submodels as soc
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
                "no_heat_flux_fuel_interconnect_eqn": nz * nt,
                "no_heat_flux_oxygen_interconnect_eqn": nz * nt,
            },
            pyo.Expression: {
                "dz": nz,
                "xface_area": nz,
                "cell_area": 1,
                "voltage_drop_contact": nz * nt,
                "voltage_drop_ohmic": nz * nt,
                "electrical_work": 1,
                "average_current_density": 1,
                "voltage_drop_interconnect": nz * nt,
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
    fuel_triple_phase_boundary_stoich_dict = {
        "H2": -0.5,
        "H2O": 0.5,
        "Vac": 0.5,
        "O^2-": -0.5,
        "e^-": 1.0,
    }
    oxygen_comps = ["O2", "N2"]
    oxygen_triple_phase_boundary_stoich_dict = {
        "O2": -0.25,
        "N2": 0,
        "Vac": -0.5,
        "O^2-": 0.5,
        "e^-": -1.0,
    }

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=time_set, time_units=pyo.units.s)
    m.fs.cell = soc.SolidOxideCell(
        has_holdup=True,
        control_volume_zfaces=zfaces,
        control_volume_xfaces_fuel_electrode=xfaces_electrode,
        control_volume_xfaces_oxygen_electrode=xfaces_electrode,
        control_volume_xfaces_electrolyte=xfaces_electrolyte,
        fuel_component_list=fuel_comps,
        fuel_triple_phase_boundary_stoich_dict=fuel_triple_phase_boundary_stoich_dict,
        inert_fuel_species_triple_phase_boundary=["N2"],
        oxygen_component_list=oxygen_comps,
        oxygen_triple_phase_boundary_stoich_dict=oxygen_triple_phase_boundary_stoich_dict,
        inert_oxygen_species_triple_phase_boundary=["N2"],
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        include_contact_resistance=True,
    )
    return m


@pytest.fixture
def model_no_contact_resistance():
    time_set = [0]
    zfaces = np.linspace(0, 1, 6).tolist()
    xfaces_electrode = np.linspace(0, 1, 4).tolist()
    xfaces_electrolyte = np.linspace(0, 1, 12).tolist()

    fuel_comps = ["H2", "H2O"]
    fuel_triple_phase_boundary_stoich_dict = {
        "H2": -0.5,
        "H2O": 0.5,
        "Vac": 0.5,
        "O^2-": -0.5,
        "e^-": 1,
    }
    oxygen_comps = ["O2", "N2", "H2O"]
    oxygen_triple_phase_boundary_stoich_dict = {
        "O2": -0.25,
        "H2O": 0,
        "Vac": -0.5,
        "O^2-": 0.5,
        "e^-": -1,
    }

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=time_set, time_units=pyo.units.s)
    m.fs.cell = soc.SolidOxideCell(
        has_holdup=False,
        control_volume_zfaces=zfaces,
        control_volume_xfaces_fuel_electrode=xfaces_electrode,
        control_volume_xfaces_oxygen_electrode=xfaces_electrode,
        control_volume_xfaces_electrolyte=xfaces_electrolyte,
        fuel_component_list=fuel_comps,
        fuel_triple_phase_boundary_stoich_dict=fuel_triple_phase_boundary_stoich_dict,
        oxygen_component_list=oxygen_comps,
        oxygen_triple_phase_boundary_stoich_dict=oxygen_triple_phase_boundary_stoich_dict,
        inert_oxygen_species_triple_phase_boundary=["N2", "H2O"],
        flow_pattern=HeatExchangerFlowPattern.cocurrent,
        include_contact_resistance=False,
    )
    return m


@pytest.mark.build
@pytest.mark.unit
def test_build(model):
    cell = model.fs.cell
    nt = len(model.fs.time)
    nz = len(cell.zfaces) - 1
    build_tester(cell, nt, nz)

    channels = [cell.fuel_channel, cell.oxygen_channel]

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
        cell.fuel_channel.heat_flux_x0
        is cell.contact_interconnect_fuel_flow_mesh.heat_flux_x1.referent
    )
    assert (
        cell.fuel_channel.temperature_deviation_x0
        is cell.contact_interconnect_fuel_flow_mesh.temperature_deviation_x.referent
    )

    assert (
        cell.fuel_channel.heat_flux_x1
        is cell.contact_flow_mesh_fuel_electrode.heat_flux_x0.referent
    )
    assert (
        cell.fuel_channel.temperature_deviation_x1
        is cell.contact_flow_mesh_fuel_electrode.temperature_deviation_x.referent
    )

    assert (
        cell.oxygen_channel.heat_flux_x1
        is cell.contact_interconnect_oxygen_flow_mesh.heat_flux_x0.referent
    )
    assert (
        cell.oxygen_channel.temperature_deviation_x1
        is cell.contact_interconnect_oxygen_flow_mesh.temperature_deviation_x.referent
    )

    assert (
        cell.oxygen_channel.heat_flux_x0
        is cell.contact_flow_mesh_oxygen_electrode.heat_flux_x1.referent
    )
    assert (
        cell.oxygen_channel.temperature_deviation_x0
        is cell.contact_flow_mesh_oxygen_electrode.temperature_deviation_x.referent
    )

    electrodes = [cell.fuel_electrode, cell.oxygen_electrode]

    for trode in electrodes:
        assert cell.temperature_z is trode.temperature_z.referent
        assert cell.current_density is trode.current_density.referent
        assert cell.length_y is trode.length_y.referent
        assert cell.length_z is trode.length_z.referent

    assert (
        cell.fuel_channel.temperature_deviation_x1
        is cell.fuel_electrode.temperature_deviation_x0.referent
    )
    assert (
        cell.contact_flow_mesh_fuel_electrode.heat_flux_x1
        is cell.fuel_electrode.heat_flux_x0.referent
    )
    assert (
        cell.fuel_channel.conc_mol_comp
        is cell.fuel_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.fuel_channel.conc_mol_comp_deviation_x1
        is cell.fuel_electrode.conc_mol_comp_deviation_x0.referent
    )
    assert (
        cell.fuel_channel.dconc_mol_compdt
        is cell.fuel_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.fuel_channel.material_flux_x1
        is cell.fuel_electrode.material_flux_x0.referent
    )

    assert (
        cell.oxygen_channel.temperature_deviation_x0
        is cell.oxygen_electrode.temperature_deviation_x1.referent
    )
    assert (
        cell.contact_flow_mesh_oxygen_electrode.heat_flux_x0
        is cell.oxygen_electrode.heat_flux_x1.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp
        is cell.oxygen_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp_deviation_x0
        is cell.oxygen_electrode.conc_mol_comp_deviation_x1.referent
    )
    assert (
        cell.oxygen_channel.dconc_mol_compdt
        is cell.oxygen_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.oxygen_channel.material_flux_x0
        is cell.oxygen_electrode.material_flux_x1.referent
    )

    tpb_list = [cell.fuel_triple_phase_boundary, cell.oxygen_triple_phase_boundary]

    for tpb in tpb_list:
        assert cell.temperature_z is tpb.temperature_z.referent
        assert cell.current_density is tpb.current_density.referent
        assert cell.length_y is tpb.length_y.referent
        assert cell.length_z is tpb.length_z.referent

    assert (
        cell.fuel_triple_phase_boundary.temperature_deviation_x.referent
        is cell.fuel_electrode.temperature_deviation_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x0.referent
        is cell.fuel_electrode.heat_flux_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.fuel_channel.conc_mol_comp
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.fuel_electrode.conc_mol_comp_deviation_x1
    )

    assert (
        cell.oxygen_triple_phase_boundary.temperature_deviation_x.referent
        is cell.oxygen_electrode.temperature_deviation_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x1.referent
        is cell.oxygen_electrode.heat_flux_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.oxygen_channel.conc_mol_comp
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.oxygen_electrode.conc_mol_comp_deviation_x0
    )

    assert cell.temperature_z is cell.electrolyte.temperature_z.referent
    assert cell.current_density is cell.electrolyte.current_density.referent
    assert cell.length_y is cell.electrolyte.length_y.referent
    assert cell.length_z is cell.electrolyte.length_z.referent
    assert (
        cell.fuel_electrode.temperature_deviation_x1
        is cell.electrolyte.temperature_deviation_x0.referent
    )
    assert (
        cell.oxygen_electrode.temperature_deviation_x0
        is cell.electrolyte.temperature_deviation_x1.referent
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x1
        is cell.electrolyte.heat_flux_x0.referent
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x0
        is cell.electrolyte.heat_flux_x1.referent
    )


@pytest.mark.build
@pytest.mark.unit
def test_build_no_contact_resistance(model_no_contact_resistance):
    cell = model_no_contact_resistance.fs.cell
    nt = len(model_no_contact_resistance.fs.time)
    nz = len(cell.zfaces) - 1
    build_tester(cell, nt, nz)

    channels = [cell.fuel_channel, cell.oxygen_channel]

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

    assert (
        cell.fuel_channel.temperature_deviation_x1
        is cell.fuel_electrode.temperature_deviation_x0.referent
    )
    assert cell.fuel_channel.heat_flux_x1 is cell.fuel_electrode.heat_flux_x0.referent
    assert (
        cell.fuel_channel.conc_mol_comp
        is cell.fuel_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.fuel_channel.conc_mol_comp_deviation_x1
        is cell.fuel_electrode.conc_mol_comp_deviation_x0.referent
    )
    assert (
        cell.fuel_channel.dconc_mol_compdt
        is cell.fuel_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.fuel_channel.material_flux_x1
        is cell.fuel_electrode.material_flux_x0.referent
    )

    assert (
        cell.oxygen_channel.temperature_deviation_x0
        is cell.oxygen_electrode.temperature_deviation_x1.referent
    )
    assert (
        cell.oxygen_channel.heat_flux_x0 is cell.oxygen_electrode.heat_flux_x1.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp
        is cell.oxygen_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp_deviation_x0
        is cell.oxygen_electrode.conc_mol_comp_deviation_x1.referent
    )
    assert (
        cell.oxygen_channel.dconc_mol_compdt
        is cell.oxygen_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.oxygen_channel.material_flux_x0
        is cell.oxygen_electrode.material_flux_x1.referent
    )

    tpb_list = [cell.fuel_triple_phase_boundary, cell.oxygen_triple_phase_boundary]

    for tpb in tpb_list:
        assert cell.temperature_z is tpb.temperature_z.referent
        assert cell.current_density is tpb.current_density.referent
        assert cell.length_y is tpb.length_y.referent
        assert cell.length_z is tpb.length_z.referent

    assert (
        cell.fuel_triple_phase_boundary.temperature_deviation_x.referent
        is cell.fuel_electrode.temperature_deviation_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x0.referent
        is cell.fuel_electrode.heat_flux_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.fuel_channel.conc_mol_comp
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.fuel_electrode.conc_mol_comp_deviation_x1
    )

    assert (
        cell.oxygen_triple_phase_boundary.temperature_deviation_x.referent
        is cell.oxygen_electrode.temperature_deviation_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x1.referent
        is cell.oxygen_electrode.heat_flux_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.oxygen_channel.conc_mol_comp
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.oxygen_electrode.conc_mol_comp_deviation_x0
    )

    assert cell.temperature_z is cell.electrolyte.temperature_z.referent
    assert cell.current_density is cell.electrolyte.current_density.referent
    assert cell.length_y is cell.electrolyte.length_y.referent
    assert cell.length_z is cell.electrolyte.length_z.referent
    assert (
        cell.fuel_electrode.temperature_deviation_x1
        is cell.electrolyte.temperature_deviation_x0.referent
    )
    assert (
        cell.oxygen_electrode.temperature_deviation_x0
        is cell.electrolyte.temperature_deviation_x1.referent
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x1
        is cell.electrolyte.heat_flux_x0.referent
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x0
        is cell.electrolyte.heat_flux_x1.referent
    )


def build_tester_interconnect(cell, nt, nz):
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
            },
            pyo.Expression: {
                "cell_area": 1,
                "dz": nz,
                "xface_area": nz,
                "voltage_drop_contact": nz * nt,
                "voltage_drop_ohmic": nz * nt,
                "electrical_work": 1,
                "average_current_density": 1,
                "voltage_drop_interconnect": nz * nt,
            },
        },
    )


@pytest.fixture
def model_contact_resistance_and_interconnect():
    time_set = [0]
    zfaces = np.linspace(0, 1, 4).tolist()
    xfaces_electrode = np.linspace(0, 1, 6).tolist()
    xfaces_electrolyte = np.linspace(0, 1, 8).tolist()

    fuel_comps = ["H2", "H2O", "N2"]
    fuel_triple_phase_boundary_stoich_dict = {
        "H2": -0.5,
        "H2O": 0.5,
        "Vac": 0.5,
        "O^2-": -0.5,
        "e^-": 1.0,
    }
    oxygen_comps = ["O2", "N2"]
    oxygen_triple_phase_boundary_stoich_dict = {
        "O2": -0.25,
        "N2": 0,
        "Vac": -0.5,
        "O^2-": 0.5,
        "e^-": -1.0,
    }

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=False,
        time_set=time_set,
        time_units=pyo.units.s,
    )
    m.fs.cell = soc.SolidOxideCell(
        has_holdup=True,
        control_volume_zfaces=zfaces,
        control_volume_xfaces_fuel_electrode=xfaces_electrode,
        control_volume_xfaces_oxygen_electrode=xfaces_electrode,
        control_volume_xfaces_electrolyte=xfaces_electrolyte,
        fuel_component_list=fuel_comps,
        fuel_triple_phase_boundary_stoich_dict=fuel_triple_phase_boundary_stoich_dict,
        inert_fuel_species_triple_phase_boundary=["N2"],
        oxygen_component_list=oxygen_comps,
        oxygen_triple_phase_boundary_stoich_dict=oxygen_triple_phase_boundary_stoich_dict,
        inert_oxygen_species_triple_phase_boundary=["N2"],
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        include_contact_resistance=True,
        flux_through_interconnect=True,
        control_volume_xfaces_interconnect=[0, 0.5, 1],
    )
    return m


@pytest.mark.build
@pytest.mark.unit
def test_build_contact_resistance_and_interconnect(
    model_contact_resistance_and_interconnect,
):
    cell = model_contact_resistance_and_interconnect.fs.cell
    nt = len(model_contact_resistance_and_interconnect.fs.time)
    nz = len(cell.zfaces) - 1
    build_tester_interconnect(cell, nt, nz)

    channels = [cell.fuel_channel, cell.oxygen_channel]

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
        cell.fuel_channel.heat_flux_x0
        is cell.contact_interconnect_fuel_flow_mesh.heat_flux_x1.referent
    )
    assert (
        cell.fuel_channel.temperature_deviation_x0
        is cell.contact_interconnect_fuel_flow_mesh.temperature_deviation_x.referent
    )

    assert (
        cell.fuel_channel.heat_flux_x1
        is cell.contact_flow_mesh_fuel_electrode.heat_flux_x0.referent
    )
    assert (
        cell.fuel_channel.temperature_deviation_x1
        is cell.contact_flow_mesh_fuel_electrode.temperature_deviation_x.referent
    )

    assert (
        cell.oxygen_channel.heat_flux_x1
        is cell.contact_interconnect_oxygen_flow_mesh.heat_flux_x0.referent
    )
    assert (
        cell.oxygen_channel.temperature_deviation_x1
        is cell.contact_interconnect_oxygen_flow_mesh.temperature_deviation_x.referent
    )

    assert (
        cell.oxygen_channel.heat_flux_x0
        is cell.contact_flow_mesh_oxygen_electrode.heat_flux_x1.referent
    )
    assert (
        cell.oxygen_channel.temperature_deviation_x0
        is cell.contact_flow_mesh_oxygen_electrode.temperature_deviation_x.referent
    )

    electrodes = [cell.fuel_electrode, cell.oxygen_electrode]

    for trode in electrodes:
        assert cell.temperature_z is trode.temperature_z.referent
        assert cell.current_density is trode.current_density.referent
        assert cell.length_y is trode.length_y.referent
        assert cell.length_z is trode.length_z.referent

    assert (
        cell.fuel_channel.temperature_deviation_x1
        is cell.fuel_electrode.temperature_deviation_x0.referent
    )
    assert (
        cell.contact_flow_mesh_fuel_electrode.heat_flux_x1
        is cell.fuel_electrode.heat_flux_x0.referent
    )
    assert (
        cell.fuel_channel.conc_mol_comp
        is cell.fuel_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.fuel_channel.conc_mol_comp_deviation_x1
        is cell.fuel_electrode.conc_mol_comp_deviation_x0.referent
    )
    assert (
        cell.fuel_channel.dconc_mol_compdt
        is cell.fuel_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.fuel_channel.material_flux_x1
        is cell.fuel_electrode.material_flux_x0.referent
    )

    assert (
        cell.oxygen_channel.temperature_deviation_x0
        is cell.oxygen_electrode.temperature_deviation_x1.referent
    )
    assert (
        cell.contact_flow_mesh_oxygen_electrode.heat_flux_x0
        is cell.oxygen_electrode.heat_flux_x1.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp
        is cell.oxygen_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp_deviation_x0
        is cell.oxygen_electrode.conc_mol_comp_deviation_x1.referent
    )
    assert (
        cell.oxygen_channel.dconc_mol_compdt
        is cell.oxygen_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.oxygen_channel.material_flux_x0
        is cell.oxygen_electrode.material_flux_x1.referent
    )

    tpb_list = [cell.fuel_triple_phase_boundary, cell.oxygen_triple_phase_boundary]

    for tpb in tpb_list:
        assert cell.temperature_z is tpb.temperature_z.referent
        assert cell.current_density is tpb.current_density.referent
        assert cell.length_y is tpb.length_y.referent
        assert cell.length_z is tpb.length_z.referent

    assert (
        cell.fuel_triple_phase_boundary.temperature_deviation_x.referent
        is cell.fuel_electrode.temperature_deviation_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x0.referent
        is cell.fuel_electrode.heat_flux_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.fuel_channel.conc_mol_comp
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.fuel_electrode.conc_mol_comp_deviation_x1
    )

    assert (
        cell.oxygen_triple_phase_boundary.temperature_deviation_x.referent
        is cell.oxygen_electrode.temperature_deviation_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x1.referent
        is cell.oxygen_electrode.heat_flux_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.oxygen_channel.conc_mol_comp
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.oxygen_electrode.conc_mol_comp_deviation_x0
    )

    assert cell.temperature_z is cell.electrolyte.temperature_z.referent
    assert cell.current_density is cell.electrolyte.current_density.referent
    assert cell.length_y is cell.electrolyte.length_y.referent
    assert cell.length_z is cell.electrolyte.length_z.referent
    assert (
        cell.fuel_electrode.temperature_deviation_x1
        is cell.electrolyte.temperature_deviation_x0.referent
    )
    assert (
        cell.oxygen_electrode.temperature_deviation_x0
        is cell.electrolyte.temperature_deviation_x1.referent
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x1
        is cell.electrolyte.heat_flux_x0.referent
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x0
        is cell.electrolyte.heat_flux_x1.referent
    )

    assert cell.temperature_z is cell.interconnect.temperature_z.referent
    assert cell.current_density is cell.interconnect.current_density.referent
    assert cell.length_y is cell.interconnect.length_y.referent
    assert cell.length_z is cell.interconnect.length_z.referent

    assert (
        cell.oxygen_channel.temperature_deviation_x1
        is cell.interconnect.temperature_deviation_x0.referent
    )
    assert (
        cell.fuel_channel.temperature_deviation_x0
        is cell.interconnect.temperature_deviation_x1.referent
    )
    assert (
        cell.contact_interconnect_oxygen_flow_mesh.heat_flux_x1
        is cell.interconnect.heat_flux_x0.referent
    )
    assert (
        cell.contact_interconnect_fuel_flow_mesh.heat_flux_x0
        is cell.interconnect.heat_flux_x1.referent
    )


@pytest.fixture
def model_no_contact_resistance_but_interconnect():
    time_set = [0]
    zfaces = np.linspace(0, 1, 6).tolist()
    xfaces_electrode = np.linspace(0, 1, 4).tolist()
    xfaces_electrolyte = np.linspace(0, 1, 12).tolist()
    xfaces_interconnect = np.linspace(0, 1, 16).tolist()

    fuel_comps = ["H2", "H2O"]
    fuel_triple_phase_boundary_stoich_dict = {
        "H2": -0.5,
        "H2O": 0.5,
        "Vac": 0.5,
        "O^2-": -0.5,
        "e^-": 1,
    }
    oxygen_comps = ["O2", "N2", "H2O"]
    oxygen_triple_phase_boundary_stoich_dict = {
        "O2": -0.25,
        "H2O": 0,
        "Vac": -0.5,
        "O^2-": 0.5,
        "e^-": -1,
    }

    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=False,
        time_set=time_set,
        time_units=pyo.units.s,
    )
    m.fs.cell = soc.SolidOxideCell(
        has_holdup=False,
        control_volume_zfaces=zfaces,
        control_volume_xfaces_fuel_electrode=xfaces_electrode,
        control_volume_xfaces_oxygen_electrode=xfaces_electrode,
        control_volume_xfaces_electrolyte=xfaces_electrolyte,
        fuel_component_list=fuel_comps,
        fuel_triple_phase_boundary_stoich_dict=fuel_triple_phase_boundary_stoich_dict,
        # inert_fuel_species_triple_phase_boundary=[], Test default
        oxygen_component_list=oxygen_comps,
        oxygen_triple_phase_boundary_stoich_dict=oxygen_triple_phase_boundary_stoich_dict,
        inert_oxygen_species_triple_phase_boundary=["N2", "H2O"],
        flow_pattern=HeatExchangerFlowPattern.cocurrent,
        flux_through_interconnect=True,
        control_volume_xfaces_interconnect=xfaces_interconnect,
    )
    return m


@pytest.mark.build
@pytest.mark.unit
def test_build_no_contact_resistance_but_interconnect(
    model_no_contact_resistance_but_interconnect,
):
    cell = model_no_contact_resistance_but_interconnect.fs.cell
    nt = len(model_no_contact_resistance_but_interconnect.fs.time)
    nz = len(cell.zfaces) - 1
    build_tester_interconnect(cell, nt, nz)

    channels = [cell.fuel_channel, cell.oxygen_channel]

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

    assert (
        cell.fuel_channel.temperature_deviation_x1
        is cell.fuel_electrode.temperature_deviation_x0.referent
    )
    assert cell.fuel_channel.heat_flux_x1 is cell.fuel_electrode.heat_flux_x0.referent
    assert (
        cell.fuel_channel.conc_mol_comp
        is cell.fuel_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.fuel_channel.conc_mol_comp_deviation_x1
        is cell.fuel_electrode.conc_mol_comp_deviation_x0.referent
    )
    assert (
        cell.fuel_channel.dconc_mol_compdt
        is cell.fuel_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.fuel_channel.material_flux_x1
        is cell.fuel_electrode.material_flux_x0.referent
    )

    assert (
        cell.oxygen_channel.temperature_deviation_x0
        is cell.oxygen_electrode.temperature_deviation_x1.referent
    )
    assert (
        cell.oxygen_channel.heat_flux_x0 is cell.oxygen_electrode.heat_flux_x1.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp
        is cell.oxygen_electrode.conc_mol_comp_ref.referent
    )
    assert (
        cell.oxygen_channel.conc_mol_comp_deviation_x0
        is cell.oxygen_electrode.conc_mol_comp_deviation_x1.referent
    )
    assert (
        cell.oxygen_channel.dconc_mol_compdt
        is cell.oxygen_electrode.dconc_mol_comp_refdt.referent
    )
    assert (
        cell.oxygen_channel.material_flux_x0
        is cell.oxygen_electrode.material_flux_x1.referent
    )

    tpb_list = [cell.fuel_triple_phase_boundary, cell.oxygen_triple_phase_boundary]

    for tpb in tpb_list:
        assert cell.temperature_z is tpb.temperature_z.referent
        assert cell.current_density is tpb.current_density.referent
        assert cell.length_y is tpb.length_y.referent
        assert cell.length_z is tpb.length_z.referent

    assert (
        cell.fuel_triple_phase_boundary.temperature_deviation_x.referent
        is cell.fuel_electrode.temperature_deviation_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x0.referent
        is cell.fuel_electrode.heat_flux_x1
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.fuel_channel.conc_mol_comp
    )
    assert (
        cell.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.fuel_electrode.conc_mol_comp_deviation_x1
    )

    assert (
        cell.oxygen_triple_phase_boundary.temperature_deviation_x.referent
        is cell.oxygen_electrode.temperature_deviation_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x1.referent
        is cell.oxygen_electrode.heat_flux_x0
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_ref.referent
        is cell.oxygen_channel.conc_mol_comp
    )
    assert (
        cell.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.referent
        is cell.oxygen_electrode.conc_mol_comp_deviation_x0
    )

    assert cell.temperature_z is cell.electrolyte.temperature_z.referent
    assert cell.current_density is cell.electrolyte.current_density.referent
    assert cell.length_y is cell.electrolyte.length_y.referent
    assert cell.length_z is cell.electrolyte.length_z.referent
    assert (
        cell.fuel_electrode.temperature_deviation_x1
        is cell.electrolyte.temperature_deviation_x0.referent
    )
    assert (
        cell.oxygen_electrode.temperature_deviation_x0
        is cell.electrolyte.temperature_deviation_x1.referent
    )
    assert (
        cell.fuel_triple_phase_boundary.heat_flux_x1
        is cell.electrolyte.heat_flux_x0.referent
    )
    assert (
        cell.oxygen_triple_phase_boundary.heat_flux_x0
        is cell.electrolyte.heat_flux_x1.referent
    )

    assert (
        cell.oxygen_channel.temperature_deviation_x1
        is cell.interconnect.temperature_deviation_x0.referent
    )
    assert (
        cell.fuel_channel.temperature_deviation_x0
        is cell.interconnect.temperature_deviation_x1.referent
    )
    assert cell.oxygen_channel.heat_flux_x1 is cell.interconnect.heat_flux_x0.referent
    assert cell.fuel_channel.heat_flux_x0 is cell.interconnect.heat_flux_x1.referent
