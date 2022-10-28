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
"""
Test SolidOxideCell model using parameters from the literature. Compare to both
literature data and cached data from the SolidOxideCell model (to detect model
changes having a significant effect on the model outputs)

This test is based on two papers. The first is Herring et al., in which results
of a series of water electrolysis experiments on a ten cell stack is reported.
Current-voltage curves for the stack are reported at a variety of flow rates,
stream compositions, and temperatures. Some of these results are then
replicated in the CFD software package FLUENT.

The second paper is Kazempoor and Braun, a modeling paper in which some of
these experimental results are replicated, and parameter values for this
replication are reported. These parameter values are used to validate the
SolidOxideCell model.

Several factors limit the ability of the SolidOxideCell model to fully
replicate the Herring et al. data:
    1) The experimental stack operates using cross flow, whereas the
    SolidOxideCell model is limited to counter and co current operation.
    2) The SolidOxideCell model is adiabatic, whereas the stack in Herring et
    al. was operated in a furnace, and received a substantial amount of heat
    3) The open cell potential in Herring et al. was somewhat below the Nernst
    potential, whereas the SolidOxideCell model does not attempt to model
    current leaks
    
However, Kazempoor and Braun use a cocurrent and (probably?) adiabatic model to
partially replicate the Herring et al. experimental data. As a result, the
SolidOxideCell can replicate their replication dataset more closely than the
actual experimental data. Because that model has a slightly different structure
than the SolidOxideCell and the open circuit potential it reports agrees with the
experimental data instead of the Nernst potential calculated both by the thermo
model used by the SolidOxideCell, which agrees with the theoretical open cell
potential reported by Herring et al., there is not a perfect match between
the Kazempoor dataset and that generated by the SolidOxideCell. 

J. Stephen Herring, James E. O’Brien, Carl M. Stoots, G.L. Hawkes, Joseph J.
Hartvigsen, Mehrdad Shahnam, Progress in high-temperature electrolysis for
hydrogen production using planar SOFC technology, International Journal of
Hydrogen Energy, Volume 32, Issue 4, 2007, Pages 440-450, ISSN 0360-3199,
https://doi.org/10.1016/j.ijhydene.2006.06.061.
(https://www.sciencedirect.com/science/article/pii/S0360319906002709)

P. Kazempoor, R.J. Braun, Model validation and performance analysis of
regenerative solid oxide cells: Electrolytic operation, International Journal
of Hydrogen Energy, Volume 39, Issue 6, 2014, Pages 2669-2684, ISSN 0360-3199,
https://doi.org/10.1016/j.ijhydene.2013.12.010.
(https://www.sciencedirect.com/science/article/pii/S036031991302908X)
"""

__author__ = "Douglas Allan"

import os
import numpy as np
from scipy.interpolate import interp1d
import pytest
import pandas as pd
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir

from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
import idaes.core.util.scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import (
    get_prop,
    EosType,
)
from idaes.models_extra.power_generation.unit_models.soc_submodels import (
    SolidOxideCell,
    SolidOxideModuleSimple,
)
from idaes.models_extra.power_generation.unit_models.soc_submodels.common import (
    _comp_enthalpy_expr,
    _comp_entropy_expr,
)

from idaes.models.unit_models.heat_exchanger import HeatExchangerFlowPattern

data_cache = os.sep.join([this_file_dir(), "data_cache"])


def cccm_to_mps(V):
    # Takes cccm, returns moles per second
    V = V * 1e-6 / 60
    P = 101300
    R = 8.31
    T = 273.15
    return V * P / (R * T)


def mps_to_cccm(n):
    P = 101300
    R = 8.31
    T = 273.15
    V = n * R * T / P
    return V * 1e6 * 60


N_cell = 10
V_min = 8
P = 85910

zfaces = np.linspace(0, 1, 11).tolist()

xfaces_electrode = [0.0, 1.0]
xfaces_electrolyte = [0.0, 1.0]
xfaces_interconnect = [0.0, 1.0]

fuel_comps = ["H2", "H2O", "N2"]
fuel_stoich_dict = {
    "H2": -0.5,
    "H2O": 0.5,
    "N2": 0,
    "Vac": 0.5,
    "O^2-": -0.5,
    "e^-": 1.0,
}
oxygen_comps = ["O2", "N2"]
oxygen_stoich_dict = {"O2": -0.25, "Vac": -0.5, "O^2-": 0.5, "e^-": -1.0}

cell_config = {
    "has_holdup": False,
    "control_volume_zfaces": zfaces,
    "control_volume_xfaces_fuel_electrode": xfaces_electrode,
    "control_volume_xfaces_oxygen_electrode": xfaces_electrode,
    "control_volume_xfaces_electrolyte": xfaces_electrolyte,
    "control_volume_xfaces_interconnect": xfaces_interconnect,
    "fuel_component_list": fuel_comps,
    "fuel_triple_phase_boundary_stoich_dict": fuel_stoich_dict,
    "inert_fuel_species_triple_phase_boundary": ["N2"],
    "oxygen_component_list": oxygen_comps,
    "oxygen_triple_phase_boundary_stoich_dict": oxygen_stoich_dict,
    "inert_oxygen_species_triple_phase_boundary": ["N2"],
    "flow_pattern": HeatExchangerFlowPattern.cocurrent,
    "include_temperature_x_thermo": True,
    "include_contact_resistance": True,
    "flux_through_interconnect": True,
}


def fix_cell_parameters(cell):
    cell.fuel_channel.length_x.fix(1.09e-3)
    cell.oxygen_channel.length_x.fix(1.09e-3)

    cell.length_y.fix(0.08)

    cell.length_z.fix(0.08)
    cell.fuel_channel.heat_transfer_coefficient.fix(100)

    cell.contact_flow_mesh_fuel_electrode.log_preexponential_factor.fix(
        pyo.log(0.55e-4)
    )
    cell.contact_flow_mesh_oxygen_electrode.log_preexponential_factor.fix(
        pyo.log(0.55e-4)
    )
    cell.contact_interconnect_fuel_flow_mesh.log_preexponential_factor.fix(
        pyo.log(0.15e-4)
    )
    cell.contact_interconnect_oxygen_flow_mesh.log_preexponential_factor.fix(
        pyo.log(0.15e-4)
    )

    cell.contact_flow_mesh_fuel_electrode.thermal_exponent_dividend.fix(0)
    cell.contact_flow_mesh_fuel_electrode.contact_fraction.fix(1)
    cell.contact_flow_mesh_oxygen_electrode.thermal_exponent_dividend.fix(0)
    cell.contact_flow_mesh_oxygen_electrode.contact_fraction.fix(1)

    cell.contact_interconnect_fuel_flow_mesh.thermal_exponent_dividend.fix(0)
    cell.contact_interconnect_fuel_flow_mesh.contact_fraction.fix(1)
    cell.contact_interconnect_oxygen_flow_mesh.thermal_exponent_dividend.fix(0)
    cell.contact_interconnect_oxygen_flow_mesh.contact_fraction.fix(1)

    cell.oxygen_channel.heat_transfer_coefficient.fix(100)

    cell.fuel_electrode.length_x.fix(23e-6)
    cell.fuel_electrode.porosity.fix(0.37)
    cell.fuel_electrode.tortuosity.fix(3)
    cell.fuel_electrode.solid_heat_capacity.fix(450)
    cell.fuel_electrode.solid_density.fix(3210.0)
    cell.fuel_electrode.solid_thermal_conductivity.fix(2.16)
    cell.fuel_electrode.resistivity_log_preexponential_factor.fix(pyo.log(8.856e-6))
    cell.fuel_electrode.resistivity_thermal_exponent_dividend.fix(0)

    cell.oxygen_electrode.length_x.fix(21e-6)
    cell.oxygen_electrode.porosity.fix(0.37)
    cell.oxygen_electrode.tortuosity.fix(3.0)
    cell.oxygen_electrode.solid_heat_capacity.fix(430)
    cell.oxygen_electrode.solid_density.fix(3030)
    cell.oxygen_electrode.solid_thermal_conductivity.fix(2.16)
    cell.oxygen_electrode.resistivity_log_preexponential_factor.fix(pyo.log(1.425e-4))
    cell.oxygen_electrode.resistivity_thermal_exponent_dividend.fix(0)

    cell.electrolyte.length_x.fix(1.4e-4)
    cell.electrolyte.heat_capacity.fix(470)
    cell.electrolyte.density.fix(5160)
    cell.electrolyte.thermal_conductivity.fix(2.16)
    cell.electrolyte.resistivity_log_preexponential_factor.fix(pyo.log(1.07e-4))
    cell.electrolyte.resistivity_thermal_exponent_dividend.fix(7237)

    cell.interconnect.length_x.fix(4.6e-4)
    cell.interconnect.heat_capacity.fix(630)
    cell.interconnect.density.fix(7700)
    cell.interconnect.thermal_conductivity.fix(25)
    cell.interconnect.resistivity_log_preexponential_factor.fix(pyo.log(1.2e-6))
    cell.interconnect.resistivity_thermal_exponent_dividend.fix(0)

    # Values from Kazempoor and Braun
    cell.fuel_triple_phase_boundary.exchange_current_log_preexponential_factor.fix(
        pyo.log(5.5e8)
    )
    cell.fuel_triple_phase_boundary.exchange_current_activation_energy.fix(100e3)
    cell.fuel_triple_phase_boundary.activation_potential_alpha1.fix(1)
    cell.fuel_triple_phase_boundary.activation_potential_alpha2.fix(1)
    cell.fuel_triple_phase_boundary.exchange_current_exponent_comp["H2"].fix(0.5)
    cell.fuel_triple_phase_boundary.exchange_current_exponent_comp["H2O"].fix(0.1)

    cell.oxygen_triple_phase_boundary.exchange_current_log_preexponential_factor.fix(
        pyo.log(7e8)
    )
    cell.oxygen_triple_phase_boundary.exchange_current_activation_energy.fix(110e3)
    cell.oxygen_triple_phase_boundary.activation_potential_alpha1.fix(1)
    cell.oxygen_triple_phase_boundary.activation_potential_alpha2.fix(1)
    cell.oxygen_triple_phase_boundary.exchange_current_exponent_comp["O2"].fix(0.5)


def model_func():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=False,
        time_set=[0],
        time_units=pyo.units.s,
    )

    m.fs.propertiesIapws95 = iapws95.Iapws95ParameterBlock()
    m.fs.prop_Iapws95 = iapws95.Iapws95StateBlock(parameters=m.fs.propertiesIapws95)

    m.fs.cell = SolidOxideCell(**cell_config)

    fix_cell_parameters(m.fs.cell)

    m.fs.cell.oxygen_inlet.flow_mol[0].fix(cccm_to_mps(3500) / N_cell)
    m.fs.cell.oxygen_inlet.mole_frac_comp[0, "O2"].fix(0.21)
    m.fs.cell.oxygen_inlet.mole_frac_comp[0, "N2"].fix(0.79)

    m.fs.cell.recursive_scaling()

    return m


@pytest.fixture
def model():
    # Want to be able to call model_func when calling this test directly
    # so that the reference dataset can be recreated if necessary
    return model_func()


@pytest.fixture
def model_stack():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        dynamic=False,
        time_set=[0],
        time_units=pyo.units.s,
    )

    m.fs.propertiesIapws95 = iapws95.Iapws95ParameterBlock()
    m.fs.prop_Iapws95 = iapws95.Iapws95StateBlock(parameters=m.fs.propertiesIapws95)

    m.fs.oxygen_params = GenericParameterBlock(
        **get_prop(oxygen_comps, {"Vap"}, eos=EosType.IDEAL),
    )
    m.fs.fuel_params = GenericParameterBlock(
        **get_prop(fuel_comps, {"Vap"}, eos=EosType.IDEAL),
    )

    m.fs.stack = SolidOxideModuleSimple(
        solid_oxide_cell_config=cell_config,
        fuel_property_package=m.fs.fuel_params,
        oxygen_property_package=m.fs.oxygen_params,
    )

    fix_cell_parameters(m.fs.stack.solid_oxide_cell)
    m.fs.stack.number_cells.fix(10)
    m.fs.stack.oxygen_inlet.flow_mol[0].fix(cccm_to_mps(3500))
    m.fs.stack.oxygen_inlet.mole_frac_comp[0, "O2"].fix(0.21)
    m.fs.stack.oxygen_inlet.mole_frac_comp[0, "N2"].fix(0.79)

    for pp in [m.fs.fuel_params, m.fs.oxygen_params]:
        pp.set_default_scaling("flow_mol", 1e3)
        pp.set_default_scaling("flow_mol_phase", 1e3)
        pp.set_default_scaling("mole_frac_comp", 10)
        pp.set_default_scaling("mole_frac_phase_comp", 10)
        pp.set_default_scaling("enth_mol_phase", 1e-3)
        pp.set_default_scaling("pressure", 1e-5)
        pp.set_default_scaling("temperature", 1e-2)

    return m


@pytest.mark.component
def test_initialization_cell(model):
    m = model
    cell = m.fs.cell

    cell.potential.fix(1.288)
    cell.fuel_inlet.temperature[0].fix(1103.15)
    cell.fuel_inlet.pressure.fix(85910)
    cell.fuel_inlet.flow_mol[0].fix(0.000484863)
    cell.fuel_inlet.mole_frac_comp[0, "H2"].fix(0.0630491)
    cell.fuel_inlet.mole_frac_comp[0, "H2O"].fix(0.6273812)
    cell.fuel_inlet.mole_frac_comp[0, "N2"].fix(0.3095697)

    cell.oxygen_inlet.temperature[0].fix(1103.15)
    cell.oxygen_inlet.pressure.fix(85910)

    cell.recursive_scaling()

    cell.initialize_build(
        optarg={"nlp_scaling_method": "user-scaling"},
        current_density_guess=-2000,
        temperature_guess=1103.15,
    )
    cell.model_check()
    # Test whether fixed degrees of freedom remain fixed
    assert degrees_of_freedom(m.fs.cell) == 0

    approx = lambda x: pytest.approx(x, 5e-3)
    assert cell.current_density[0, 1].value == approx(-2394.77)
    assert cell.current_density[0, 3].value == approx(-2326.71)
    assert cell.current_density[0, 5].value == approx(-2268.31)
    assert cell.current_density[0, 8].value == approx(-2191.66)
    assert cell.current_density[0, 10].value == approx(-2145.27)

    assert cell.fuel_outlet.temperature[0].value == approx(1103.40)
    assert pyo.value(cell.fuel_outlet.pressure[0]) == approx(85910)
    assert approx(484.863e-6) == pyo.value(cell.fuel_outlet.flow_mol[0])
    assert approx(0.472735) == pyo.value(cell.fuel_outlet.mole_frac_comp[0, "H2O"])
    assert approx(0.217695) == pyo.value(cell.fuel_outlet.mole_frac_comp[0, "H2"])
    assert approx(0.309570) == pyo.value(cell.fuel_outlet.mole_frac_comp[0, "N2"])

    assert cell.oxygen_outlet.temperature[0].value == approx(1103.42)
    assert pyo.value(cell.oxygen_outlet.pressure[0]) == approx(85910)
    assert approx(297.821e-6) == pyo.value(cell.oxygen_outlet.flow_mol[0])
    assert approx(0.3094487) == pyo.value(cell.oxygen_outlet.mole_frac_comp[0, "O2"])
    assert approx(0.690551) == pyo.value(cell.oxygen_outlet.mole_frac_comp[0, "N2"])

    # Test whether unfixed degrees of freedom remain unfixed
    cell.potential.unfix()
    cell.fuel_inlet.temperature[0].unfix()
    cell.fuel_inlet.pressure.unfix()
    cell.fuel_inlet.flow_mol[0].unfix()
    cell.fuel_inlet.mole_frac_comp[0, "H2"].unfix()
    cell.fuel_inlet.mole_frac_comp[0, "H2O"].unfix()
    cell.fuel_inlet.mole_frac_comp[0, "N2"].unfix()

    cell.oxygen_inlet.temperature[0].unfix()
    cell.oxygen_inlet.pressure.unfix()
    cell.oxygen_inlet.mole_frac_comp[0, "O2"].unfix()
    cell.oxygen_inlet.mole_frac_comp[0, "N2"].unfix()

    assert degrees_of_freedom(cell) == 11

    cell.initialize(current_density_guess=-1500, temperature_guess=1103.15)

    assert degrees_of_freedom(cell) == 11

    # Clean up side effects for other tests
    cell.oxygen_inlet.mole_frac_comp[0, "O2"].fix()
    cell.oxygen_inlet.mole_frac_comp[0, "N2"].fix()


@pytest.mark.component
def test_initialization_stack(model_stack):
    m = model_stack
    stack = m.fs.stack
    cell = m.fs.stack.solid_oxide_cell

    stack.potential_cell[0].fix(1.288)
    stack.fuel_inlet.temperature[0].fix(1103.15)
    stack.fuel_inlet.pressure.fix(85910)
    stack.fuel_inlet.flow_mol[0].fix(0.00484863)
    stack.fuel_inlet.mole_frac_comp[0, "H2"].fix(0.0630491)
    stack.fuel_inlet.mole_frac_comp[0, "H2O"].fix(0.6273812)
    stack.fuel_inlet.mole_frac_comp[0, "N2"].fix(0.3095697)

    stack.oxygen_inlet.temperature[0].fix(1103.15)
    stack.oxygen_inlet.pressure.fix(85910)

    iscale.calculate_scaling_factors(m)

    assert degrees_of_freedom(m.fs.stack) == 0

    stack.initialize_build(
        optarg={"nlp_scaling_method": "user-scaling"},
        current_density_guess=-2000,
        temperature_guess=1103.15,
    )
    cell.model_check()
    # Test whether fixed degrees of freedom remain fixed
    assert degrees_of_freedom(m.fs.stack) == 0

    approx = lambda x: pytest.approx(x, 5e-3)
    assert cell.current_density[0, 1].value == approx(-2394.77)
    assert cell.current_density[0, 3].value == approx(-2326.71)
    assert cell.current_density[0, 5].value == approx(-2268.31)
    assert cell.current_density[0, 8].value == approx(-2191.66)
    assert cell.current_density[0, 10].value == approx(-2145.27)

    assert stack.fuel_outlet.temperature[0].value == approx(1103.40)
    assert pyo.value(stack.fuel_outlet.pressure[0]) == approx(85910)
    assert approx(484.863e-5) == pyo.value(stack.fuel_outlet.flow_mol[0])
    assert approx(0.472735) == pyo.value(stack.fuel_outlet.mole_frac_comp[0, "H2O"])
    assert approx(0.217695) == pyo.value(stack.fuel_outlet.mole_frac_comp[0, "H2"])
    assert approx(0.309570) == pyo.value(stack.fuel_outlet.mole_frac_comp[0, "N2"])

    assert stack.oxygen_outlet.temperature[0].value == approx(1103.42)
    assert pyo.value(stack.oxygen_outlet.pressure[0]) == approx(85910)
    assert approx(297.821e-5) == pyo.value(stack.oxygen_outlet.flow_mol[0])
    assert approx(0.3094487) == pyo.value(stack.oxygen_outlet.mole_frac_comp[0, "O2"])
    assert approx(0.690551) == pyo.value(stack.oxygen_outlet.mole_frac_comp[0, "N2"])

    # Test whether unfixed degrees of freedom remain unfixed
    stack.potential_cell.unfix()
    stack.fuel_inlet.temperature[0].unfix()
    stack.fuel_inlet.pressure.unfix()
    stack.fuel_inlet.flow_mol[0].unfix()
    stack.fuel_inlet.mole_frac_comp[0, "H2"].unfix()
    stack.fuel_inlet.mole_frac_comp[0, "H2O"].unfix()
    stack.fuel_inlet.mole_frac_comp[0, "N2"].unfix()

    stack.oxygen_inlet.temperature[0].unfix()
    stack.oxygen_inlet.pressure.unfix()
    stack.oxygen_inlet.mole_frac_comp[0, "O2"].unfix()
    stack.oxygen_inlet.mole_frac_comp[0, "N2"].unfix()

    assert degrees_of_freedom(stack) == 11

    stack.initialize_build(
        optarg={"nlp_scaling_method": "user-scaling"},
        current_density_guess=-1500,
        temperature_guess=1103.15,
    )

    assert degrees_of_freedom(stack) == 11

    # Clean up side effects for other tests
    stack.oxygen_inlet.mole_frac_comp[0, "O2"].fix()
    stack.oxygen_inlet.mole_frac_comp[0, "N2"].fix()


def kazempoor_braun_replication(model):
    m = model
    cell = m.fs.cell
    case = 5

    solver = pyo.SolverFactory("ipopt")

    N_voltage = 11
    N_case = 5
    _constF = 96485
    df = pd.read_csv(os.sep.join([data_cache, "herring-et-al-data.csv"]), index_col=0)

    fuel_out = cell.fuel_outlet
    oxygen_out = cell.oxygen_outlet
    tracked_vars = {
        "voltage": cell.potential[0],
        "average_current_density": sum(cell.current_density[0, j] for j in cell.iznodes)
        / len(cell.iznodes),
        "average_temperature_z": sum(cell.temperature_z[0, j] for j in cell.iznodes)
        / len(cell.iznodes),
        "fuel_out_flow_mol": fuel_out.flow_mol[0],
        "fuel_out_temperature": fuel_out.temperature[0],
        "fuel_out_pressure": fuel_out.pressure[0],
        "fuel_out_mol_frac_comp[H2]": fuel_out.mole_frac_comp[0, "H2"],
        "fuel_out_mol_frac_comp[H2O]": fuel_out.mole_frac_comp[0, "H2O"],
        "fuel_out_mol_frac_comp[N2]": fuel_out.mole_frac_comp[0, "N2"],
        "oxygen_out_flow_mol": oxygen_out.flow_mol[0],
        "oxygen_out_temperature": oxygen_out.temperature[0],
        "oxygen_out_pressure": oxygen_out.pressure[0],
        "oxygen_out_mol_frac_comp[O2]": oxygen_out.mole_frac_comp[0, "O2"],
        "oxygen_out_mol_frac_comp[N2]": oxygen_out.mole_frac_comp[0, "N2"],
    }

    results = dict()
    for key in tracked_vars.keys():
        results[key] = []

    for case in range(1, N_case + 1):

        T_in = df["T_in"][case] + 273.15
        T_dew = df["T_dew"][case] + 273.15
        N_N2 = cccm_to_mps(df["sccm_N2"][case])
        N_H2 = cccm_to_mps(df["sccm_H2"][case])

        # IAPWS95 returns psat in kPa
        p_H2O = 1e3 * pyo.value(m.fs.prop_Iapws95.p_sat_func("H2O", 647.096 / T_dew))
        y_H2O = p_H2O / P
        N_H2O = (N_N2 + N_H2) * y_H2O / (1 - y_H2O)

        N_tot = N_N2 + N_H2 + N_H2O

        y_H2 = N_H2 / N_tot

        assert abs(y_H2O * N_tot - N_H2O) < 1e-15

        dHrxn = (
            _comp_enthalpy_expr(T_in, "H2")
            + 0.5 * _comp_enthalpy_expr(T_in, "O2")
            - _comp_enthalpy_expr(T_in, "H2O")
        )
        dSrxn = (
            _comp_entropy_expr(T_in, "H2")
            + 0.5 * _comp_entropy_expr(T_in, "O2")
            - _comp_entropy_expr(T_in, "H2O")
        )
        dGrxn = dHrxn - T_in * dSrxn

        V_nernst = pyo.value(
            1
            / (2 * _constF)
            * (
                dGrxn
                - 8.31
                * T_in
                * pyo.log(y_H2O / (y_H2 * 0.21**0.5) * (P / 101300) ** -0.5)
            )
        )

        m.fs.cell.potential.fix(V_nernst)
        m.fs.cell.fuel_inlet.temperature[0].fix(T_in)
        m.fs.cell.fuel_inlet.pressure.fix(P)
        m.fs.cell.fuel_inlet.flow_mol[0].fix(N_tot / N_cell)
        m.fs.cell.fuel_inlet.mole_frac_comp[0, "H2"].fix(y_H2)
        m.fs.cell.fuel_inlet.mole_frac_comp[0, "H2O"].fix(y_H2O)
        m.fs.cell.fuel_inlet.mole_frac_comp[0, "N2"].fix(N_N2 / N_tot)

        m.fs.cell.oxygen_inlet.temperature[0].fix(T_in)
        m.fs.cell.oxygen_inlet.pressure.fix(P)

        m.fs.cell.recursive_scaling()

        m.fs.cell.initialize_build(
            optarg={"nlp_scaling_method": "user-scaling"},
            current_density_guess=0,
            temperature_guess=T_in,
        )

        results["voltage"].append(
            np.linspace(V_nernst, 0.9 * df["V_max"][case] / N_cell, N_voltage)
        )
        for key, lst in results.items():
            if key == "voltage":
                continue
            else:
                lst.append([])

        for i in range(N_voltage):
            m.fs.cell.potential.fix(results["voltage"][case - 1][i])
            assert degrees_of_freedom(m) == 0
            status = solver.solve(
                m,
                tee=False,
                options={
                    "tol": 1e-6,
                    "max_iter": 100,
                    "nlp_scaling_method": "user-scaling",
                },
            )
            assert (
                status.solver.termination_condition == pyo.TerminationCondition.optimal
            )
            assert status.solver.status == pyo.SolverStatus.ok
            m.fs.cell.model_check()

            for key, lst in results.items():
                if key == "voltage":
                    continue
                else:
                    lst[case - 1].append(pyo.value(tracked_vars[key]))

    df = pd.DataFrame(results)
    out = []

    for row in df.iterrows():
        tmp_dict = {}
        for key, val in row[1].items():
            tmp_dict[key] = np.array(val)
        tmp_df = pd.DataFrame(tmp_dict)
        out.append(tmp_df)

    return out


@pytest.mark.integration
def test_model_replication(model):
    out = kazempoor_braun_replication(model)
    cached_results = []
    for i in range(len(out)):
        cached_results.append(
            pd.read_csv(
                os.sep.join([data_cache, f"case_{i+1}_interconnect.csv"]), index_col=0
            )
        )

    for df, cached_df in zip(out, cached_results):
        pd.testing.assert_frame_equal(
            df, cached_df, check_dtype=False, check_exact=False, rtol=3e-3
        )

    df = out[4]
    data = pd.read_csv(
        os.sep.join([data_cache, "sweep_5_kazempoor_replication.csv"]),
        names=["average_current_density", "voltage"],
    )
    intr = interp1d(
        x=data["voltage"].to_numpy(), y=data["average_current_density"].to_numpy()
    )
    for row in df.iterrows():
        assert intr(row[1]["voltage"]) == pytest.approx(
            -row[1]["average_current_density"], rel=3e-2, abs=50
        )


if __name__ == "__main__":
    m = model_func()
    out = kazempoor_braun_replication(m)
    # Uncomment to recreate cached data
    # for i, df in enumerate(out):
    #     df.to_csv(os.sep.join([data_cache, f"case_{i+1}_interconnect.csv"]))
