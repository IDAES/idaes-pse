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
from idaes.core.util.constants import Constants
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError
import idaes.models_extra.power_generation.unit_models.soc_submodels as soc
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
import idaes.models_extra.power_generation.unit_models.soc_submodels.testing as soc_testing

solver = pyo.SolverFactory("ipopt")


def common_components(nt, nz, ncomp, nreact):
    return {
        pyo.Var: {
            "temperature_z": nz * nt,
            "temperature_deviation_x": nz * nt,
            "conc_mol_comp_ref": nz * nt * ncomp,
            "conc_mol_comp_deviation_x": nz * nt * ncomp,
            "material_flux_x": nz * nt * ncomp,
            "heat_flux_x0": nz * nt,
            "current_density": nz * nt,
            "length_z": 1,
            "length_y": 1,
            "heat_flux_x1": nz * nt,
            "mole_frac_comp": nz * nt * ncomp,
            "log_mole_frac_comp": nz * nt * nreact,
            "activation_potential": nz * nt,
            "activation_potential_alpha1": 1,
            "activation_potential_alpha2": 1,
            "exchange_current_exponent_comp": nreact,
            "exchange_current_log_preexponential_factor": 1,
            "exchange_current_activation_energy": 1,
        },
        pyo.Constraint: {
            "mole_frac_comp_eqn": nz * nt * ncomp,
            "log_mole_frac_comp_eqn": nz * nt * nreact,
            "material_flux_x_eqn": nz * nt * ncomp,
            "activation_potential_eqn": nz * nt,
            "heat_flux_x_eqn": nz * nt,
        },
        pyo.Expression: {
            "temperature": nz * nt,
            "conc_mol_comp": nz * nt * ncomp,
            "pressure": nz * nt,
            "ds_rxn": nz * nt,
            "dh_rxn": nz * nt,
            "dg_rxn": nz * nt,
            "potential_nernst": nz * nt,
            "log_exchange_current_density": nz * nt,
            "reaction_rate_per_unit_area": nz * nt,
            "voltage_drop_total": nz * nt,
        },
    }


@pytest.fixture
def modelFuel():
    time_set = [0]
    zfaces = np.linspace(0, 1, 6).tolist()
    fuel_comps = ["H2", "N2", "H2O"]
    m = soc_testing._cell_flowsheet_model(
        dynamic=False, time_set=time_set, zfaces=zfaces
    )
    iznodes = m.fs.iznodes
    # time_units = m.fs.time_units
    tset = m.fs.config.time
    comps = m.fs.comps = pyo.Set(initialize=fuel_comps)

    m.fs.temperature_deviation_x = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.K
    )
    m.fs.conc_mol_comp_ref = pyo.Var(
        tset, iznodes, comps, initialize=1, units=pyo.units.mol / pyo.units.m**3
    )
    m.fs.conc_mol_comp_deviation_x = pyo.Var(
        tset, iznodes, comps, initialize=0, units=pyo.units.mol / pyo.units.m**3
    )
    m.fs.material_flux_x = pyo.Var(
        tset,
        iznodes,
        comps,
        initialize=0,
        units=pyo.units.mol / (pyo.units.s * pyo.units.m**2),
    )
    m.fs.heat_flux_x0 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.W / pyo.units.m**2
    )
    m.fs.fuel_triple_phase_boundary = soc.SocTriplePhaseBoundary(
        control_volume_zfaces=zfaces,
        length_z=m.fs.length_z,
        length_y=m.fs.length_y,
        component_list=fuel_comps,
        reaction_stoichiometry={"H2": -0.5, "H2O": 0.5, "N2": 0, "e^-": 1.0},
        inert_species=["N2"],
        below_electrolyte=True,
        current_density=m.fs.current_density,
        temperature_z=m.fs.temperature_z,
        temperature_deviation_x=m.fs.temperature_deviation_x,
        heat_flux_x0=m.fs.heat_flux_x0,
        conc_mol_comp_ref=m.fs.conc_mol_comp_ref,
        conc_mol_comp_deviation_x=m.fs.conc_mol_comp_deviation_x,
        material_flux_x=m.fs.material_flux_x,
    )
    m.fs.temperature_deviation_x.fix(0)
    m.fs.heat_flux_x0.fix(0)
    m.fs.conc_mol_comp_ref.fix(1)
    m.fs.conc_mol_comp_deviation_x.fix(0)

    m.fs.fuel_triple_phase_boundary.exchange_current_log_preexponential_factor.fix(
        pyo.log(1.375e10)
    )
    m.fs.fuel_triple_phase_boundary.exchange_current_activation_energy.fix(120e3)
    m.fs.fuel_triple_phase_boundary.activation_potential_alpha1.fix(0.5)
    m.fs.fuel_triple_phase_boundary.activation_potential_alpha2.fix(0.5)

    m.fs.fuel_triple_phase_boundary.exchange_current_exponent_comp["H2"].fix(1)
    m.fs.fuel_triple_phase_boundary.exchange_current_exponent_comp["H2O"].fix(1)

    return m


@pytest.fixture
def modelOxygen():
    time_set = [0, 1, 2]
    zfaces = np.linspace(0, 1, 8).tolist()
    o2_comps = ["O2", "N2"]
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=time_set, time_units=pyo.units.s)
    m.fs.oxygen_triple_phase_boundary = soc.SocTriplePhaseBoundary(
        control_volume_zfaces=zfaces,
        component_list=o2_comps,
        reaction_stoichiometry={"O2": -0.25, "e^-": -1.0},
        inert_species=["N2"],
    )
    m.fs.oxygen_triple_phase_boundary.temperature_z.fix(0)
    m.fs.oxygen_triple_phase_boundary.current_density.fix(0)

    m.fs.oxygen_triple_phase_boundary.temperature_deviation_x.fix(0)
    m.fs.oxygen_triple_phase_boundary.heat_flux_x0.fix(0)
    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_ref.fix(1)
    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.fix(0)

    m.fs.oxygen_triple_phase_boundary.exchange_current_log_preexponential_factor.fix(
        pyo.log(1.375e10)
    )
    m.fs.oxygen_triple_phase_boundary.exchange_current_activation_energy.fix(120e3)
    m.fs.oxygen_triple_phase_boundary.activation_potential_alpha1.fix(0.5)
    m.fs.oxygen_triple_phase_boundary.activation_potential_alpha2.fix(0.5)

    m.fs.oxygen_triple_phase_boundary.exchange_current_exponent_comp["O2"].fix(0.25)

    return m


@pytest.mark.build
@pytest.mark.unit
def test_build_fuel(modelFuel):
    tpb = modelFuel.fs.fuel_triple_phase_boundary
    nz = len(modelFuel.fs.iznodes)
    nt = len(modelFuel.fs.time)
    ncomp = len(tpb.component_list)
    nreact = 2
    soc_testing._build_test_utility(
        block=tpb,
        comp_dict=common_components(nt, nz, ncomp, nreact),
        references=[
            "temperature_z",
            "temperature_deviation_x",
            "conc_mol_comp_ref",
            "conc_mol_comp_deviation_x",
            "material_flux_x",
            "heat_flux_x0",
            "current_density",
            "length_z",
            "length_y",
        ],
    )
    assert degrees_of_freedom(tpb) == 0


@pytest.mark.build
@pytest.mark.unit
def test_build_oxygen(modelOxygen):
    tpb = modelOxygen.fs.oxygen_triple_phase_boundary
    nz = len(tpb.iznodes)
    nt = len(modelOxygen.fs.time)
    ncomp = len(tpb.component_list)
    nreact = 1
    soc_testing._build_test_utility(
        block=tpb,
        comp_dict=common_components(nt, nz, ncomp, nreact),
    )

    assert degrees_of_freedom(tpb) == 0


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialization_fuel(modelFuel):
    modelFuel.fs.fuel_triple_phase_boundary.initialize(
        fix_x0=True, optarg={"nlp_scaling_method": "user-scaling"}
    )


@pytest.mark.unit
def test_extra_inert():
    time_set = [0, 1, 2]
    zfaces = np.linspace(0, 1, 8).tolist()
    o2_comps = ["O2", "N2"]
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=time_set, time_units=pyo.units.s)
    with pytest.raises(
        ConfigurationError,
        match="fs.oxygen_triple_phase_boundary invalid component in inert_species "
        "argument. H2O is not in the provided component list.",
    ):
        m.fs.oxygen_triple_phase_boundary = soc.SocTriplePhaseBoundary(
            control_volume_zfaces=zfaces,
            component_list=o2_comps,
            reaction_stoichiometry={"O2": -0.25, "Vac": -0.5, "O^2-": 0.5, "e^-": -1.0},
            inert_species=["N2", "H2O"],
        )


def modelFuelAndOxygen(include_solid_species):
    time_set = [0]
    zfaces = [0, 1]
    fuel_comps = ["H2", "N2", "H2O"]
    oxygen_comps = ["O2", "N2"]
    fuel_stoich = {
        "H2": -1.0,
        "H2O": 1.0,
        "N2": 0.0,
        "e^-": 2.0,
    }
    oxygen_stoich = {"O2": -1.0, "e^-": -4.0}
    if include_solid_species:
        fuel_stoich["Vac"] = 1.0
        fuel_stoich["O^2-"] = -1.0
        oxygen_stoich["Vac"] = -2.0
        oxygen_stoich["O^2-"] = 2.0

    m = soc_testing._cell_flowsheet_model(
        dynamic=False, time_set=time_set, zfaces=zfaces
    )
    iznodes = m.fs.iznodes
    # time_units = m.fs.time_units
    tset = m.fs.config.time
    comps = m.fs.comps = pyo.Set(initialize=fuel_comps)

    m.fs.temperature_deviation_x = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.K
    )
    m.fs.fuel_triple_phase_boundary = soc.SocTriplePhaseBoundary(
        control_volume_zfaces=zfaces,
        length_z=m.fs.length_z,
        length_y=m.fs.length_y,
        component_list=fuel_comps,
        reaction_stoichiometry=fuel_stoich,
        inert_species=["N2"],
        below_electrolyte=True,
        current_density=m.fs.current_density,
        temperature_z=m.fs.temperature_z,
        temperature_deviation_x=m.fs.temperature_deviation_x,
    )
    m.fs.oxygen_triple_phase_boundary = soc.SocTriplePhaseBoundary(
        control_volume_zfaces=zfaces,
        length_z=m.fs.length_z,
        length_y=m.fs.length_y,
        component_list=oxygen_comps,
        reaction_stoichiometry=oxygen_stoich,
        inert_species=["N2"],
        below_electrolyte=False,
        current_density=m.fs.current_density,
        temperature_z=m.fs.temperature_z,
        temperature_deviation_x=m.fs.temperature_deviation_x,
        heat_flux_x0=m.fs.fuel_triple_phase_boundary.heat_flux_x1,
    )

    m.fs.fuel_triple_phase_boundary.exchange_current_log_preexponential_factor.fix(
        pyo.log(1.375e10)
    )
    m.fs.fuel_triple_phase_boundary.exchange_current_activation_energy.fix(120e3)
    m.fs.fuel_triple_phase_boundary.activation_potential_alpha1.fix(0.5)
    m.fs.fuel_triple_phase_boundary.activation_potential_alpha2.fix(0.5)

    m.fs.fuel_triple_phase_boundary.exchange_current_exponent_comp["H2"].fix(1)
    m.fs.fuel_triple_phase_boundary.exchange_current_exponent_comp["H2O"].fix(1)

    m.fs.oxygen_triple_phase_boundary.exchange_current_log_preexponential_factor.fix(
        pyo.log(1.375e10)
    )
    m.fs.oxygen_triple_phase_boundary.exchange_current_activation_energy.fix(120e3)
    m.fs.oxygen_triple_phase_boundary.activation_potential_alpha1.fix(0.5)
    m.fs.oxygen_triple_phase_boundary.activation_potential_alpha2.fix(0.5)

    T = 1000
    m.fs.temperature_z.fix(T)
    m.fs.temperature_deviation_x.fix(0)
    m.fs.fuel_triple_phase_boundary.conc_mol_comp_deviation_x.fix(0)
    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_deviation_x.fix(0)

    m.fs.oxygen_triple_phase_boundary.exchange_current_exponent_comp["O2"].fix(0.25)

    C_tot = 1.2e5 / pyo.value(Constants.gas_constant * T)

    m.fs.fuel_triple_phase_boundary.conc_mol_comp_ref[0, :, "H2O"].fix(0.5 * C_tot)
    m.fs.fuel_triple_phase_boundary.conc_mol_comp_ref[0, :, "H2"].fix(0.05 * C_tot)
    m.fs.fuel_triple_phase_boundary.conc_mol_comp_ref[0, :, "N2"].fix(0.45 * C_tot)

    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_ref[0, :, "O2"].fix(0.21 * C_tot)
    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_ref[0, :, "N2"].fix(0.79 * C_tot)

    m.fs.fuel_triple_phase_boundary.heat_flux_x1.fix(0)
    m.fs.current_density.fix(1000)

    @m.fs.Expression(m.fs.time, m.fs.fuel_triple_phase_boundary.iznodes)
    def net_energy_flux_out(b, t, iz):
        return (
            -sum(
                b.fuel_triple_phase_boundary.material_flux_x[t, iz, comp]
                * common._comp_enthalpy_expr(
                    b.fuel_triple_phase_boundary.temperature[t, iz], comp
                )
                for comp in b.fuel_triple_phase_boundary.component_list
            )
            - b.fuel_triple_phase_boundary.heat_flux_x0[t, iz]
            + sum(
                b.oxygen_triple_phase_boundary.material_flux_x[t, iz, comp]
                * common._comp_enthalpy_expr(
                    b.oxygen_triple_phase_boundary.temperature[t, iz], comp
                )
                for comp in b.oxygen_triple_phase_boundary.component_list
            )
            + b.oxygen_triple_phase_boundary.heat_flux_x1[t, iz]
        )

    @m.fs.Expression(m.fs.time, m.fs.fuel_triple_phase_boundary.iznodes)
    def voltage_difference(b, t, iz):
        return (
            b.fuel_triple_phase_boundary.potential_nernst[t, iz]
            + b.oxygen_triple_phase_boundary.potential_nernst[t, iz]
            - (
                b.fuel_triple_phase_boundary.voltage_drop_total[t, iz]
                + b.oxygen_triple_phase_boundary.voltage_drop_total[t, iz]
            )
        )

    @m.fs.Expression(m.fs.time, m.fs.fuel_triple_phase_boundary.iznodes)
    def electric_work_flux(b, t, iz):
        return b.voltage_difference[t, iz] * b.current_density[t, iz]

    @m.fs.Expression(m.fs.time, m.fs.fuel_triple_phase_boundary.iznodes)
    def energy_created(b, t, iz):
        return b.net_energy_flux_out[t, iz] + b.electric_work_flux[t, iz]

    m.fs.element_list = pyo.Set(initialize=["H", "O", "N"], ordered=True)

    @m.fs.Expression(
        m.fs.time, m.fs.fuel_triple_phase_boundary.iznodes, m.fs.element_list
    )
    def element_created(b, t, iz, i):
        return -sum(
            common._element_dict[i][j]
            * b.fuel_triple_phase_boundary.material_flux_x[t, iz, j]
            for j in b.fuel_triple_phase_boundary.component_list
        ) + sum(
            common._element_dict[i][j]
            * b.oxygen_triple_phase_boundary.material_flux_x[t, iz, j]
            for j in b.oxygen_triple_phase_boundary.component_list
        )

    return m


def conservation_tester(m):
    for P_fuel in np.linspace(1e5, 5e5, 3):
        for T in np.linspace(900, 1100, 4):
            m.fs.temperature_z.fix(T)
            C_tot = P_fuel / pyo.value(Constants.gas_constant * T)
            for y_H2 in np.linspace(0.1, 0.5, 3):
                m.fs.fuel_triple_phase_boundary.conc_mol_comp_ref[0, :, "H2"].fix(
                    y_H2 * C_tot
                )
                m.fs.fuel_triple_phase_boundary.conc_mol_comp_ref[0, :, "N2"].fix(
                    0.4 * C_tot
                )
                m.fs.fuel_triple_phase_boundary.conc_mol_comp_ref[0, :, "H2O"].fix(
                    (0.6 - y_H2) * C_tot
                )
                for y_O2 in np.linspace(0.1, 0.3, 3):
                    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_ref[0, :, "O2"].fix(
                        y_O2 * C_tot
                    )
                    m.fs.oxygen_triple_phase_boundary.conc_mol_comp_ref[0, :, "N2"].fix(
                        (1 - y_O2) * C_tot
                    )
                    for J in np.linspace(-1000, 1000, 5):
                        m.fs.current_density.fix(J)
                        solver.solve(m)
                        assert pyo.value(m.fs.energy_created[0, 1]) == pytest.approx(
                            0, rel=1e-4
                        )
                        for i in m.fs.element_list:
                            assert pyo.value(
                                m.fs.element_created[0, 1, i]
                            ) == pytest.approx(0, rel=1e-4)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.integration
def test_conservation_no_solid_species():
    m = modelFuelAndOxygen(False)
    solver.solve(m)
    conservation_tester(m)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.integration
def test_conservation_solid_species():
    m = modelFuelAndOxygen(True)
    solver.solve(m)
    conservation_tester(m)
