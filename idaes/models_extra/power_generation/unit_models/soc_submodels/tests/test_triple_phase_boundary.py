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
            "log_mole_frac_comp": nz * nt * ncomp,
            "activation_potential": nz * nt,
            "activation_potential_alpha1": 1,
            "activation_potential_alpha2": 1,
            "exchange_current_exponent_comp": nreact,
            "exchange_current_log_preexponential_factor": 1,
            "exchange_current_activation_energy": 1,
        },
        pyo.Constraint: {
            "mole_frac_comp_eqn": nz * nt * ncomp,
            "log_mole_frac_comp_eqn": nz * nt * ncomp,
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
            "nernst_potential": nz * nt,
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
    m.fs.fuel_tpb = soc.SocTriplePhaseBoundary(
        default={
            "control_volume_zfaces": zfaces,
            "length_z": m.fs.length_z,
            "length_y": m.fs.length_y,
            "component_list": fuel_comps,
            "tpb_stoich_dict": {
                "H2": -0.5,
                "H2O": 0.5,
                "N2": 0,
                "Vac": 0.5,
                "O^2-": -0.5,
            },
            "inert_species": ["N2"],
            "current_density": m.fs.current_density,
            "temperature_z": m.fs.temperature_z,
            "temperature_deviation_x": m.fs.temperature_deviation_x,
            "heat_flux_x0": m.fs.heat_flux_x0,
            "conc_mol_comp_ref": m.fs.conc_mol_comp_ref,
            "conc_mol_comp_deviation_x": m.fs.conc_mol_comp_deviation_x,
            "material_flux_x": m.fs.material_flux_x,
        }
    )
    m.fs.temperature_deviation_x.fix(0)
    m.fs.heat_flux_x0.fix(0)
    m.fs.conc_mol_comp_ref.fix(1)
    m.fs.conc_mol_comp_deviation_x.fix(0)

    m.fs.fuel_tpb.exchange_current_log_preexponential_factor.fix(pyo.log(1.375e10))
    m.fs.fuel_tpb.exchange_current_activation_energy.fix(120e3)
    m.fs.fuel_tpb.activation_potential_alpha1.fix(0.5)
    m.fs.fuel_tpb.activation_potential_alpha2.fix(0.5)

    m.fs.fuel_tpb.exchange_current_exponent_comp["H2"].fix(1)
    m.fs.fuel_tpb.exchange_current_exponent_comp["H2O"].fix(1)

    return m


@pytest.fixture
def modelOxygen():
    time_set = [0, 1, 2]
    zfaces = np.linspace(0, 1, 8).tolist()
    o2_comps = ["O2", "N2"]
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": False,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    m.fs.oxygen_tpb = soc.SocTriplePhaseBoundary(
        default={
            "control_volume_zfaces": zfaces,
            "component_list": o2_comps,
            "tpb_stoich_dict": {"O2": -0.25, "Vac": -0.5, "O^2-": 0.5},
            "inert_species": ["N2"],
        }
    )
    m.fs.oxygen_tpb.temperature_z.fix(0)
    m.fs.oxygen_tpb.current_density.fix(0)

    m.fs.oxygen_tpb.temperature_deviation_x.fix(0)
    m.fs.oxygen_tpb.heat_flux_x0.fix(0)
    m.fs.oxygen_tpb.conc_mol_comp_ref.fix(1)
    m.fs.oxygen_tpb.conc_mol_comp_deviation_x.fix(0)

    m.fs.oxygen_tpb.exchange_current_log_preexponential_factor.fix(pyo.log(1.375e10))
    m.fs.oxygen_tpb.exchange_current_activation_energy.fix(120e3)
    m.fs.oxygen_tpb.activation_potential_alpha1.fix(0.5)
    m.fs.oxygen_tpb.activation_potential_alpha2.fix(0.5)

    m.fs.oxygen_tpb.exchange_current_exponent_comp["O2"].fix(0.25)

    return m


@pytest.mark.build
@pytest.mark.unit
def test_build_fuel(modelFuel):
    tpb = modelFuel.fs.fuel_tpb
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
    tpb = modelOxygen.fs.oxygen_tpb
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
    modelFuel.fs.fuel_tpb.initialize(
        fix_x0=True, optarg={"nlp_scaling_method": "user-scaling"}
    )


@pytest.mark.unit
def test_extra_inert():
    time_set = [0, 1, 2]
    zfaces = np.linspace(0, 1, 8).tolist()
    o2_comps = ["O2", "N2"]
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": False,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    with pytest.raises(
        ConfigurationError,
        match="fs.oxygen_tpb invalid component in inert_species "
        "argument. H2O is not in the provided component list.",
    ):
        m.fs.oxygen_tpb = soc.SocTriplePhaseBoundary(
            default={
                "control_volume_zfaces": zfaces,
                "component_list": o2_comps,
                "tpb_stoich_dict": {"O2": -0.25, "Vac": -0.5, "O^2-": 0.5},
                "inert_species": ["N2", "H2O"],
            }
        )
