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
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
import idaes.models_extra.power_generation.unit_models.soc_submodels.testing as soc_testing

def common_components(nt, nz, ncomp):
    return {
        pyo.Var: {
            "temperature_z": nz * nt,
            "Dtemp_x0": nz * nt,
            "Dtemp": nz * nt,
            "Dtemp_x1": nz * nt,
            "enth_mol": nz * nt,
            "heat_transfer_coefficient": nt * nz,
            "velocity": nz * nt,
            "pressure": nz * nt,
            "flow_mol": nz * nt,
            "mole_frac_comp": nz * nt * ncomp,
            "conc": nz * nt * ncomp,
            "qflux_x0": nz * nt,
            "qflux_x1": nz * nt,
            "length_x": 1,
            "length_z": 1,
            "length_y": 1,
            "flow_mol_inlet": nt,
            "pressure_inlet": nt,
            "temperature_inlet": nt,
            "temperature_outlet": nt,
            "mole_frac_comp_inlet": nt * ncomp,
        },
        pyo.Constraint: {
            "flow_mol_eqn": nz * nt,
            "constant_pressure_eqn": nz * nt,
            "conc_eqn": nz * nt * ncomp,
            "enth_mol_eqn": nt * nz,
            "mole_frac_eqn": nt * nz,
            "material_balance_eqn": nt * nz * ncomp,
            "energy_balance_eqn": nt * nz,
            "temperature_x0_eqn": nt * nz,
            "temperature_x1_eqn": nt * nz,
            "temperature_outlet_eqn": nt,
        },
        pyo.Expression: {
            "temperature_x0": nz * nt,
            "temperature": nz * nt,
            "temperature_x1": nz * nt,
            "flow_area": 1,
            "flow_mol_comp_inlet": nt * ncomp,
            "dz": nz,
            "node_volume": nz,
            "xface_area": nz,
            "enth_mol_inlet": nt,
            "diff_eff_coeff": nt * nz * ncomp,
            "volume_molar": nt * nz,
            "volume_molar_inlet": nt,
            "mass_transfer_coeff": nt * nz * ncomp,
            "zflux_inlet": nt * ncomp,
            "zflux_enth_inlet": nt,
            "zflux": nt * (nz + 1) * ncomp,
            "zflux_enth": nt * (nz + 1),
            "pressure_face": nt * (nz + 1),
            "material_balance_eqn": nt * nz * ncomp,
            "flow_mol_comp_outlet": nt * ncomp,
            "enth_mol_outlet": nt,
            # Pyomo considers these VarLikeExpressions as Expressions
            "flow_mol_outlet": nt,
            "pressure_outlet": nt,
            "mole_frac_comp_outlet": nt * ncomp,
        },
    }

def fix_boundary_conditions(channel):
    channel.qflux_x0.fix(0)
    channel.qflux_x1.fix(0)
    if channel.config.below_electrode:
        channel.xflux_x1.fix(0)
    else:
        channel.xflux_x0.fix(0)
    channel.flow_mol_inlet.fix()
    channel.temperature_inlet.fix()
    channel.pressure_inlet.fix()
    channel.mole_frac_comp_inlet.fix()

@pytest.fixture
def modelNoHoldup():
    time_set = [0, 1]
    zfaces = np.linspace(0, 1, 4).tolist()
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(
        default={
            "dynamic": False,
            "time_set": time_set,
            "time_units": pyo.units.s,
        }
    )
    m.fs.oxygen_chan = soc.SocChannel(
        default={
            "has_holdup": False,
            "cv_zfaces": zfaces,
            "interpolation_scheme": common.CV_Interpolation.UDS,
            "opposite_flow": True,
            "below_electrode": False,
            "comp_list": ["O2","H2O"],
        }
    )
    m.fs.oxygen_chan.length_x.fix(0.002)
    m.fs.oxygen_chan.heat_transfer_coefficient.fix(100)
    m.fs.oxygen_chan.length_y.fix(0.08)
    m.fs.oxygen_chan.length_z.fix(0.08)
    m.fs.oxygen_chan.temperature_z.fix(1000)
    fix_boundary_conditions(m.fs.oxygen_chan)
    return m

@pytest.fixture
def modelHoldupNotDynamic():
    time_set = [0, 1]
    zfaces = np.linspace(0, 1, 6).tolist()
    m = soc_testing._cell_flowsheet_model(
        dynamic=False,
        time_set=time_set,
        zfaces=zfaces
    )

    m.fs.fuel_chan = soc.SocChannel(
        default={
            "has_holdup": True,
            "cv_zfaces": zfaces,
            "length_z": m.fs.length_z,
            "length_y": m.fs.length_y,
            "temperature_z": m.fs.temperature_z,
            "interpolation_scheme": common.CV_Interpolation.UDS,
            "opposite_flow": False,
            "below_electrode": True,
            "comp_list": ["H2", "H2O", "N2"],
        }
    )
    m.fs.fuel_chan.length_x.fix(0.002)
    m.fs.fuel_chan.heat_transfer_coefficient.fix(100)
    fix_boundary_conditions(m.fs.fuel_chan)
    return m

@pytest.mark.build
@pytest.mark.unit
def test_build_modelNoHoldup(modelNoHoldup):
    channel = modelNoHoldup.fs.oxygen_chan
    nz = len(channel.znodes)
    nt = len(channel.flowsheet().time)
    ncomp = len(channel.comps)

    comp_dict = common_components(nt, nz, ncomp)
    comp_dict[pyo.Var]["Dconc_x0"] = nt * nz * ncomp
    comp_dict[pyo.Var]["xflux_x0"] = nt * nz * ncomp
    comp_dict[pyo.Constraint]["xflux_x0_eqn"] = nt * nz * ncomp
    comp_dict[pyo.Expression]["Dconc_x1"] = nt * nz * ncomp
    comp_dict[pyo.Expression]["xflux_x1"] = nt * nz * ncomp
    comp_dict[pyo.Param] = {}
    comp_dict[pyo.Param]["dcdt"] = nt * nz * ncomp
    comp_dict[pyo.Param]["dcedt"] = nt * nz

    soc_testing._build_test_utility(
        block=channel,
        comp_dict=comp_dict,
    )
    assert degrees_of_freedom(channel) == 0

@pytest.mark.build
@pytest.mark.unit
def test_build_modelHoldupNotDynamic(modelHoldupNotDynamic):
    channel = modelHoldupNotDynamic.fs.fuel_chan
    nz = len(channel.znodes)
    nt = len(channel.flowsheet().time)
    ncomp = len(channel.comps)

    comp_dict = common_components(nt, nz, ncomp)
    comp_dict[pyo.Expression]["Dconc_x0"] = nt * nz * ncomp
    comp_dict[pyo.Expression]["xflux_x0"] = nt * nz * ncomp
    comp_dict[pyo.Var]["Dconc_x1"] = nt * nz * ncomp
    comp_dict[pyo.Var]["xflux_x1"] = nt * nz * ncomp
    comp_dict[pyo.Constraint]["xflux_x1_eqn"] = nt * nz * ncomp
    comp_dict[pyo.Var]["int_energy_mol"] = nt * nz
    comp_dict[pyo.Var]["int_energy_density"] = nt * nz
    comp_dict[pyo.Param] = {}
    comp_dict[pyo.Param]["dcdt"] = nt * nz * ncomp
    comp_dict[pyo.Param]["dcedt"] = nt * nz
    comp_dict[pyo.Constraint]["int_energy_mol_eqn"] = nt * nz
    comp_dict[pyo.Constraint]["int_energy_density_eqn"] = nt * nz

    soc_testing._build_test_utility(
        block=channel,
        comp_dict=comp_dict,
        references=[
            "temperature_z",
            "length_z",
            "length_y",
        ],
    )
    assert degrees_of_freedom(channel) == 0