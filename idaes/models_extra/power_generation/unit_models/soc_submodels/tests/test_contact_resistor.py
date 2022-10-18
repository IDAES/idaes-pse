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

solver = pyo.SolverFactory("ipopt")


@pytest.fixture
def model():
    time_set = [0]
    zfaces = np.linspace(0, 1, 4).tolist()
    m = soc_testing._cell_flowsheet_model(
        dynamic=False, time_set=time_set, zfaces=zfaces
    )
    iznodes = m.fs.iznodes
    tset = m.fs.config.time
    m.fs.temperature_deviation_x = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.K
    )
    m.fs.heat_flux_x0 = pyo.Var(
        tset, iznodes, initialize=0, units=pyo.units.W / pyo.units.m**2
    )

    m.fs.contact = soc.SocContactResistor(
        control_volume_zfaces=zfaces,
        length_z=m.fs.length_z,
        length_y=m.fs.length_y,
        current_density=m.fs.current_density,
        temperature_z=m.fs.temperature_z,
        temperature_deviation_x=m.fs.temperature_deviation_x,
        heat_flux_x0=m.fs.heat_flux_x0,
    )
    m.fs.temperature_deviation_x.fix(0)
    m.fs.heat_flux_x0.fix(0)

    m.fs.contact.log_preexponential_factor.fix(pyo.log(0.46e-4))
    m.fs.contact.thermal_exponent_dividend.fix(0)
    m.fs.contact.contact_fraction.fix(1)

    return m


@pytest.fixture
def model2():
    time_set = [0, 1]
    zfaces = np.linspace(0, 1, 8).tolist()
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False, time_set=time_set, time_units=pyo.units.s)
    m.fs.contact = soc.SocContactResistor(control_volume_zfaces=zfaces)
    m.fs.contact.current_density.fix(0)
    m.fs.contact.temperature_z.fix(0)
    m.fs.contact.temperature_deviation_x.fix(0)
    m.fs.contact.heat_flux_x0.fix(0)

    m.fs.contact.log_preexponential_factor.fix(pyo.log(0.46e-4))
    m.fs.contact.thermal_exponent_dividend.fix(0)
    m.fs.contact.contact_fraction.fix(1)

    return m


@pytest.mark.build
@pytest.mark.unit
def test_build(model):
    contact = model.fs.contact
    nz = len(contact.iznodes)
    nt = len(model.fs.time)
    soc_testing._build_test_utility(
        contact,
        comp_dict={
            pyo.Var: {
                "temperature_z": nz * nt,
                "temperature_deviation_x": nz * nt,
                "heat_flux_x0": nz * nt,
                "current_density": nz * nt,
                "length_z": 1,
                "length_y": 1,
                "heat_flux_x1": nz * nt,
                "log_preexponential_factor": 1,
                "thermal_exponent_dividend": 1,
                "contact_fraction": 1,
            },
            pyo.Constraint: {"heat_flux_x_eqn": nz * nt},
            pyo.Expression: {
                "temperature": nz * nt,
                "contact_resistance": nz * nt,
                "voltage_drop_total": nz * nt,
                "joule_heating_flux": nz * nt,
            },
        },
        references=[
            "temperature_z",
            "temperature_deviation_x",
            "heat_flux_x0",
            "current_density",
            "length_z",
            "length_y",
        ],
    )
    assert degrees_of_freedom(model.fs.contact) == 0


@pytest.mark.build
@pytest.mark.unit
def test_build2(model2):
    contact = model2.fs.contact
    nz = len(contact.iznodes)
    nt = len(model2.fs.time)
    soc_testing._build_test_utility(
        contact,
        comp_dict={
            pyo.Var: {
                "temperature_z": nz * nt,
                "temperature_deviation_x": nz * nt,
                "heat_flux_x0": nz * nt,
                "current_density": nz * nt,
                "length_z": 1,
                "length_y": 1,
                "heat_flux_x1": nz * nt,
                "log_preexponential_factor": 1,
                "thermal_exponent_dividend": 1,
                "contact_fraction": 1,
            },
            pyo.Constraint: {"heat_flux_x_eqn": nz * nt},
            pyo.Expression: {
                "temperature": nz * nt,
                "contact_resistance": nz * nt,
                "voltage_drop_total": nz * nt,
                "joule_heating_flux": nz * nt,
            },
        },
    )
    assert degrees_of_freedom(model2.fs.contact) == 0


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialization(model):
    model.fs.contact.initialize_build(optarg={"nlp_scaling_method": "user-scaling"})
    model.fs.heat_flux_x0.unfix()
    model.fs.contact.heat_flux_x1.fix()
    model.fs.contact.initialize_build(
        fix_heat_flux_x0=False, optarg={"nlp_scaling_method": "user-scaling"}
    )
