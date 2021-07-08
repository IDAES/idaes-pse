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

__author__ = "John Eslick"

import pytest
import pyomo.environ as pyo
from idaes.power_generation.flowsheets.gas_turbine.gas_turbine import (
    main, run_full_load)
from idaes.core.util.model_statistics import degrees_of_freedom

solver_available = pyo.SolverFactory('ipopt').available()

@pytest.mark.component
def test_build():
    comps = { # components present
        "CH4", "C2H6", "C2H4", "CO", "H2S", "H2", "O2", "H2O", "CO2", "N2",
        "Ar", "SO2"}
    rxns = { # reactions and key compoents for conversion
        "ch4_cmb":"CH4",
        "c2h6_cmb":"C2H6",
        "c2h4_cmb":"C2H4",
        "co_cmb":"CO",
        "h2s_cmb":"H2S",
        "h2_cmb":"H2"}
    phases = ["Vap"]
    air_comp = {
        "CH4":0.0,
        "C2H6":0.0,
        "C2H4":0.0,
        "CO":0.0,
        "H2S":0.0,
        "H2":0.0,
        "O2":0.2074,
        "H2O":0.0099,
        "CO2":0.0003,
        "N2":0.7732,
        "Ar":0.0092,
        "SO2":0.0}
    ng_comp = {
        "CH4":0.87,
        "C2H6":0.0846,
        "C2H4":0.0003,
        "CO":0.0009,
        "H2S":0.0004,
        "H2":0.0036,
        "O2":0.0007,
        "H2O":0.0,
        "CO2":0.0034,
        "N2":0.0361,
        "Ar":0.0,
        "SO2":0.0}

    m, solver = main(
        comps=comps,
        rxns=rxns,
        phases=phases,
        air_comp=air_comp,
        ng_comp=ng_comp,
        initialize=False)
    assert degrees_of_freedom(m) == 0


@pytest.mark.integration
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_initialize():
    # Speed this up by using a simpler composition and reaction set
    comps = {"CH4", "O2", "H2O", "CO2", "N2", "Ar"}
    rxns = {"ch4_cmb":"CH4"}
    phases = ["Vap"]
    air_comp = {
        "CH4":0.0,
        "O2":0.2074,
        "H2O":0.0099,
        "CO2":0.0003,
        "N2":0.7732,
        "Ar":0.0092}
    ng_comp = { # simplified composition to make it run faster
        "CH4":1.0,
        "O2":0.0,
        "H2O":0.0,
        "CO2":0.0,
        "N2":0.0,
        "Ar":0.0}

    m, solver = main(
        comps=comps,
        rxns=rxns,
        phases=phases,
        air_comp=air_comp,
        ng_comp=ng_comp,
        initialize=True)
    res = run_full_load(m, solver)
    assert res.solver.status == pyo.SolverStatus.ok
