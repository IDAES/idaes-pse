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
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.fileutils import this_file_dir
from idaes.models_extra.power_generation.properties import (
    FlueGasParameterBlock,
    FlueGasStateBlock,
)
from idaes.core import FlowsheetBlock
import csv
import os
from math import log
from idaes.core.solvers import get_solver

# Set up solver
solver = get_solver()


def read_data(fname, params):
    dfile = os.path.join(this_file_dir(), fname)
    # the data format is data[component][temperature][property]
    data = {
        "N2": {},
        "O2": {},
        "H2O": {},
        "CO2": {},
        "NO": {},
        "SO2": {},
    }

    with open(dfile, "r") as csvfile:
        dat = csv.reader(csvfile, delimiter="\t")
        for i in range(7):
            next(dat)  # skip header
        for row in dat:
            data[row[4]][int(row[0])] = {}
            d = data[row[4]][int(row[0])]
            d["Cp"] = float(row[1])
            d["S"] = float(row[2])
            H = pyo.value(params.cp_mol_ig_comp_coeff_H[(row[4])] * 1000)
            d["H"] = float(row[3]) + H  # H = enthalpy of formation
            d["comp"] = {row[4]: 1.0}

    # Add a mixture to test
    data["mix1"] = {}
    for T in data["N2"]:
        data["mix1"][T] = {}
        d = data["mix1"][T]
        comp = {"N2": 0.2, "O2": 0.2, "H2O": 0.2, "CO2": 0.2, "NO": 0.1, "SO2": 0.1}
        d["Cp"] = sum(data[i][T]["Cp"] * comp[i] for i in comp)
        d["H"] = sum(data[i][T]["H"] * comp[i] for i in comp)
        d["S"] = sum((data[i][T]["S"] + 8.314 * log(comp[i])) * comp[i] for i in comp)
        d["comp"] = comp
    return data


@pytest.mark.integration
def test_thermo():
    # Read in test data and add mixtures
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.params = FlueGasParameterBlock()
    m.fs.state = FlueGasStateBlock(parameters=m.fs.params)
    m.fs.state.pressure.fix(1e5)  # ideal gas properties are pressure independent

    assert_units_consistent(m)

    data = read_data("pure-prop-nist-webbook.csv", m.fs.params)

    assert hasattr(m.fs.state, "cp_mol")
    assert hasattr(m.fs.state, "enth_mol")
    assert hasattr(m.fs.state, "entr_mol")

    m.fs.state.cp_mol = 30
    m.fs.state.enth_mol = 30000
    m.fs.state.entr_mol = 30

    for i in data:
        for T in data[i]:
            m.fs.state.temperature.fix(T)
            test_entropy = True
            for j, f in m.fs.state.flow_mol_comp.items():
                cf = data[i][T]["comp"].get(j, 0)
                if cf <= 0:
                    test_entropy = False
                f.fix(cf)

            if test_entropy:
                m.fs.state.entr_mol.unfix()
                m.fs.state.entropy_correlation.deactivate()
            else:
                m.fs.state.entr_mol.unfix()
                m.fs.state.entropy_correlation.deactivate()
            m.fs.state.initialize()
            solver.solve(m)
            assert data[i][T]["H"] == pytest.approx(
                pyo.value(m.fs.state.enth_mol), rel=1e-2
            )
            assert data[i][T]["Cp"] == pytest.approx(
                pyo.value(m.fs.state.cp_mol), rel=1e-2
            )
            if test_entropy:
                assert data[i][T]["S"] == pytest.approx(
                    pyo.value(m.fs.state.entr_mol), rel=1e-2
                )
