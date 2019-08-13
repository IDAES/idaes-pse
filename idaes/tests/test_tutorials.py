##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests that tutorials run.

Author: Andrew Lee
"""
# stdlib
# third-party
import pytest
from pyomo.environ import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
# package
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.examples.tutorials import (
    Tutorial_1_Basic_Flowsheets,
    Tutorial_2_Basic_Flowsheet_Optimization,
    Tutorial_3_Dynamic_Flowsheets,
)

# See if ipopt is available and set up solver
if not SolverFactory("ipopt").available():
    solver = None
else:
    solver = True


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_tutorial_1():
    m, results = Tutorial_1_Basic_Flowsheets.main()

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert degrees_of_freedom(m) == 0

    assert m.fs.Tank2.outlet.flow_vol[0].value == pytest.approx(1.0, abs=1e-2)
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "Ethanol"].value == pytest.approx(
        89.628, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "EthylAcetate"].value == pytest.approx(
        10.372, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "NaOH"].value == pytest.approx(
        10.372, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "SodiumAcetate"].value == pytest.approx(
        89.628, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "H2O"].value == pytest.approx(
        55388.0, abs=1
    )
    assert m.fs.Tank2.outlet.pressure[0].value == pytest.approx(101325, abs=1)
    assert m.fs.Tank2.outlet.temperature[0].value == pytest.approx(304.20, abs=1e-1)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_tutorial_2():
    m, results = Tutorial_2_Basic_Flowsheet_Optimization.main()

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert degrees_of_freedom(m) == 1

    assert m.fs.Tank2.outlet.flow_vol[0].value == pytest.approx(1.0, abs=1e-2)
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "Ethanol"].value == pytest.approx(
        92.106, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "EthylAcetate"].value == pytest.approx(
        7.894, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "NaOH"].value == pytest.approx(
        7.894, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "SodiumAcetate"].value == pytest.approx(
        92.106, abs=1e-2
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[0, "H2O"].value == pytest.approx(
        55388.0, abs=1
    )
    assert m.fs.Tank2.outlet.pressure[0].value == pytest.approx(101325, abs=1)
    assert m.fs.Tank2.outlet.temperature[0].value == pytest.approx(304.23, abs=1e-1)

    assert m.fs.Tank1.volume[0].value == pytest.approx(1.215, abs=1e-2)
    assert m.fs.Tank2.volume[0].value == pytest.approx(1.785, abs=1e-2)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_tutorial_3():
    m, results = Tutorial_3_Dynamic_Flowsheets.main()

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    assert degrees_of_freedom(m) == 0

    assert m.fs.time.last() == 20.0
    assert len(m.fs.time) == 51

    assert m.fs.Tank2.outlet.flow_vol[20.0].value == pytest.approx(1.0, abs=1e-2)
    assert m.fs.Tank2.outlet.conc_mol_comp[20.0, "Ethanol"].value == pytest.approx(
        43.62, abs=1e-1
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[20.0, "EthylAcetate"].value == pytest.approx(
        1.65, abs=1e-1
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[20.0, "NaOH"].value == pytest.approx(
        6.38, abs=1e-1
    )
    assert m.fs.Tank2.outlet.conc_mol_comp[
        20.0, "SodiumAcetate"
    ].value == pytest.approx(43.62, abs=1e-1)
    assert m.fs.Tank2.outlet.conc_mol_comp[20.0, "H2O"].value == pytest.approx(
        55388.0, abs=1
    )
    assert m.fs.Tank2.outlet.pressure[20.0].value == pytest.approx(101325, abs=1)
    assert m.fs.Tank2.outlet.temperature[20.0].value == pytest.approx(303.66, abs=1e-1)
