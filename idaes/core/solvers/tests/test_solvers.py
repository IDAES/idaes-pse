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
import pyomo.environ as pyo
import pytest
from idaes.core.solvers.features import lp, milp, nlp, minlp, nle, dae
from idaes.core.solvers import ipopt_has_linear_solver
from idaes.core.solvers import petsc


@pytest.mark.unit
def test_petsc_available():
    if not pyo.SolverFactory("petsc_snes").available():
        raise RuntimeError("Could not find petsc (petsc is an optional extra).")


@pytest.mark.unit
def test_couenne_available():
    if not pyo.SolverFactory("couenne").available():
        raise RuntimeError("Could not find couenne.")


@pytest.mark.unit
def test_bonmin_available():
    if not pyo.SolverFactory("bonmin").available():
        raise RuntimeError("Could not find bonmin.")


@pytest.mark.unit
def test_sipopt_available():
    if not pyo.SolverFactory("ipopt_sens").available():
        raise RuntimeError("Could not find ipopt_sens.")


@pytest.mark.unit
def test_ipopt_idaes_available():
    """
    Tries to set-up the IPOPT with the IDAES SolverFactory wrapper
    """
    if not pyo.SolverFactory("ipopt").available():
        raise RuntimeError(
            "Could not find IPOPT. Users are strongly encouraged to have a "
            "version of IPOPT available, as it is the default solver assumed "
            "by many IDAES examples and tests. See the IDAES install "
            "documentation for instructions on how to get IPOPT."
        )


@pytest.mark.unit
def test_cbc_available():
    if not pyo.SolverFactory("cbc").available():
        raise RuntimeError("Could not find cbc.")


@pytest.mark.unit
def test_clp_available():
    if not pyo.SolverFactory("clp").available():
        raise RuntimeError("Could not find clp.")


@pytest.mark.unit
def test_sipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = nlp()
    solver = pyo.SolverFactory("ipopt_sens")
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)


@pytest.mark.unit
def test_ipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = nlp()
    solver = pyo.SolverFactory("ipopt")
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)


@pytest.mark.unit
def test_ipopt_has_ma27():
    if not ipopt_has_linear_solver("ma27"):
        raise RuntimeError(
            "The ma27 linear solver is not available to Ipopt. Models may solve"
            " more reliably with HSL linear solvers see https://www.hsl.rl.ac.uk/,"
            " or use solvers distributed by the IDAES project. See IDAES install"
            " guide."
        )


@pytest.mark.unit
def test_ipopt_has_ma57():
    if not ipopt_has_linear_solver("ma57"):
        raise RuntimeError("The ma57 linear solver is not available to Ipopt.")


@pytest.mark.skip
@pytest.mark.unit
def test_ipopt_has_ma77():
    if not ipopt_has_linear_solver("ma77"):
        raise RuntimeError("The ma77 linear solver is not available to Ipopt.")


@pytest.mark.skip
@pytest.mark.unit
def test_ipopt_has_ma86():
    if not ipopt_has_linear_solver("ma86"):
        raise RuntimeError("The ma86 linear solver is not available to Ipopt.")


@pytest.mark.unit
def test_ipopt_has_ma97():
    if not ipopt_has_linear_solver("ma97"):
        raise RuntimeError("The ma97 linear solver is not available to Ipopt.")


@pytest.mark.unit
def test_ipopt_has_mumps():
    if not ipopt_has_linear_solver("mumps"):
        raise RuntimeError("The mumps linear solver is not available to Ipopt.")


@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_petsc_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = nle()
    solver = pyo.SolverFactory("petsc_snes")
    solver.solve(m, tee=True)
    assert pytest.approx(x) == pyo.value(m.x)


@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_petsc_dae_idaes_solve():
    """
    Check that the PETSc DAE solver works.
    """
    m, y1, y2, y3, y4, y5, y6 = dae()
    petsc.petsc_dae_by_time_element(
        m,
        time=m.t,
        ts_options={
            "--ts_type": "cn",  # Crank–Nicolson
            "--ts_adapt_type": "basic",
            "--ts_dt": 0.1,
        },
    )
    assert pytest.approx(y1, rel=1e-3) == pyo.value(m.y[m.t.last(), 1])
    assert pytest.approx(y2, rel=1e-3) == pyo.value(m.y[m.t.last(), 2])
    assert pytest.approx(y3, rel=1e-3) == pyo.value(m.y[m.t.last(), 3])
    assert pytest.approx(y4, rel=1e-3) == pyo.value(m.y[m.t.last(), 4])
    assert pytest.approx(y5, rel=1e-3) == pyo.value(m.y[m.t.last(), 5])
    assert pytest.approx(y6, rel=1e-3) == pyo.value(m.y6[m.t.last()])


@pytest.mark.unit
def test_bonmin_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x, i = minlp()
    solver = pyo.SolverFactory("bonmin")
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)
    assert i == pyo.value(m.i)


@pytest.mark.unit
def test_couenne_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x, i = minlp()
    solver = pyo.SolverFactory("couenne")
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)
    assert i == pyo.value(m.i)


@pytest.mark.unit
def test_cbc_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = milp()
    solver = pyo.SolverFactory("cbc")
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)


@pytest.mark.unit
def test_clp_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = lp()
    solver = pyo.SolverFactory("clp")
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)
