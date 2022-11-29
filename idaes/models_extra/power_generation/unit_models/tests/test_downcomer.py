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
Downcomer model
Main assumptions:
* Heat is given (zero if adiabatic)
* Calculate pressure change due to friction and gravity
* Assume enthalpy_in == enthalpy_out + heat

Created on Aug 27, 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.models_extra.power_generation.unit_models.downcomer import Downcomer

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_downcomer():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = Downcomer(
        property_package=m.fs.properties, has_holdup=False, has_heat_transfer=True
    )
    return m


@pytest.mark.unit
def test_basic_build(build_downcomer):
    """Make a model and make sure it doesn't throw exception"""
    m = build_downcomer
    assert degrees_of_freedom(m) == 7
    # Check unit config arguments
    assert len(m.fs.unit.config) == 8
    assert m.fs.unit.config.has_heat_transfer
    assert isinstance(m.fs.unit.heat_duty, pyo.Var)
    assert isinstance(m.fs.unit.deltaP, pyo.Var)
    assert m.fs.unit.config.property_package is m.fs.properties


@pytest.mark.integration
def test_units(build_downcomer):
    assert_units_consistent(build_downcomer)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_downcomer):
    m = build_downcomer
    m.fs.unit.diameter.fix(0.3)
    m.fs.unit.number_downcomers.fix(6)  # number of downcomers
    m.fs.unit.height.fix(40.1)
    m.fs.unit.heat_duty[:].fix(0.0)

    # NETL Baseline values
    m.fs.unit.inlet.enth_mol.fix(24944.7)
    m.fs.unit.inlet.flow_mol.fix(49552.8)
    m.fs.unit.inlet.pressure.fix(12025072.9)

    state_args = {"flow_mol": 49552.8, "pressure": 12025072.9, "enth_mol": 24944.7}
    initialization_tester(build_downcomer, dof=0, state_args=state_args)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve_unit(build_downcomer):
    m = build_downcomer

    m.fs.unit.inlet.enth_mol.fix()
    m.fs.unit.inlet.flow_mol.fix()
    m.fs.unit.inlet.pressure.fix()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    # check material balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].flow_mol
                - m.fs.unit.control_volume.properties_out[0].flow_mol
            ),
            abs=1e-3,
        )
        == 0
    )

    # pressure drop
    assert pytest.approx(273583.120306, rel=1e-5) == pyo.value(m.fs.unit.deltaP[0])
    # check energy balance
    assert pytest.approx(
        pyo.value(m.fs.unit.control_volume.properties_in[0].enth_mol), abs=1e-3
    ) == pyo.value(m.fs.unit.control_volume.properties_out[0].enth_mol)

    assert pytest.approx(580.83299, abs=1e-3) == pyo.value(
        m.fs.unit.control_volume.properties_out[0].temperature
    )
