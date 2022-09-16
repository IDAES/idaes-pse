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
"""Tests that helmholtz specific NTU condenser model."""


import pytest
import pyomo.environ as pyo
import idaes.core
from idaes.models_extra.power_generation.unit_models.helm import HelmNtuCondenser
from idaes.models.properties import iapws95
from idaes.core.solvers import get_solver

solver = get_solver()


@pytest.mark.unit
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_condenser_steady():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmNtuCondenser(
        dynamic=False,
        shell={"has_pressure_change": False, "property_package": m.fs.properties},
        tube={"has_pressure_change": False, "property_package": m.fs.properties},
    )

    m.fs.unit.shell_inlet.flow_mol.fix(100)
    m.fs.unit.shell_inlet.pressure[:] = 9000
    m.fs.unit.shell_inlet.enth_mol.fix(30000)

    m.fs.unit.tube_inlet.flow_mol.fix(1200)
    m.fs.unit.tube_inlet.pressure.fix(101325)
    m.fs.unit.tube_inlet.enth_mol.fix(1800)

    m.fs.unit.area.fix(1000)
    m.fs.unit.overall_heat_transfer_coefficient.fix(500)

    m.fs.unit.initialize(unfix="pressure")
    solver.solve(m)

    assert pyo.value(m.fs.unit.shell_inlet.pressure[0]) == pytest.approx(
        14036.3360, rel=1e-2
    )


@pytest.mark.unit
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_condenser_dynamic():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(
        dynamic=True, time_set=[0, 3], time_units=pyo.units.s
    )
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmNtuCondenser(
        dynamic=False,
        shell={"has_pressure_change": False, "property_package": m.fs.properties},
        tube={"has_pressure_change": False, "property_package": m.fs.properties},
    )

    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m.fs, nfe=4, wrt=m.fs.time, scheme="BACKWARD"
    )

    m.fs.unit.shell_inlet.flow_mol.fix(100)
    m.fs.unit.shell_inlet.pressure[:] = 9000
    m.fs.unit.shell_inlet.enth_mol.fix(30000)

    m.fs.unit.tube_inlet.flow_mol.fix(1200)
    m.fs.unit.tube_inlet.pressure.fix(101325)
    m.fs.unit.tube_inlet.enth_mol.fix(1800)

    m.fs.unit.area.fix(1000)
    m.fs.unit.overall_heat_transfer_coefficient.fix(500)

    m.fs.unit.initialize(unfix="pressure")
    solver.solve(m)

    for t in m.fs.unit.shell_inlet.pressure:
        assert pyo.value(m.fs.unit.shell_inlet.pressure[t]) == pytest.approx(
            14036.3360, rel=1e-2
        )
