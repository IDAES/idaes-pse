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
"""Tests that helmholtz specific NTU condenser model."""


import pytest
import pyomo.environ as pyo
import idaes.core
from idaes.power_generation.unit_models.helm import HelmNtuCondenser
from idaes.generic_models.properties import iapws95
from idaes.core.util import get_solver

solver = get_solver()

@pytest.mark.unit
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_condenser_steady():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmNtuCondenser(
        default={
            "dynamic": False,
            "shell": {
                "has_pressure_change": False,
                "property_package": m.fs.properties,
            },
            "tube": {
                "has_pressure_change": False,
                "property_package": m.fs.properties,
            },
        }
    )

    m.fs.unit.shell_inlet.flow_mol.fix(100)
    m.fs.unit.shell_inlet.pressure[:] = 9000
    m.fs.unit.shell_inlet.enth_mol.fix(30000)

    m.fs.unit.tube_inlet.flow_mol.fix(1200)
    m.fs.unit.tube_inlet.pressure.fix(101325)
    m.fs.unit.tube_inlet.enth_mol.fix(1800)

    m.fs.unit.area.fix(1000)
    m.fs.unit.overall_heat_transfer_coefficient.fix(500)

    m.fs.unit.initialize(unfix='pressure')
    solver.solve(m)

    assert (pyo.value(m.fs.unit.shell_inlet.pressure[0]) ==
        pytest.approx(14036.3360, rel=1e-2))


@pytest.mark.unit
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_condenser_dynamic():
    m = pyo.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock(
        default={"dynamic": True, "time_set":[0,3]})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmNtuCondenser(
        default={
            "dynamic": False,
            "shell": {
                "has_pressure_change": False,
                "property_package": m.fs.properties,
            },
            "tube": {
                "has_pressure_change": False,
                "property_package": m.fs.properties,
            },
        }
    )

    pyo.TransformationFactory('dae.finite_difference').apply_to(
        m.fs, nfe=4, wrt=m.fs.time, scheme='BACKWARD')

    m.fs.unit.shell_inlet.flow_mol.fix(100)
    m.fs.unit.shell_inlet.pressure[:] = 9000
    m.fs.unit.shell_inlet.enth_mol.fix(30000)

    m.fs.unit.tube_inlet.flow_mol.fix(1200)
    m.fs.unit.tube_inlet.pressure.fix(101325)
    m.fs.unit.tube_inlet.enth_mol.fix(1800)

    m.fs.unit.area.fix(1000)
    m.fs.unit.overall_heat_transfer_coefficient.fix(500)

    m.fs.unit.initialize(unfix='pressure')
    solver.solve(m)

    for t in m.fs.unit.shell_inlet.pressure:
        assert (pyo.value(m.fs.unit.shell_inlet.pressure[t]) ==
            pytest.approx(14036.3360, rel=1e-2))
