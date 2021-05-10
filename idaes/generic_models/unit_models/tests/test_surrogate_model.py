##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Generic template for a surrogate unit model.
"""

import pytest
from pyomo.environ import ConcreteModel, Constraint, Var, value, SolverFactory
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import SurrogateModel
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import get_default_solver

__author__ = "Jaffer Ghouse"


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()
# -----------------------------------------------------------------------------

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.surrogate = SurrogateModel()
m.fs.surrogate.x_1 = Var(initialize=0.5)
m.fs.surrogate.x_2 = Var(initialize=0.5)
m.fs.surrogate.y_1 = Var(initialize=0.5)
m.fs.surrogate.y_2 = Var(initialize=0.5)
m.fs.surrogate.c1 = Constraint(
    expr=m.fs.surrogate.y_1 == m.fs.surrogate.x_1)
m.fs.surrogate.c2 = Constraint(
    expr=m.fs.surrogate.y_2 == m.fs.surrogate.x_2**2)

inlet_dict = {"x_1": m.fs.surrogate.x_1,
              "x_2": m.fs.surrogate.x_2}
outlet_dict = {"y_1": m.fs.surrogate.y_1,
               "y_2": m.fs.surrogate.y_2}

m.fs.surrogate.add_ports(name="inlet", member_list=inlet_dict)
m.fs.surrogate.add_ports(name="outlet", member_list=outlet_dict)


def my_initialize():
    m.fs.surrogate.c2.deactivate()
    m.fs.surrogate.x_1.fix(1)
    m.fs.surrogate.x_2.fix(2)

    opt = SolverFactory("ipopt")

    opt.solve(m)
    m.fs.surrogate.c2.activate()
    opt.solve(m)


def test_ports():
    assert hasattr(m.fs.surrogate, "inlet")
    assert hasattr(m.fs.surrogate.inlet, "x_1")
    assert hasattr(m.fs.surrogate.inlet, "x_2")

    assert hasattr(m.fs.surrogate, "outlet")
    assert hasattr(m.fs.surrogate.outlet, "y_1")
    assert hasattr(m.fs.surrogate.outlet, "y_2")


def test_build():
    assert hasattr(m.fs.surrogate, "x_1")
    assert hasattr(m.fs.surrogate, "x_2")

    assert hasattr(m.fs.surrogate, "c1")
    assert hasattr(m.fs.surrogate, "c2")

    assert hasattr(m.fs.surrogate, "y_1")
    assert hasattr(m.fs.surrogate, "y_2")


def test_default_initialize():
    m.fs.surrogate.x_1.fix()
    m.fs.surrogate.x_2.fix()

    assert degrees_of_freedom(m) == 0

    m.fs.surrogate.initialize()

    assert value(m.fs.surrogate.x_1) == pytest.approx(0.5, abs=1e-3)
    assert value(m.fs.surrogate.x_2) == pytest.approx(0.5, abs=1e-3)

    assert value(m.fs.surrogate.y_1) == pytest.approx(0.5, abs=1e-3)
    assert value(m.fs.surrogate.y_2) == pytest.approx(0.25, abs=1e-3)


def test_custom_initialize():
    m.fs.surrogate.initialize(custom_initialize=my_initialize())
    assert value(m.fs.surrogate.x_1) == pytest.approx(1, abs=1e-3)
    assert value(m.fs.surrogate.x_2) == pytest.approx(2, abs=1e-3)

    assert value(m.fs.surrogate.y_1) == pytest.approx(1, abs=1e-3)
    assert value(m.fs.surrogate.y_2) == pytest.approx(4, abs=1e-3)
