#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pytest
from pyomo.environ import ConcreteModel, Constraint, Expression, Var
from idaes.apps.grid_integration.pricetaker.design_and_operation_models import (
    DesignModel,
    OperationModel,
)


@pytest.mark.unit
def test_design_model_class(caplog):
    """Tests the DesignModel class"""

    def _my_design_model(m, p_min, p_max, cost):
        m.power = Var()
        m.min_capacity = Constraint(expr=p_min * m.install_unit <= m.power)
        m.max_capacity = Constraint(expr=m.power <= p_max * m.install_unit)

        # capex and fom must either be a constant, or Var, or Expression
        m.capex = Expression(expr=cost["capex"] * m.power)
        m.fom = Expression(expr=cost["fom"] * m.power)

    blk = ConcreteModel()
    blk.unit_1 = DesignModel(
        model_func=_my_design_model,
        model_args={
            "p_min": 150,
            "p_max": 600,
            "cost": {"capex": 10, "fom": 1},
        },
    )

    # Begin checks
    for attr in ["install_unit", "power", "min_capacity", "max_capacity"]:
        assert hasattr(blk.unit_1, attr)

    assert isinstance(blk.unit_1.capex, Expression)
    assert isinstance(blk.unit_1.fom, Expression)
    assert "Setting the capital cost of the unit to zero." not in caplog.text
    assert "Setting the fixed O&M cost of the unit to zero." not in caplog.text

    # Test empty design model
    blk.unit_2 = DesignModel()
    assert hasattr(blk.unit_2, "install_unit")
    for attr in ["power", "min_capacity", "max_capacity", "capex", "fom"]:
        assert not hasattr(blk.unit_2, attr)

    assert "Setting the capital cost of the unit to zero." not in caplog.text
    assert "Setting the fixed O&M cost of the unit to zero." not in caplog.text

    # Test warning messages
    def dummy_func(_):
        pass

    blk.unit_3 = DesignModel(model_func=dummy_func)
    assert hasattr(blk.unit_3, "install_unit")
    assert (
        "'capex' attribute is not set for the design model "
        "unit_3. Setting the capital cost of the unit to zero."
    ) in caplog.text
    assert (
        "'fom' attribute is not set for the design model "
        "unit_3. Setting the fixed O&M cost of the unit to zero."
    ) in caplog.text
    assert blk.unit_3.capex == 0
    assert blk.unit_3.fom == 0


@pytest.mark.unit
def test_design_model_class_logger_message1(caplog):
    caplog.clear()

    blk = ConcreteModel()
    blk.unit_1 = DesignModel(
        model_args={
            "p_min": 150,
            "p_max": 600,
            "cost": {"capex": 10, "fom": 1},
        },
    )

    assert (
        "The function that builds the design model is not specified."
        "model_func must declare all the necessary design variables,"
        "relations among design variables, capital cost correlations,"
        "and fixed operating and maintenance cost correlations." in caplog.text
    )


@pytest.mark.unit
def test_operation_model_class():
    """Tests the OperationModel class"""

    def des_model(m):
        m.power = Var()
        m.capex = 0
        m.fom = 0

    def op_model(m, des_blk):
        m.power = Var()
        m.scaled_power = Expression(
            expr=m.power / des_blk.power,
            doc="Instantaneous power scaled with the design power",
        )

    blk = ConcreteModel()
    blk.unit_1_design = DesignModel(model_func=des_model)
    blk.unit_1_op = OperationModel(
        model_func=op_model, model_args={"des_blk": blk.unit_1_design}
    )

    # Begin checks
    for attr in ["op_mode", "startup", "shutdown", "power", "LMP"]:
        assert hasattr(blk.unit_1_op, attr)

    assert isinstance(blk.unit_1_op.scaled_power, Expression)

    # Empty operation model
    blk.unit_2_op = OperationModel(declare_op_vars=False, declare_lmp_param=False)
    for attr in ["op_mode", "startup", "shutdown", "power", "LMP"]:
        assert not hasattr(blk.unit_2_op, attr)


@pytest.mark.unit
def test_operation_model_class_logger_message1(caplog):
    caplog.clear()

    def des_model(m):
        m.power = Var()
        m.capex = 0
        m.fom = 0

    blk = ConcreteModel()
    blk.unit_1_design = DesignModel(model_func=des_model)
    blk.unit_1_op = OperationModel(model_args={"des_blk": blk.unit_1_design})

    assert (
        "The function that builds the operation model is not specified."
        "model_func must declare all the necessary operation variables,"
        "relations among operation variables, and variable"
        "operating and maintenance cost correlations." in caplog.text
    )
