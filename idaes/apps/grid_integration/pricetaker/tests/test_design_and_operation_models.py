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
from pyomo.environ import ConcreteModel, Constraint, Expression, Param, Var

# pylint: disable = no-name-in-module,
from idaes.apps.grid_integration.pricetaker.design_and_operation_models import (
    DesignModel,
    OperationModel,
    StorageModel,
    _format_data,
    _is_valid_data_type_for_storage_model,
    is_valid_variable_design_data,
    is_valid_polynomial_surrogate_data,
)
from idaes.core.util.config import ConfigurationError


@pytest.mark.unit
def test_format_data():
    """Tests the _format_data function"""
    m = ConcreteModel()
    m.x = Var()
    m.p1 = Param(initialize=42)

    # Test scalar argument input
    assert _format_data(5) == [0, 5]
    assert _format_data(4.6) == [0, 4.6]
    assert _format_data(m.x) == [0, m.x]
    assert _format_data(m.p1) == [0, m.p1]

    # Test array-type argument input
    assert _format_data([0, 2, 4]) == [0, 2, 4]
    assert _format_data((0, 2.2, 4.6, 8.9)) == (0, 2.2, 4.6, 8.9)

    # Test error for unsupported data structures
    with pytest.raises(
        ConfigurationError,
        match="Unrecognized data structure <class 'str'> for auxiliary variable coefficients.",
    ):
        _format_data("hello!")


@pytest.mark.unit
def test_is_valid_variable_design_data():
    """Tests the is_valid_variable_design_data function"""

    # Test design_var missing error
    with pytest.raises(ConfigurationError, match="design_var is not specified"):
        is_valid_variable_design_data({"a": 5, "b": [2.2, 4.4]})

    # Test design_var_bounds missing error
    with pytest.raises(ConfigurationError, match="design_var_bounds is not specified"):
        is_valid_variable_design_data({"design_var": "power", "a": 5, "b": [2.2, 4.4]})

    # Test if the correct data structure is returned
    data = is_valid_variable_design_data(
        {
            "design_var": "power",
            "design_var_bounds": [40, 100],
            "a": 5.4,
            "b": [2.2, 4.4],
        }
    )
    assert data == {
        "design_var": "power",
        "design_var_bounds": [40, 100],
        "auxiliary_vars": {
            "a": [0, 5.4],
            "b": [2.2, 4.4],
        },
    }


@pytest.mark.unit
def test_is_valid_polynomial_surrogate_data():
    """Tests the is_valid_surrogate_data function"""
    # Test the operation_var missing error
    with pytest.raises(ConfigurationError, match="operation_var is not specified"):
        is_valid_polynomial_surrogate_data({"a": 5.4, "b": [2.2, 4.4]})

    # Test if the correct data structure is returned
    data = is_valid_polynomial_surrogate_data(
        {
            "operation_var": "power",
            "a": 5.4,
            "b": [2.2, 4.4],
        }
    )
    assert data == {
        "operation_var": "power",
        "auxiliary_vars": {
            "a": [0, 5.4],
            "b": [2.2, 4.4],
        },
    }


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
    """
    Tests the DesignModel class logger message when no
    input is specified
    """
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
def test_design_model_with_fixed_design():
    """Tests the DesignModel class with fixed_design_data argument"""
    blk = ConcreteModel()
    blk.d1 = DesignModel(fixed_design_data={"capacity": 650, "capex": 40, "fom": 35})
    assert isinstance(blk.d1.capacity, Param)
    assert isinstance(blk.d1.capex, Param)
    assert isinstance(blk.d1.fom, Param)
    assert blk.d1.capacity.value == 650
    assert blk.d1.capex.value == 40
    assert blk.d1.fom.value == 35


@pytest.mark.unit
def test_design_model_with_variable_design():
    """Tests the DesignModel class with variable_design_data argument"""
    blk = ConcreteModel()
    blk.d1 = DesignModel(
        variable_design_data={
            "design_var": "power",
            "design_var_bounds": [150, 650],
            "capex": [4, 5, 6, 7],
            "fom": [2, 3],
            "ng_flow": 42,
        }
    )

    assert isinstance(blk.d1.power, Var)
    assert blk.d1.power.lb == 0
    assert blk.d1.power.ub == 650

    assert isinstance(blk.d1.design_lb_constraint, Constraint)
    assert isinstance(blk.d1.design_ub_constraint, Constraint)
    assert isinstance(blk.d1.capex, Expression)
    assert isinstance(blk.d1.fom, Expression)
    assert isinstance(blk.d1.ng_flow, Expression)

    assert str(blk.d1.design_lb_constraint.expr) == "150*d1.install_unit  <=  d1.power"
    assert str(blk.d1.design_ub_constraint.expr) == "d1.power  <=  650*d1.install_unit"
    assert str(blk.d1.capex.expr) == (
        "5*d1.power + 6*d1.power**2 + 7*d1.power**3 + 4*d1.install_unit"
    )
    assert str(blk.d1.fom.expr) == "2*d1.install_unit + 3*d1.power"
    assert str(blk.d1.ng_flow.expr) == "0*d1.install_unit + 42*d1.power"


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

    # test the multiple startup types
    blk.unit_3_op = OperationModel(
        model_func=op_model,
        model_args={"des_blk": blk.unit_1_design},
        startup_types={'hot': 4, 'warm': 8, 'cold': 12}
    )
    assert hasattr(blk.unit_3_op, "startup_type_vars")


@pytest.mark.unit
def test_operation_model_class_logger_message1(caplog):
    """
    Testst the OperationModel class when no input
    is specified while instantiating.
    """
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


@pytest.mark.unit
def test_operation_model_with_surrogates():
    """Tests the OperationModel with additional surrogate models"""
    m = ConcreteModel()
    m.d1 = DesignModel(
        variable_design_data={
            "design_var": "capacity",
            "design_var_bounds": [150, 650],
            "capex": 100,
            "fom": 50,
        }
    )
    m.op_blk = OperationModel(
        polynomial_surrogate_data={
            "operation_var": "power",
            "vom": [4, 5, 6, 7],
        }
    )
    assert isinstance(m.op_blk.power, Var)
    assert isinstance(m.op_blk.vom, Var)
    assert isinstance(m.op_blk.compute_vom, Constraint)
    assert str(m.op_blk.compute_vom.expr) == (
        "op_blk.vom  ==  5*op_blk.power + 6*op_blk.power**2"
        " + 7*op_blk.power**3 + 4*op_blk.op_mode"
    )

    # Test the build_polynomial_surrogates method
    m.op_blk.build_polynomial_surrogates(
        surrogates={"ng_flow": [2 * m.d1.capacity, 42]},
        op_var=m.op_blk.power,
    )
    assert isinstance(m.op_blk.ng_flow, Var)
    assert isinstance(m.op_blk.compute_ng_flow, Constraint)
    assert str(m.op_blk.compute_ng_flow.expr) == (
        "op_blk.ng_flow  ==  2*d1.capacity*op_blk.op_mode + 42*op_blk.power"
    )

    # Test the build_expressions method with auxiliary variables
    m.op_blk.build_expressions(
        expressions={"e1": 5 * m.op_blk.power + 56, "e2": m.op_blk.power**5.26},
        declare_variables=True,
    )
    assert isinstance(m.op_blk.e1, Var)
    assert isinstance(m.op_blk.e2, Var)
    assert isinstance(m.op_blk.compute_e1, Constraint)
    assert isinstance(m.op_blk.compute_e2, Constraint)
    assert str(m.op_blk.compute_e1.expr) == "op_blk.e1  ==  5*op_blk.power + 56"
    assert str(m.op_blk.compute_e2.expr) == "op_blk.e2  ==  op_blk.power**5.26"

    # Test the build_expressions method without auxiliary variables
    m.op_blk.build_expressions(
        expressions={"e3": 5 * m.op_blk.power + 56, "e4": m.op_blk.power**5.26}
    )
    assert isinstance(m.op_blk.e3, Expression)
    assert isinstance(m.op_blk.e4, Expression)
    assert str(m.op_blk.e3.expr) == "5*op_blk.power + 56"
    assert str(m.op_blk.e4.expr) == "op_blk.power**5.26"


@pytest.mark.unit
def test_is_valid_data_type_for_storage_model():
    """Tests the _is_valid_data_type_for_storage_model function"""
    m = ConcreteModel()
    m.blk = DesignModel()
    m.blk.x = Var([1, 2, 3])
    m.y = Var()
    m.blk.p1 = Param(initialize=42)
    m.p2 = Param([2, 3, 4], initialize={2: 1, 3: 2, 4: 3})
    m.blk.e1 = Expression(expr=m.blk.x[1] + m.blk.x[2])
    m.e2 = Expression([1, 2, 3], rule=lambda b, i: b.blk.x[i] ** i)

    # Test the invalid pyomo component error
    with pytest.raises(
        ConfigurationError,
        match="Received an invalid Pyomo component as an argument for StorageModel",
    ):
        _is_valid_data_type_for_storage_model(m.blk)

    # Test the unsupported data type error
    with pytest.raises(
        ConfigurationError, match="Received an unsupported data type for StorageModel"
    ):
        _is_valid_data_type_for_storage_model("hello!")

    # Test the function
    assert _is_valid_data_type_for_storage_model(5) == 5
    assert _is_valid_data_type_for_storage_model(2.5) == 2.5
    assert _is_valid_data_type_for_storage_model(m.blk.p1) is m.blk.p1
    assert _is_valid_data_type_for_storage_model(m.p2[3]) is m.p2[3]
    assert _is_valid_data_type_for_storage_model(m.blk.x[1]) is m.blk.x[1]
    assert _is_valid_data_type_for_storage_model(m.y) is m.y
    assert _is_valid_data_type_for_storage_model(m.blk.e1) is m.blk.e1
    assert _is_valid_data_type_for_storage_model(m.e2[2]) is m.e2[2]

    # Test unnamed expressions
    assert (
        str(_is_valid_data_type_for_storage_model(m.blk.x[1] + m.blk.x[2] + m.blk.x[3]))
        == "blk.x[1] + blk.x[2] + blk.x[3]"
    )


@pytest.mark.unit
def test_storage_model_fixed_design():
    """Tests the StorageModel class with int/float-type arguments"""
    m = ConcreteModel()
    m.sb = StorageModel(
        time_interval=0.14,
        discharge_efficiency=0.8,
        min_holdup=50.5,
        max_holdup=100,
        max_charge_rate=10.5,
        max_discharge_rate=20,
    )

    assert m.sb.config.charge_efficiency == 1
    assert m.sb.config.discharge_efficiency == 0.8
    assert not hasattr(m.sb, "charge_rate_ub_con")
    assert not hasattr(m.sb, "discharge_rate_ub_con")
    assert not hasattr(m.sb, "final_holdup_ub_con")
    assert not hasattr(m.sb, "initial_holdup_ub_con")
    assert not hasattr(m.sb, "final_holdup_lb_con")
    assert not hasattr(m.sb, "initial_holdup_lb_con")

    assert m.sb.final_holdup.ub == 100
    assert m.sb.initial_holdup.ub == 100
    assert m.sb.final_holdup.lb == 50.5
    assert m.sb.initial_holdup.lb == 50.5
    assert m.sb.charge_rate.ub == 10.5
    assert m.sb.discharge_rate.ub == 20
    assert m.sb.charge_rate.lb == 0
    assert m.sb.discharge_rate.lb == 0

    assert str(m.sb.track_holdup_constraint.expr) == (
        "sb.final_holdup - sb.initial_holdup  ==  "
        "(sb.charge_rate - 1.25*sb.discharge_rate)*0.14"
    )


@pytest.mark.unit
def test_storage_model_variable_design():
    """Tests the storage model with a variable design"""
    m = ConcreteModel()
    m.capacity = Var()
    m.charge_limit = Var()

    m.sb = StorageModel(
        time_interval=0.5,
        charge_efficiency=0.9,
        discharge_efficiency=0.8,
        min_holdup=0.2 * m.capacity,
        max_holdup=m.capacity,
        max_charge_rate=m.charge_limit,
        max_discharge_rate=0.8 * m.charge_limit,
    )

    assert m.sb.config.charge_efficiency == 0.9
    assert m.sb.config.discharge_efficiency == 0.8
    assert str(m.sb.charge_rate_ub_con.expr) == "sb.charge_rate  <=  charge_limit"
    assert str(m.sb.discharge_rate_ub_con.expr) == (
        "sb.discharge_rate  <=  0.8*charge_limit"
    )
    assert str(m.sb.final_holdup_ub_con.expr) == "sb.final_holdup  <=  capacity"
    assert str(m.sb.initial_holdup_ub_con.expr) == "sb.initial_holdup  <=  capacity"
    assert str(m.sb.final_holdup_lb_con.expr) == "0.2*capacity  <=  sb.final_holdup"
    assert str(m.sb.initial_holdup_lb_con.expr) == "0.2*capacity  <=  sb.initial_holdup"

    assert m.sb.final_holdup.ub is None
    assert m.sb.initial_holdup.ub is None
    assert m.sb.charge_rate.ub is None
    assert m.sb.discharge_rate.ub is None
    assert m.sb.final_holdup.lb == 0
    assert m.sb.initial_holdup.lb == 0
    assert m.sb.charge_rate.lb == 0
    assert m.sb.discharge_rate.lb == 0

    assert str(m.sb.track_holdup_constraint.expr) == (
        "sb.final_holdup - sb.initial_holdup  ==  "
        "(0.9*sb.charge_rate - 1.25*sb.discharge_rate)*0.5"
    )
