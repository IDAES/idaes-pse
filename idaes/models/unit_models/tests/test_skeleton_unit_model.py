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
Test for skeleton unit model.
"""

import pytest
from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Var,
    Set,
    value,
)
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import SkeletonUnitModel
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError

__author__ = "Jaffer Ghouse"


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
# -----------------------------------------------------------------------------


class TestSkeletonDefault(object):
    """
    Test for SkeletonUnitModel without callback for initialization

    """

    @pytest.fixture(scope="class")
    def skeleton_default(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.skeleton = SkeletonUnitModel(dynamic=False)
        m.fs.skeleton.comp_list = Set(initialize=["c1", "c2"])

        # input vars for skeleton
        m.fs.skeleton.flow_comp_in = Var(
            m.fs.time, m.fs.skeleton.comp_list, initialize=1.0
        )
        m.fs.skeleton.temperature_in = Var(m.fs.time, initialize=298.15)
        m.fs.skeleton.pressure_in = Var(m.fs.time, initialize=101)

        # output vars for skeleton
        m.fs.skeleton.flow_comp_out = Var(
            m.fs.time, m.fs.skeleton.comp_list, initialize=1.0
        )
        m.fs.skeleton.temperature_out = Var(m.fs.time, initialize=298.15)
        m.fs.skeleton.pressure_out = Var(m.fs.time, initialize=101)

        # Surrogate model equations
        # Flow equation

        def rule_flow(m, t, i):
            return m.flow_comp_out[t, i] == m.flow_comp_in[t, i]

        m.fs.skeleton.eq_flow = Constraint(
            m.fs.time, m.fs.skeleton.comp_list, rule=rule_flow
        )

        m.fs.skeleton.eq_temperature = Constraint(
            expr=m.fs.skeleton.temperature_out[0]
            == m.fs.skeleton.temperature_in[0] + 10
        )
        m.fs.skeleton.eq_pressure = Constraint(
            expr=m.fs.skeleton.pressure_out[0] == m.fs.skeleton.pressure_in[0] - 10
        )

        inlet_dict = {
            "flow_mol_comp": m.fs.skeleton.flow_comp_in,
            "temperature": m.fs.skeleton.temperature_in,
            "pressure": m.fs.skeleton.pressure_in,
        }
        outlet_dict = {
            "flow_mol_comp": m.fs.skeleton.flow_comp_out,
            "temperature": m.fs.skeleton.temperature_out,
            "pressure": m.fs.skeleton.pressure_out,
        }

        m.fs.skeleton.add_ports(name="inlet", member_dict=inlet_dict)
        m.fs.skeleton.add_ports(name="outlet", member_dict=outlet_dict)

        return m

    @pytest.mark.unit
    def test_ports(self, skeleton_default):

        # Check inlet port and port members
        assert hasattr(skeleton_default.fs.skeleton, "inlet")
        assert len(skeleton_default.fs.skeleton.inlet.vars) == 3
        assert hasattr(skeleton_default.fs.skeleton.inlet, "flow_mol_comp")
        assert hasattr(skeleton_default.fs.skeleton.inlet, "temperature")
        assert hasattr(skeleton_default.fs.skeleton.inlet, "pressure")

        # Check outlet port and port members
        assert hasattr(skeleton_default.fs.skeleton, "outlet")
        assert len(skeleton_default.fs.skeleton.outlet.vars) == 3
        assert hasattr(skeleton_default.fs.skeleton.outlet, "flow_mol_comp")
        assert hasattr(skeleton_default.fs.skeleton.outlet, "temperature")
        assert hasattr(skeleton_default.fs.skeleton.outlet, "pressure")

    @pytest.mark.unit
    def test_build(self, skeleton_default):
        # Check build of model
        assert hasattr(skeleton_default.fs.skeleton, "flow_comp_in")
        assert hasattr(skeleton_default.fs.skeleton, "temperature_in")
        assert hasattr(skeleton_default.fs.skeleton, "pressure_in")

        assert hasattr(skeleton_default.fs.skeleton, "eq_flow")

        assert hasattr(skeleton_default.fs.skeleton, "flow_comp_out")
        assert hasattr(skeleton_default.fs.skeleton, "temperature_out")
        assert hasattr(skeleton_default.fs.skeleton, "pressure_out")

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_exception(self, skeleton_default):
        # Check default initialize method raises ConfigurationError

        with pytest.raises(
            ConfigurationError,
            match="Degrees of freedom is not zero during start of "
            "initialization. Fix/unfix appropriate number of variables "
            "to result in zero degrees of freedom for initialization.",
        ):
            skeleton_default.fs.skeleton.initialize()

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_default_initialize(self, skeleton_default):

        skeleton_default.fs.skeleton.inlet.flow_mol_comp[0, "c1"].fix(2)
        skeleton_default.fs.skeleton.inlet.flow_mol_comp[0, "c2"].fix(2)
        skeleton_default.fs.skeleton.inlet.temperature.fix(325)
        skeleton_default.fs.skeleton.inlet.pressure.fix(200)

        assert degrees_of_freedom(skeleton_default) == 0

        skeleton_default.fs.skeleton.initialize()

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, skeleton_default):
        results = solver.solve(skeleton_default)

        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, skeleton_default):

        assert value(
            skeleton_default.fs.skeleton.outlet.flow_mol_comp[0, "c1"]
        ) == pytest.approx(2, abs=1e-3)
        assert value(
            skeleton_default.fs.skeleton.outlet.flow_mol_comp[0, "c2"]
        ) == pytest.approx(2, abs=1e-3)

        assert value(
            skeleton_default.fs.skeleton.outlet.temperature[0]
        ) == pytest.approx(335, abs=1e-3)
        assert value(skeleton_default.fs.skeleton.outlet.pressure[0]) == pytest.approx(
            190, abs=1e-3
        )


class TestSkeletonCustom(object):
    """
    Test for SkeletonUnitModel with callback for initialization

    """

    @pytest.fixture(scope="class")
    def skeleton_custom(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.skeleton = SkeletonUnitModel(dynamic=False)
        m.fs.skeleton.comp_list = Set(initialize=["c1", "c2"])

        # input vars for skeleton
        m.fs.skeleton.flow_comp_in = Var(
            m.fs.time, m.fs.skeleton.comp_list, initialize=1.0
        )
        m.fs.skeleton.temperature_in = Var(m.fs.time, initialize=298.15)
        m.fs.skeleton.pressure_in = Var(m.fs.time, initialize=101)

        # output vars for skeleton
        m.fs.skeleton.flow_comp_out = Var(
            m.fs.time, m.fs.skeleton.comp_list, initialize=1.0
        )
        m.fs.skeleton.temperature_out = Var(m.fs.time, initialize=298.15)
        m.fs.skeleton.pressure_out = Var(m.fs.time, initialize=101)

        # Model equations
        # Flow equation

        def rule_flow(m, t, i):
            return m.flow_comp_out[t, i] == m.flow_comp_in[t, i]

        m.fs.skeleton.eq_flow = Constraint(
            m.fs.time, m.fs.skeleton.comp_list, rule=rule_flow
        )

        m.fs.skeleton.eq_temperature = Constraint(
            expr=m.fs.skeleton.temperature_out[0]
            == m.fs.skeleton.temperature_in[0] + 10
        )
        m.fs.skeleton.eq_pressure = Constraint(
            expr=m.fs.skeleton.pressure_out[0] == m.fs.skeleton.pressure_in[0] - 10
        )

        inlet_dict = {
            "flow_mol_comp": m.fs.skeleton.flow_comp_in,
            "temperature": m.fs.skeleton.temperature_in,
            "pressure": m.fs.skeleton.pressure_in,
        }
        outlet_dict = {
            "flow_mol_comp": m.fs.skeleton.flow_comp_out,
            "temperature": m.fs.skeleton.temperature_out,
            "pressure": m.fs.skeleton.pressure_out,
        }

        m.fs.skeleton.add_ports(name="inlet", member_dict=inlet_dict)
        m.fs.skeleton.add_ports(name="outlet", member_dict=outlet_dict)

        return m

    @pytest.mark.unit
    def test_ports(self, skeleton_custom):

        # Check inlet port and port members
        assert hasattr(skeleton_custom.fs.skeleton, "inlet")
        assert len(skeleton_custom.fs.skeleton.inlet.vars) == 3
        assert hasattr(skeleton_custom.fs.skeleton.inlet, "flow_mol_comp")
        assert hasattr(skeleton_custom.fs.skeleton.inlet, "temperature")
        assert hasattr(skeleton_custom.fs.skeleton.inlet, "pressure")

        # Check outlet port and port members
        assert hasattr(skeleton_custom.fs.skeleton, "outlet")
        assert len(skeleton_custom.fs.skeleton.outlet.vars) == 3
        assert hasattr(skeleton_custom.fs.skeleton.outlet, "flow_mol_comp")
        assert hasattr(skeleton_custom.fs.skeleton.outlet, "temperature")
        assert hasattr(skeleton_custom.fs.skeleton.outlet, "pressure")

    @pytest.mark.unit
    def test_build(self, skeleton_custom):
        # Check build of model
        assert hasattr(skeleton_custom.fs.skeleton, "flow_comp_in")
        assert hasattr(skeleton_custom.fs.skeleton, "temperature_in")
        assert hasattr(skeleton_custom.fs.skeleton, "pressure_in")

        assert hasattr(skeleton_custom.fs.skeleton, "eq_flow")

        assert hasattr(skeleton_custom.fs.skeleton, "flow_comp_out")
        assert hasattr(skeleton_custom.fs.skeleton, "temperature_out")
        assert hasattr(skeleton_custom.fs.skeleton, "pressure_out")

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_exception(self, skeleton_custom):
        # Check default initialize method raises ConfigurationError
        # when config argument for custom_initializer is True but a method
        # was not set to the attribute custom_initializer
        with pytest.raises(
            ConfigurationError,
            match="Degrees of freedom is not zero during start of "
            "initialization. Fix/unfix appropriate number of variables "
            "to result in zero degrees of freedom for initialization.",
        ):
            skeleton_custom.fs.skeleton.initialize()

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_custom_initialize(self, skeleton_custom):
        # Check default initialize method raises ConfigurationError
        # when config argument for custom_initializer is True but a method
        # was not set to the attribute custom_initializer
        with pytest.raises(ConfigurationError):
            skeleton_custom.fs.skeleton.initialize()

        skeleton_custom.fs.skeleton.inlet.flow_mol_comp[0, "c1"].fix(2)
        skeleton_custom.fs.skeleton.inlet.flow_mol_comp[0, "c2"].fix(2)
        skeleton_custom.fs.skeleton.inlet.temperature.fix(325)
        skeleton_custom.fs.skeleton.inlet.pressure.fix(200)

        assert degrees_of_freedom(skeleton_custom) == 0

        def my_initialize(unit, **kwargs):
            # Callback for user provided initialization sequence
            unit.eq_temperature.deactivate()
            unit.eq_pressure.deactivate()
            unit.outlet.temperature.fix(325)
            unit.outlet.pressure.fix(200)

            if degrees_of_freedom == 0:
                solver.solve(unit)

            unit.eq_temperature.activate()
            unit.eq_pressure.activate()
            unit.outlet.temperature.unfix()
            unit.outlet.pressure.unfix()

        skeleton_custom.fs.skeleton.config.initializer = my_initialize
        skeleton_custom.fs.skeleton.initialize()

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, skeleton_custom):
        results = solver.solve(skeleton_custom)

        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, skeleton_custom):

        assert value(
            skeleton_custom.fs.skeleton.outlet.flow_mol_comp[0, "c1"]
        ) == pytest.approx(2, abs=1e-3)
        assert value(
            skeleton_custom.fs.skeleton.outlet.flow_mol_comp[0, "c2"]
        ) == pytest.approx(2, abs=1e-3)

        assert value(
            skeleton_custom.fs.skeleton.outlet.temperature[0]
        ) == pytest.approx(335, abs=1e-3)
        assert value(skeleton_custom.fs.skeleton.outlet.pressure[0]) == pytest.approx(
            190, abs=1e-3
        )
