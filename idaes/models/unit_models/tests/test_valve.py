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
Tests for valve model
"""
import math
import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, units, value
from idaes.core import FlowsheetBlock
from idaes.models.unit_models import Valve, ValveFunctionType
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    initialization_tester,
)
from idaes.core.solvers import get_solver
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.models.properties import iapws95
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(options={"nlp_scaling_method": "user-scaling"})


class GenericValve(object):
    @pytest.fixture(scope="class")
    def valve_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = iapws95.Iapws95ParameterBlock()
        m.fs.valve = Valve(
            valve_function_callback=self.type, property_package=m.fs.properties
        )
        fin = 1000  # mol/s
        pin = 200000  # Pa
        pout = 100000  # Pa
        tin = 300  # K
        hin = iapws95.htpx(T=tin * units.K, P=pin * units.Pa)  # J/mol
        # Calculate the flow coefficient to give 1000 mol/s flow with given P
        if self.type == ValveFunctionType.linear:
            cv = 1000 / math.sqrt(pin - pout) / 0.5
        elif self.type == ValveFunctionType.quick_opening:
            cv = 1000 / math.sqrt(pin - pout) / math.sqrt(0.5)
        elif self.type == ValveFunctionType.equal_percentage:
            cv = 1000 / math.sqrt(pin - pout) / 100 ** (0.5 - 1)
        # set inlet
        m.fs.valve.inlet.enth_mol[0].fix(hin)
        m.fs.valve.inlet.flow_mol[0].fix(fin)
        m.fs.valve.inlet.pressure[0].fix(pin)
        m.fs.valve.outlet.pressure[0].set_value(pout)
        m.fs.valve.Cv.fix(cv)
        m.fs.valve.valve_opening.fix(0.5)
        iscale.calculate_scaling_factors(m)
        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, valve_model):
        assert hasattr(valve_model.fs.valve, "inlet")
        assert len(valve_model.fs.valve.inlet.vars) == 3
        assert hasattr(valve_model.fs.valve.inlet, "flow_mol")
        assert hasattr(valve_model.fs.valve.inlet, "enth_mol")
        assert hasattr(valve_model.fs.valve.inlet, "pressure")

        assert hasattr(valve_model.fs.valve, "outlet")
        assert len(valve_model.fs.valve.outlet.vars) == 3
        assert hasattr(valve_model.fs.valve.outlet, "flow_mol")
        assert hasattr(valve_model.fs.valve.outlet, "enth_mol")
        assert hasattr(valve_model.fs.valve.outlet, "pressure")

        assert hasattr(valve_model.fs.valve, "pressure_flow_equation")
        assert hasattr(valve_model.fs.valve, "valve_opening")
        assert hasattr(valve_model.fs.valve, "valve_function")
        assert hasattr(valve_model.fs.valve, "flow_var")

        if self.type == ValveFunctionType.equal_percentage:
            # alpha valve function parameter is one more than others
            assert number_variables(valve_model) == 12
        else:
            assert number_variables(valve_model) == 11

        assert number_total_constraints(valve_model) == 6
        assert number_unused_variables(valve_model) == 0

    @pytest.mark.component
    def test_units(self, valve_model):
        assert_units_consistent(valve_model)
        assert_units_equivalent(valve_model.fs.valve.flow_var[0], units.mol / units.s)

    @pytest.mark.unit
    def test_dof(self, valve_model):
        assert degrees_of_freedom(valve_model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, valve_model):
        initialization_tester(valve_model, unit=valve_model.fs.valve)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, valve_model):
        results = solver.solve(valve_model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, valve_model):
        # calculated Cv to yeild this solution
        assert pytest.approx(1000, rel=1e-3) == value(
            valve_model.fs.valve.outlet.flow_mol[0]
        )


class TestLinearValve(GenericValve):
    type = ValveFunctionType.linear


class TestEPValve(GenericValve):
    type = ValveFunctionType.equal_percentage


class TestQuickValve(GenericValve):
    type = ValveFunctionType.quick_opening
