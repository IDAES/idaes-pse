#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for valve model
"""

import math
import pytest
import re

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    units,
    value,
)
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    initialization_tester,
)
from idaes.core.scaling.util import jacobian_cond
from idaes.core.solvers import get_solver
from idaes.core.util import DiagnosticsToolbox
from idaes.core.util.exceptions import ConfigurationError

from idaes.models.properties import iapws95
from idaes.models.properties.general_helmholtz import helmholtz_available
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmValve,
    ValveFunctionType,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(solver="ipopt_v2")


@pytest.mark.parametrize(
    "phase",
    ["Liq", "Vap"],
    scope="class",
)
@pytest.mark.parametrize(
    "callback_type",
    [
        ValveFunctionType.linear,
        ValveFunctionType.quick_opening,
        ValveFunctionType.equal_percentage,
        ValveFunctionType.custom,
    ],
    scope="class",
)
@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
class TestHelmholtzValve(object):
    @pytest.fixture(scope="class")
    @classmethod
    def valve_model(cls, phase, callback_type):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # Build parameter block with the correct phase
        if phase == "Liq":
            pp = iapws95.PhaseType.L
        elif phase == "Vap":
            pp = iapws95.PhaseType.G
        else:
            raise AssertionError()
        m.fs.properties = iapws95.Iapws95ParameterBlock(phase_presentation=pp)

        if callback_type is ValveFunctionType.custom:

            def valve_function_callback(blk):
                @blk.Expression(blk.flowsheet().time)
                def valve_function(b, t):
                    return b.valve_opening[t] ** 2

        else:
            valve_function_callback = None

        m.fs.valve = HelmValve(
            valve_function=callback_type,
            valve_function_callback=valve_function_callback,
            property_package=m.fs.properties,
            phase=phase,
        )
        fin = 1000  # mol/s
        pin = 200000  # Pa
        pout = 100000  # Pa
        if phase == "Liq":
            tin = 300  # K
        elif phase == "Vap":
            tin = 400  # K
        else:
            raise AssertionError()
        hin = iapws95.htpx(T=tin * units.K, P=pin * units.Pa)  # J/mol
        # Calculate the flow coefficient to give 1000 mol/s flow with given P
        if phase == "Liq":
            dP = pin - pout
        elif phase == "Vap":
            dP = pin**2 - pout**2
        else:
            raise AssertionError()
        if callback_type == ValveFunctionType.linear:
            cv = 1000 / math.sqrt(dP) / 0.5
        elif callback_type == ValveFunctionType.quick_opening:
            cv = 1000 / math.sqrt(dP) / math.sqrt(0.5)
        elif callback_type == ValveFunctionType.equal_percentage:
            cv = 1000 / math.sqrt(dP) / 100 ** (0.5 - 1)
        elif callback_type == ValveFunctionType.custom:
            cv = 1000 / math.sqrt(dP) / 0.5**2

        # Set inlet conditions
        m.fs.valve.inlet.enth_mol[0].fix(hin)
        m.fs.valve.inlet.flow_mol[0].fix(fin)
        m.fs.valve.inlet.pressure[0].fix(pin)
        m.fs.valve.outlet.pressure[0].set_value(pout)
        m.fs.valve.Cv.fix(cv)
        m.fs.valve.valve_opening.fix(0.5)

        # Scale valves
        helmholtz_scaler = m.fs.properties.default_state_scaler_class()
        helmholtz_scaler.default_scaling_factors["flow_mol"] = 1e-3
        m.fs.properties.default_state_scaler_object = helmholtz_scaler

        valve_scaler = m.fs.valve.default_scaler()
        valve_scaler.scale_model(m.fs.valve)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, valve_model, callback_type):
        unit = valve_model.fs.valve

        assert hasattr(unit, "inlet")
        assert len(unit.inlet.vars) == 3
        assert hasattr(unit.inlet, "flow_mol")
        assert hasattr(unit.inlet, "enth_mol")
        assert hasattr(unit.inlet, "pressure")

        assert hasattr(unit, "outlet")
        assert len(unit.outlet.vars) == 3
        assert hasattr(unit.outlet, "flow_mol")
        assert hasattr(unit.outlet, "enth_mol")
        assert hasattr(unit.outlet, "pressure")

        assert hasattr(unit, "pressure_flow_equation")
        assert hasattr(unit, "valve_opening")
        assert hasattr(unit, "valve_function")

        if callback_type is ValveFunctionType.equal_percentage:
            # alpha valve function parameter is one more than others
            assert number_variables(unit) == 10
        else:
            assert number_variables(unit) == 9

        assert number_total_constraints(unit) == 4
        assert number_unused_variables(unit) == 0

    @pytest.mark.component
    def test_structural_issues(self, valve_model):
        dt = DiagnosticsToolbox(valve_model.fs.valve)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, valve_model):
        initialization_tester(valve_model, unit=valve_model.fs.valve, solver="ipopt_v2")

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, valve_model):
        results = solver.solve(valve_model, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, valve_model):
        # calculated Cv to yield this solution
        assert pytest.approx(1000, rel=1e-3) == value(
            valve_model.fs.valve.outlet.flow_mol[0]
        )
        assert pytest.approx(100000, rel=1e-3) == value(
            valve_model.fs.valve.outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, valve_model, phase, callback_type):
        dt = DiagnosticsToolbox(valve_model)
        dt.assert_no_numerical_warnings()
        expected_condition_numbers = {
            # phase, callback_type: (unscaled, scaled)
            ("Liq", ValveFunctionType.linear): (5612.48, 234.796),
            ("Vap", ValveFunctionType.linear): (2.40365e6, 257.844),
            ("Liq", ValveFunctionType.quick_opening): (5612.48, 234.796),
            ("Vap", ValveFunctionType.quick_opening): (2.40365e6, 257.844),
            ("Liq", ValveFunctionType.equal_percentage): (5612.48, 234.796),
            ("Vap", ValveFunctionType.equal_percentage): (2.40365e6, 257.844),
            ("Liq", ValveFunctionType.custom): (5612.48, 234.796),
            ("Vap", ValveFunctionType.custom): (2.40365e6, 257.844),
        }
        assert (phase, callback_type) in expected_condition_numbers
        unscaled, scaled = expected_condition_numbers[(phase, callback_type)]
        assert jacobian_cond(valve_model, scaled=False) == pytest.approx(
            unscaled, rel=1e-3
        )
        assert jacobian_cond(valve_model, scaled=True) == pytest.approx(
            scaled, rel=1e-3
        )


@pytest.mark.unit
def test_custom_callback_exceptions():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()

    def valve_function_callback(blk):
        @blk.Expression(blk.flowsheet().time)
        def valve_function(b, t):
            return b.valve_opening[t] ** 2

    with pytest.raises(
        ConfigurationError,
        match=re.escape(
            "A valve function callback was provided but the valve "
            "function type is not custom."
        ),
    ):
        m.fs.valve = HelmValve(
            valve_function_callback=valve_function_callback,
            property_package=m.fs.properties,
            phase="Liq",
        )
    with pytest.raises(
        ConfigurationError,
        match=re.escape("No custom valve function callback provided"),
    ):
        m.fs.valve = HelmValve(
            valve_function=ValveFunctionType.custom,
            property_package=m.fs.properties,
            phase="Liq",
        )
