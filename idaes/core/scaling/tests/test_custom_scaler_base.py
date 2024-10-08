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
"""
Tests for CustomScalerBase.

Author: Andrew Lee
"""
import pytest
import re

from pyomo.environ import (
    Block,
    ComponentMap,
    ConcreteModel,
    Constraint,
    Set,
    Suffix,
    units,
    value,
    Var,
)
from pyomo.common.config import ConfigDict

from idaes.core.scaling.custom_scaler_base import (
    CustomScalerBase,
    ConstraintScalingScheme,
)
from idaes.core.util.constants import Constants
from idaes.core.util.testing import PhysicalParameterTestBlock
import idaes.logger as idaeslog


class DummyScaler(CustomScalerBase):
    """
    Create some dummy methods that will record that they ran in lists
    """

    def variable_scaling_routine(self, model, overwrite, submodel_scalers):
        # Create some lists to show that each method was run
        model._verification = []
        model._subscalers = {}
        model.overwrite = []

        model._verification.append("variables")
        model._subscalers["variables"] = submodel_scalers
        model.overwrite.append(overwrite)

    def constraint_scaling_routine(self, model, overwrite, submodel_scalers):
        model._verification.append("constraints")
        model._subscalers["constraints"] = submodel_scalers
        model.overwrite.append(overwrite)

    def fill_in_1(self, model):
        model._verification.append("fill_in_1")

    def fill_in_2(self, model):
        model._verification.append("fill_in_2")

    def dummy_method(self, model, overwrite):
        model._dummy_scaler_test = overwrite


@pytest.fixture
def model():
    m = ConcreteModel()

    m.pressure = Var(units=units.Pa, bounds=(1e4, None))
    m.temperature = Var(units=units.K, bounds=(250, 350))
    m.volume_mol = Var(units=units.m**3 / units.mol, bounds=(None, 10))
    m.enth_mol = Var(units=units.J / units.mol)

    m.ideal_gas = Constraint(
        expr=m.pressure * m.volume_mol == Constants.gas_constant * m.temperature
    )
    m.enthalpy_eq = Constraint(
        expr=4.81 * units.J / units.mol / units.K * m.temperature == m.enth_mol
    )

    m.scaling_factor = Suffix(direction=Suffix.EXPORT)

    return m


class TestCustomScalerBase:
    @pytest.mark.unit
    def test_init(self):
        sb = CustomScalerBase()

        assert isinstance(sb.config, ConfigDict)
        assert sb.config.zero_tolerance == 1e-12
        assert sb.config.max_variable_scaling_factor == 1e10
        assert sb.config.min_variable_scaling_factor == 1e-10
        assert sb.config.max_constraint_scaling_factor == 1e10
        assert sb.config.min_constraint_scaling_factor == 1e-10
        assert not sb.config.overwrite

        assert sb.default_scaling_factors == {}
        assert sb.unit_scaling_factors == {
            "Temperature": (units.K, 1e-2),
            "Pressure": (units.Pa, 1e-5),
        }

    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        sb = CustomScalerBase()

        with pytest.raises(
            NotImplementedError,
            match="Custom Scaler has not implemented a variable_scaling_routine method.",
        ):
            sb.variable_scaling_routine("foo")

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        sb = CustomScalerBase()

        with pytest.raises(
            NotImplementedError,
            match="Custom Scaler has not implemented a constraint_scaling_routine method.",
        ):
            sb.constraint_scaling_routine("foo")

    @pytest.mark.unit
    def test_scale_model(self):
        # Dummy object to hold testing data
        model = ConcreteModel()

        sb = DummyScaler()
        sb.scale_model(
            model,
            first_stage_fill_in=[sb.fill_in_1],
            second_stage_fill_in=[sb.fill_in_2],
            submodel_scalers={"foo", "bar"},
        )

        # Check order methods were run in
        assert model._verification == [
            "variables",
            "fill_in_1",
            "constraints",
            "fill_in_2",
        ]

        # Check that overwrite was passed on
        model.overwrite == [False, False]

        # Check that submodel scalers were passed on to methods
        assert model._subscalers == {
            "variables": {"foo", "bar"},
            "constraints": {"foo", "bar"},
        }

    @pytest.mark.unit
    def test_get_default_scaling_factor(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # No defaults defined yet
        assert sb.get_default_scaling_factor(model.pressure) is None
        assert "No default scaling factor found for pressure" in caplog.text

        # Set a default
        sb.default_scaling_factors["pressure"] = 1e-4
        assert sb.get_default_scaling_factor(model.pressure) == 1e-4

    @pytest.mark.unit
    def test_get_default_scaling_factor_indexed(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        m = ConcreteModel()
        m.v = Var([1, 2, 3, 4])

        sb = CustomScalerBase()

        # No defaults defined yet
        assert sb.get_default_scaling_factor(m.v[1]) is None
        assert "No default scaling factor found for v[1]" in caplog.text

        # Set a default for the indexed var
        sb.default_scaling_factors["v"] = 1e-4
        assert sb.get_default_scaling_factor(m.v[1]) == 1e-4

        # Set a default for the specific element
        sb.default_scaling_factors["v[1]"] = 1e-8
        assert sb.get_default_scaling_factor(m.v[1]) == 1e-8

    @pytest.mark.unit
    def test_scale_variable_by_component(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # No scaling factors set
        sb.scale_variable_by_component(model.pressure, model.temperature)
        assert sb.get_scaling_factor(model.pressure) is None
        assert (
            "Could not set scaling factor for pressure, no scaling factor set for temperature"
            in caplog.text
        )

        # Set a scaling factor for temperature
        model.scaling_factor = Suffix(direction=Suffix.EXPORT)
        model.scaling_factor[model.temperature] = 1e-2

        sb.scale_variable_by_component(model.pressure, model.temperature)
        assert model.scaling_factor[model.pressure] == 1e-2

        # Change temperature scaling and check for overwrite
        model.scaling_factor[model.temperature] = 1e-3

        sb.scale_variable_by_component(
            model.pressure, model.temperature, overwrite=False
        )
        assert model.scaling_factor[model.pressure] == 1e-2

        sb.scale_variable_by_component(
            model.pressure, model.temperature, overwrite=True
        )
        assert model.scaling_factor[model.pressure] == 1e-3

    @pytest.mark.unit
    def test_scale_variable_by_bounds(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # Set a scaling factor for temperature to test overwrite
        model.scaling_factor[model.temperature] = 1e-2

        sb.scale_variable_by_bounds(model.temperature, overwrite=False)
        assert model.scaling_factor[model.temperature] == 1e-2

        # Both bounds
        sb.scale_variable_by_bounds(model.temperature, overwrite=True)
        assert model.scaling_factor[model.temperature] == 1 / 300

        # Lower bound only
        sb.scale_variable_by_bounds(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-4

        # Upper bound only
        sb.scale_variable_by_bounds(model.volume_mol)
        assert model.scaling_factor[model.volume_mol] == 0.1

        # No bounds
        sb.scale_variable_by_bounds(model.enth_mol)
        assert (
            "No scaling factor set for enth_mol; variable has no bounds." in caplog.text
        )
        assert model.enth_mol not in model.scaling_factor

        # Test case where bounds give a magnitude of 0
        model.pressure.setlb(-1e4)
        model.pressure.setub(1e4)
        sb.scale_variable_by_bounds(model.pressure, overwrite=True)
        assert model.scaling_factor[model.pressure] == 1

    @pytest.mark.unit
    def test_scale_variable_by_default(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # No defaults defined yet
        sb.scale_variable_by_default(model.pressure)
        assert model.pressure not in model.scaling_factor
        assert (
            "Could not set scaling factor for pressure, no default scaling factor set."
            in caplog.text
        )

        # Set a default
        sb.default_scaling_factors["pressure"] = 1e-4
        sb.scale_variable_by_default(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-4

        # Change default to check overwrite
        sb.default_scaling_factors["pressure"] = 1e-5
        sb.scale_variable_by_default(model.pressure, overwrite=False)
        assert model.scaling_factor[model.pressure] == 1e-4
        sb.scale_variable_by_default(model.pressure, overwrite=True)
        assert model.scaling_factor[model.pressure] == 1e-5

    @pytest.mark.unit
    def test_scale_variable_by_units(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # Units in dict, no conversion required
        sb.scale_variable_by_units(model.pressure)
        assert model.scaling_factor[model.pressure] == 1e-5

        # Equivalent units in dict, need conversion
        model.p2 = Var(units=units.bar)
        sb.scale_variable_by_units(model.p2)
        assert model.scaling_factor[model.p2] == pytest.approx(1, rel=1e-8)

        # No match - no sf assigned
        sb.scale_variable_by_units(model.enth_mol)
        assert model.enth_mol not in model.scaling_factor
        assert (
            "No scaling factor set for enth_mol; no match for units J/mol found "
            "in self.unit_scaling_factors" in caplog.text
        )

        # Test for overwrite
        model.scaling_factor[model.temperature] = 42
        sb.scale_variable_by_units(model.temperature, overwrite=False)
        assert model.scaling_factor[model.temperature] == 42
        sb.scale_variable_by_units(model.temperature, overwrite=True)
        assert model.scaling_factor[model.temperature] == 1e-2

    @pytest.mark.unit
    def test_scale_constraint_by_default(self, model, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")
        sb = CustomScalerBase()

        # No defaults defined yet
        sb.scale_constraint_by_default(model.ideal_gas)
        assert model.ideal_gas not in model.scaling_factor
        assert (
            "Could not set scaling factor for ideal_gas, no default scaling factor set."
            in caplog.text
        )

        # Set a default
        sb.default_scaling_factors["ideal_gas"] = 1e-3
        sb.scale_constraint_by_default(model.ideal_gas)
        assert model.scaling_factor[model.ideal_gas] == 1e-3

        # Change default to check overwrite
        sb.default_scaling_factors["ideal_gas"] = 1e-5
        sb.scale_constraint_by_default(model.ideal_gas, overwrite=False)
        assert model.scaling_factor[model.ideal_gas] == 1e-3
        sb.scale_constraint_by_default(model.ideal_gas, overwrite=True)
        assert model.scaling_factor[model.ideal_gas] == 1e-5

    @pytest.mark.unit
    def test_get_expression_nominal_values(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1

        nominal_values = sb.get_expression_nominal_values(model.ideal_gas.expr)

        # Nominal values will be (R*T, P*V)
        assert nominal_values == [
            pytest.approx(831.446, rel=1e-5),
            pytest.approx(1e6, rel=1e-5),
        ]

        # Check redirection for ConstraintData objects
        nominal_values = sb.get_expression_nominal_values(model.ideal_gas)
        assert nominal_values == [
            pytest.approx(831.446, rel=1e-5),
            pytest.approx(1e6, rel=1e-5),
        ]

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_harmonic(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas,
            scheme=ConstraintScalingScheme.harmonicMean,
            overwrite=False,
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.harmonicMean, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            (1 / 831.446 + 1e-6), rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_max(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(model.ideal_gas, overwrite=False)
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(model.ideal_gas, overwrite=True)
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(1e-6, rel=1e-5)

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_min(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas,
            scheme=ConstraintScalingScheme.inverseMinimum,
            overwrite=False,
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas,
            scheme=ConstraintScalingScheme.inverseMinimum,
            overwrite=True,
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / 831.446, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_inv_sum(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseSum, overwrite=False
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseSum, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / (831.446 + 1e6), rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_rss(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseRSS, overwrite=False
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_value(
            model.ideal_gas, scheme=ConstraintScalingScheme.inverseRSS, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / (831.446**2 + 1e6**2) ** 0.5, rel=1e-5
        )

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_value_invalid_scheme(self, model):
        sb = CustomScalerBase()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Invalid value for 'scheme' argument (foo) in "
                "scale_constraint_by_nominal_value."
            ),
        ):
            sb.scale_constraint_by_nominal_value(model.ideal_gas, scheme="foo")

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_derivative_1norm(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_derivative_norm(
            model.ideal_gas, norm=1, overwrite=False
        )
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_derivative_norm(
            model.ideal_gas, norm=1, overwrite=True
        )
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            4.99792e-7, rel=1e-5  # (1/(8.314+1e6+1e6)
        )

        # Check for clean up
        assert model.pressure.value is None
        assert model.temperature.value is None
        assert model.pressure.value is None

    @pytest.mark.unit
    def test_scale_constraint_by_nominal_derivative_2norm(self, model):
        sb = CustomScalerBase()

        # Set variable scaling factors for testing
        model.scaling_factor[model.pressure] = 1e-5
        model.scaling_factor[model.temperature] = 1e-2
        model.scaling_factor[model.volume_mol] = 1e-1
        model.scaling_factor[model.ideal_gas] = 1

        # overwrite = False, no change
        sb.scale_constraint_by_nominal_derivative_norm(model.ideal_gas, overwrite=False)
        assert model.scaling_factor[model.ideal_gas] == 1
        # overwrite = True
        sb.scale_constraint_by_nominal_derivative_norm(model.ideal_gas, overwrite=True)
        assert model.scaling_factor[model.ideal_gas] == pytest.approx(
            1 / (8.314**2 + 1e6**2 + 1e6**2) ** 0.5, rel=1e-5
        )

        # Check for clean up
        assert model.pressure.value is None
        assert model.temperature.value is None
        assert model.pressure.value is None

    @pytest.mark.unit
    def test_propagate_state_scaling(self):
        # Dummy up two state blocks
        m = ConcreteModel()

        m.properties = PhysicalParameterTestBlock()

        m.state1 = m.properties.build_state_block([1, 2, 3])
        m.state2 = m.properties.build_state_block([1, 2, 3])

        # Set scaling factors on state1
        for t, sd in m.state1.items():
            sd.scaling_factor = Suffix(direction=Suffix.EXPORT)
            sd.scaling_factor[sd.temperature] = 100 * t
            sd.scaling_factor[sd.pressure] = 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 10 * t * count
                count += 1

        sb = CustomScalerBase()
        sb.propagate_state_scaling(m.state2, m.state1)

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 100 * t
            assert sd.scaling_factor[sd.pressure] == 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * t * count
                count += 1

        # Test for overwrite
        for t, sd in m.state1.items():
            sd.scaling_factor[sd.temperature] = 200 * t
            sd.scaling_factor[sd.pressure] = 2e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                sd.scaling_factor[j] = 20 * t * count
                count += 1

        sb.propagate_state_scaling(m.state2, m.state1, overwrite=False)

        for t, sd in m.state2.items():
            assert sd.scaling_factor[sd.temperature] == 100 * t
            assert sd.scaling_factor[sd.pressure] == 1e5 * t

            count = 1
            for j in sd.flow_mol_phase_comp.values():
                assert sd.scaling_factor[j] == 10 * t * count
                count += 1

    @pytest.mark.unit
    def test_call_submodel_scaler_method_no_scaler(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(m, "b", method="dummy_method", overwrite=True)

        for bd in m.b.values():
            assert not hasattr(bd, "_dummy_scaler_test")

        assert (
            "No default Scaler set for unknown.b. Cannot call dummy_method."
            in caplog.text
        )

    @pytest.mark.unit
    def test_call_submodel_scaler_method_default_scaler(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])
        for bd in m.b.values():
            bd.default_scaler = DummyScaler

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(m, "b", method="dummy_method", overwrite=True)

        for bd in m.b.values():
            assert bd._dummy_scaler_test

        assert "Using default Scaler for unknown.b." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_user_scaler(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        scaler_map = ComponentMap()
        scaler_map[m.b] = DummyScaler()

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(
            m,
            "b",
            method="dummy_method",
            submodel_scalers=scaler_map,
            overwrite=False,
        )

        for bd in m.b.values():
            assert not bd._dummy_scaler_test

        assert "Using user-defined Scaler for unknown.b." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_user_scaler_class(self, caplog):
        caplog.set_level(idaeslog.DEBUG, logger="idaes")

        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        scaler_map = ComponentMap()
        scaler_map[m.b] = DummyScaler()

        sb = CustomScalerBase()
        sb.call_submodel_scaler_method(
            m,
            "b",
            method="dummy_method",
            submodel_scalers=scaler_map,
            overwrite=False,
        )

        for bd in m.b.values():
            assert not bd._dummy_scaler_test

        assert "Using user-defined Scaler for unknown.b." in caplog.text

    @pytest.mark.unit
    def test_call_submodel_scaler_method_invalid_method(self):
        # Dummy up a nested model
        m = ConcreteModel()
        m.b = Block([1, 2, 3])

        scaler_map = ComponentMap()
        scaler_map[m.b] = DummyScaler()

        sb = CustomScalerBase()

        with pytest.raises(
            AttributeError,
            match="Scaler for unknown.b does not have a method named foo.",
        ):
            sb.call_submodel_scaler_method(
                m,
                "b",
                method="foo",
                submodel_scalers=scaler_map,
                overwrite=False,
            )
