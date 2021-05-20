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
Tests for FPhx state formulation.
Authors: Andrew Lee
"""

import pytest
from sys import modules

from pyomo.environ import ConcreteModel, Constraint, Var, units as pyunits
from pyomo.util.check_units import (
    check_units_equivalent, assert_units_consistent)

# Need define_default_scaling_factors, even though it is not used directly
from idaes.generic_models.properties.core.state_definitions.FPhx import \
    FPhx, define_state, set_metadata, define_default_scaling_factors
from idaes.core import (MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        declare_process_block_class)
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)
from idaes.generic_models.properties.core.generic.tests.dummy_eos import DummyEoS
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


# Mark module as all unit tests
pytestmark = pytest.mark.unit


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    pass


@pytest.mark.unit
def test_set_metadata():
    m = ConcreteModel()
    m.props = PhysicalParameterTestBlock()

    # Set metadata to make sure it is overwritten
    m.props.get_metadata().properties["enth_mol"] = {'method': "test"}

    set_metadata(m.props)

    assert m.props.get_metadata().properties["enth_mol"] == {'method': None}


class TestInvalidBounds(object):

    @pytest.mark.unit
    def test_bad_name(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K},
                "state_bounds": {"foo": (None, None, None)}})

        with pytest.raises(
                ConfigurationError,
                match="props\[1\] - found unexpected state_bounds key foo. "
                "Please ensure bounds are provided only for expected state "
                "variables and that you have typed the variable names "
                "correctly."):
            # Build state block
            m.props = m.params.build_state_block(
                [1], default={"defined_state": False})

    @pytest.mark.unit
    def test_mole_frac(self, caplog):
        m = ConcreteModel()

        caplog.set_level(
            idaeslog.WARNING,
            logger=("idaes.generic_models.properties.core."))

        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K},
                "state_bounds": {"mole_frac_comp": (None, None, None)}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": False})
        m.props[1].enth_mol_phase["p1"] = 1

        define_state(m.props[1])

        assert ("props[1] - found state_bounds argument for mole_frac_comp."
                " Mole fraction bounds are set automatically and "
                "this argument will be ignored." in caplog.text)


class Test1PhaseDefinedStateFalseNoBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": False})

        # Add necessary variables that would be built by other methods
        m.props[1].enth_mol_phase = Var(m.params.phase_list,
                                        initialize=1,
                                        units=pyunits.J/pyunits.mol)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol, Var)
        assert frame.props[1].flow_mol.value is None
        assert check_units_equivalent(frame.props[1].flow_mol,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1/3
            assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure,
                                      pyunits.Pa)

        assert isinstance(frame.props[1].enth_mol, Var)
        assert frame.props[1].enth_mol.value is None
        assert check_units_equivalent(frame.props[1].enth_mol,
                                      pyunits.J/pyunits.mol)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature,
                                      pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 1
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value is None
        assert check_units_equivalent(frame.props[1].flow_mol_phase,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 1
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1
        assert check_units_equivalent(frame.props[1].phase_frac,
                                      None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 3
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1/3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp,
                                      None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
            frame.props[1].flow_mol -
            frame.props[1].flow_mol_phase[frame.params.phase_list[1]])

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].mole_frac_comp[i] -
                frame.props[1].mole_frac_phase_comp[
                    frame.params.phase_list[1], i])

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)
        assert len(frame.props[1].sum_mole_frac_out) == 1
        assert str(frame.props[1].sum_mole_frac_out.body) == str(
                sum(frame.props[1].mole_frac_comp[i]
                    for i in frame.props[1].params.component_list))

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 1
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == \
                str(frame.props[1].phase_frac[i])

        assert_units_consistent(frame.props[1])


class Test1PhaseDefinedStateTrueWithBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "state_bounds": {"flow_mol": (0, 100, 200),
                                 "temperature": (290, 345, 400),
                                 "pressure": (1e5, 3e5, 5e5),
                                 "enth_mol": (0, 500, 1000)},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": True})

        # Add necessary variables that would be built by other methods
        m.props[1].enth_mol_phase = Var(m.params.phase_list,
                                        initialize=1,
                                        units=pyunits.J/pyunits.mol)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol, Var)
        assert frame.props[1].flow_mol.value == 100
        assert frame.props[1].flow_mol.lb == 0
        assert frame.props[1].flow_mol.ub == 200

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1/3

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5

        assert isinstance(frame.props[1].enth_mol, Var)
        assert frame.props[1].enth_mol.value == 500
        assert frame.props[1].enth_mol.lb == 0
        assert frame.props[1].enth_mol.ub == 1000

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 1
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value == 100
            assert frame.props[1].flow_mol_phase[i].lb == 0
            assert frame.props[1].flow_mol_phase[i].ub == 200

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 1
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 3
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1/3

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
            frame.props[1].flow_mol -
            frame.props[1].flow_mol_phase[frame.params.phase_list[1]])

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].mole_frac_comp[i] -
                frame.props[1].mole_frac_phase_comp[
                    frame.params.phase_list[1], i])

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 1
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == \
                str(frame.props[1].phase_frac[i])

        assert_units_consistent(frame.props[1])


class Test2PhaseDefinedStateFalseNoBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS},
                    "p2": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": False})

        # Add necessary variables that would be built by other methods
        m.props[1].enth_mol_phase = Var(m.params.phase_list,
                                        initialize=1,
                                        units=pyunits.J/pyunits.mol)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol, Var)
        assert frame.props[1].flow_mol.value is None
        assert check_units_equivalent(frame.props[1].flow_mol,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1/3
            assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure,
                                      pyunits.Pa)

        assert isinstance(frame.props[1].enth_mol, Var)
        assert frame.props[1].enth_mol.value is None
        assert check_units_equivalent(frame.props[1].enth_mol,
                                      pyunits.J/pyunits.mol)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature,
                                      pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 2
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value is None
        assert check_units_equivalent(frame.props[1].flow_mol_phase,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 2
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1/2
        assert check_units_equivalent(frame.props[1].phase_frac,
                                      None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3"),
                         ("p2", "c1"), ("p2", "c2"), ("p2", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1/3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp,
                                      None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
                sum(frame.props[1].flow_mol_phase[p]
                    for p in frame.props[1].params.phase_list) -
                frame.props[1].flow_mol)

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].flow_mol * frame.props[1].mole_frac_comp[i] -
                sum(frame.props[1].flow_mol_phase[p] *
                    frame.props[1].mole_frac_phase_comp[p, i]
                    for p in frame.props[1].params.phase_list))

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 1
        assert str(frame.props[1].sum_mole_frac.body) == str(
                sum(frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list[1], i]
                    for i in frame.props[1].params.component_list) -
                sum(frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list[2], i]
                    for i in frame.props[1].params.component_list))

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)
        assert len(frame.props[1].sum_mole_frac_out) == 1
        assert str(frame.props[1].sum_mole_frac_out.body) == str(
                sum(frame.props[1].mole_frac_comp[i]
                    for i in frame.props[1].params.component_list))

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 2
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == \
                str(frame.props[1].phase_frac[i]*frame.props[1].flow_mol -
                    frame.props[1].flow_mol_phase[i])

        assert_units_consistent(frame.props[1])


class Test2PhaseDefinedStateTrueWithBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS},
                    "p2": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "state_bounds": {"flow_mol": (0, 100, 200),
                                 "temperature": (290, 345, 400),
                                 "pressure": (1e5, 3e5, 5e5),
                                 "enth_mol": (0, 500, 1000)},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": True})

        # Add necessary variables that would be built by other methods
        m.props[1].enth_mol_phase = Var(m.params.phase_list,
                                        initialize=1,
                                        units=pyunits.J/pyunits.mol)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol, Var)
        assert frame.props[1].flow_mol.value == 100
        assert frame.props[1].flow_mol.lb == 0
        assert frame.props[1].flow_mol.ub == 200

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1/3

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5

        assert isinstance(frame.props[1].enth_mol, Var)
        assert frame.props[1].enth_mol.value == 500
        assert frame.props[1].enth_mol.lb == 0
        assert frame.props[1].enth_mol.ub == 1000

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 2
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value == 50
            assert frame.props[1].flow_mol_phase[i].lb == 0
            assert frame.props[1].flow_mol_phase[i].ub == 200

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 2
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 0.5

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3"),
                         ("p2", "c1"), ("p2", "c2"), ("p2", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1/3

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
                sum(frame.props[1].flow_mol_phase[p]
                    for p in frame.props[1].params.phase_list) -
                frame.props[1].flow_mol)

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].flow_mol * frame.props[1].mole_frac_comp[i] -
                sum(frame.props[1].flow_mol_phase[p] *
                    frame.props[1].mole_frac_phase_comp[p, i]
                    for p in frame.props[1].params.phase_list))

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 1
        assert str(frame.props[1].sum_mole_frac.body) == str(
                sum(frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list[1], i]
                    for i in frame.props[1].params.component_list) -
                sum(frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list[2], i]
                    for i in frame.props[1].params.component_list))

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 2
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == \
                str(frame.props[1].phase_frac[i]*frame.props[1].flow_mol -
                    frame.props[1].flow_mol_phase[i])

        assert_units_consistent(frame.props[1])


class Test3PhaseDefinedStateFalseNoBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS},
                    "p2": {"equation_of_state": DummyEoS},
                    "p3": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": False})

        # Add necessary variables that would be built by other methods
        m.props[1].enth_mol_phase = Var(m.params.phase_list,
                                        initialize=1,
                                        units=pyunits.J/pyunits.mol)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol, Var)
        assert frame.props[1].flow_mol.value is None
        assert check_units_equivalent(frame.props[1].flow_mol,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1/3
            assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure,
                                      pyunits.Pa)

        assert isinstance(frame.props[1].enth_mol, Var)
        assert frame.props[1].enth_mol.value is None
        assert check_units_equivalent(frame.props[1].enth_mol,
                                      pyunits.J/pyunits.mol)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature,
                                      pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 3
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value is None
        assert check_units_equivalent(frame.props[1].flow_mol_phase,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 3
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1/3
        assert check_units_equivalent(frame.props[1].phase_frac,
                                      None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 9
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3"),
                         ("p2", "c1"), ("p2", "c2"), ("p2", "c3"),
                         ("p3", "c1"), ("p3", "c2"), ("p3", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1/3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp,
                                      None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for j in frame.props[1].component_flow_balances:
            assert j in frame.params.component_list
            assert str(frame.props[1].component_flow_balances[j].body) == str(
                frame.props[1].flow_mol*frame.props[1].mole_frac_comp[j] -
                sum(frame.props[1].flow_mol_phase[p] *
                    frame.props[1].mole_frac_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list))

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 3
        for p in frame.props[1].sum_mole_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].sum_mole_frac[p].body) == str(
                    sum(frame.props[1].mole_frac_phase_comp[p, i]
                        for i in frame.props[1].params.component_list))

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)
        assert len(frame.props[1].sum_mole_frac_out) == 1
        assert str(frame.props[1].sum_mole_frac_out.body) == str(
                sum(frame.props[1].mole_frac_comp[i]
                    for i in frame.props[1].params.component_list))

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 3
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == \
                str(frame.props[1].phase_frac[i]*frame.props[1].flow_mol -
                    frame.props[1].flow_mol_phase[i])

        assert_units_consistent(frame.props[1])


class Test3PhaseDefinedStateTrueWithBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "p1": {"equation_of_state": DummyEoS},
                    "p2": {"equation_of_state": DummyEoS},
                    "p3": {"equation_of_state": DummyEoS}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "state_bounds": {"flow_mol": (0, 100, 200),
                                 "temperature": (290, 345, 400),
                                 "pressure": (1e5, 3e5, 5e5),
                                 "enth_mol": (0, 500, 1000)},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": True})

        # Add necessary variables that would be built by other methods
        m.props[1].enth_mol_phase = Var(m.params.phase_list,
                                        initialize=1,
                                        units=pyunits.J/pyunits.mol)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol, Var)
        assert frame.props[1].flow_mol.value == 100
        assert frame.props[1].flow_mol.lb == 0
        assert frame.props[1].flow_mol.ub == 200
        assert check_units_equivalent(frame.props[1].flow_mol,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1/3
            assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure,
                                      pyunits.Pa)

        assert isinstance(frame.props[1].enth_mol, Var)
        assert frame.props[1].enth_mol.value == 500
        assert frame.props[1].enth_mol.lb == 0
        assert frame.props[1].enth_mol.ub == 1000
        assert check_units_equivalent(frame.props[1].enth_mol,
                                      pyunits.J/pyunits.mol)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature,
                                      pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 3
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value == 100/3
            assert frame.props[1].flow_mol_phase[i].lb == 0
            assert frame.props[1].flow_mol_phase[i].ub == 200
        assert check_units_equivalent(frame.props[1].flow_mol_phase,
                                      pyunits.mol/pyunits.s)

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 3
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1/3
        assert check_units_equivalent(frame.props[1].phase_frac,
                                      None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 9
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3"),
                         ("p2", "c1"), ("p2", "c2"), ("p2", "c3"),
                         ("p3", "c1"), ("p3", "c2"), ("p3", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1/3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp,
                                      None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for j in frame.props[1].component_flow_balances:
            assert j in frame.params.component_list
            assert str(frame.props[1].component_flow_balances[j].body) == str(
                frame.props[1].flow_mol*frame.props[1].mole_frac_comp[j] -
                sum(frame.props[1].flow_mol_phase[p] *
                    frame.props[1].mole_frac_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list))

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 3
        for p in frame.props[1].sum_mole_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].sum_mole_frac[p].body) == str(
                    sum(frame.props[1].mole_frac_phase_comp[p, i]
                        for i in frame.props[1].params.component_list))

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 3
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == \
                str(frame.props[1].phase_frac[i]*frame.props[1].flow_mol -
                    frame.props[1].flow_mol_phase[i])

        assert_units_consistent(frame.props[1])


class TestCommon(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
                "components": {"c1": {}, "c2": {}, "c3": {}},
                "phases": {
                    "a": {"equation_of_state": DummyEoS},
                    "b": {"equation_of_state": DummyEoS}},
                "state_definition": FPhx,
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "state_bounds": {
                    "flow_mol": (0, 0.1, 0.2, pyunits.kmol/pyunits.s),
                    "temperature": (522, 621, 720, pyunits.degR),
                    "pressure": (1, 3, 5, pyunits.bar),
                    "enth_mol": (0, 0.5, 1, pyunits.kJ/pyunits.mol)},
                "base_units": {"time": pyunits.s,
                               "length": pyunits.m,
                               "mass": pyunits.kg,
                               "amount": pyunits.mol,
                               "temperature": pyunits.K}})

        # Build state block
        m.props = m.params.build_state_block(
            [1], default={"defined_state": False})

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_convert_vars(self, frame):
        # Check that all state var values and bounds were converted correctly
        assert frame.props[1].flow_mol.value == 100
        assert frame.props[1].flow_mol.lb == 0
        assert frame.props[1].flow_mol.ub == 200
        assert check_units_equivalent(frame.props[1].flow_mol,
                                      pyunits.mol/pyunits.s)

        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure,
                                      pyunits.Pa)

        assert frame.props[1].enth_mol.value == 500
        assert frame.props[1].enth_mol.lb == 0
        assert frame.props[1].enth_mol.ub == 1000
        assert check_units_equivalent(frame.props[1].enth_mol,
                                      pyunits.J/pyunits.mol)

        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature,
                                      pyunits.K)

        # Check supporting variables
        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 2

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 2

        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 1

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 2

    @pytest.mark.unit
    def test_calculate_scaling_factors(self, frame):
        frame.props[1].calculate_scaling_factors()

        assert len(frame.props[1].scaling_factor) == 23

        assert frame.props[1].scaling_factor[
            frame.props[1].dens_mol_phase["a"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].dens_mol_phase["b"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].enth_mol] == 1e-4
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase["a"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase["b"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase_comp["a", "c1"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase_comp["a", "c2"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase_comp["a", "c3"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase_comp["b", "c1"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase_comp["b", "c2"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].flow_mol_phase_comp["b", "c3"]] == 1e-2
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_comp["c1"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_comp["c2"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_comp["c3"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_phase_comp["a", "c1"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_phase_comp["a", "c2"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_phase_comp["a", "c3"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_phase_comp["b", "c1"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_phase_comp["b", "c2"]] == 1e3
        assert frame.props[1].scaling_factor[
            frame.props[1].mole_frac_phase_comp["b", "c3"]] == 1e3
        assert frame.props[1].scaling_factor[frame.props[1].pressure] == 1e-5
        assert frame.props[1].scaling_factor[
            frame.props[1].temperature] == 1e-2

    # Test General Methods
    @pytest.mark.unit
    def test_get_material_flow_terms(self, frame):
        for (p, j) in frame.params._phase_component_set:
            assert str(frame.props[1].get_material_flow_terms(p, j)) == str(
                frame.props[1].flow_mol_phase_comp[p, j])

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, frame):
        for p in frame.params.phase_list:
            assert str(frame.props[1].get_enthalpy_flow_terms(p)) == str(
                frame.props[1]._enthalpy_flow_term[p])
            assert str(frame.props[1]._enthalpy_flow_term[p].expr) == str(
                frame.props[1].flow_mol_phase[p] *
                frame.props[1].enth_mol_phase[p])

    @pytest.mark.unit
    def test_get_material_density_terms(self, frame):
        for (p, j) in frame.params._phase_component_set:
            assert str(frame.props[1].get_material_density_terms(p, j)) == str(
                frame.props[1]._material_density_term[p, j])
            assert str(
                frame.props[1]._material_density_term[p, j].expr) == str(
                frame.props[1].dens_mol_phase[p] *
                frame.props[1].mole_frac_phase_comp[p, j])

    @pytest.mark.unit
    def test_get_energy_density_terms(self, frame):
        for p in frame.params.phase_list:
            assert str(frame.props[1].get_energy_density_terms(p)) == str(
                frame.props[1]._energy_density_term[p])
            assert str(
                frame.props[1]._energy_density_term[p].expr) == str(
                frame.props[1].dens_mol_phase[p] *
                frame.props[1].energy_internal_mol_phase[p])

    @pytest.mark.unit
    def test_default_material_balance_type(self, frame):
        assert frame.props[1].default_material_balance_type() == \
            MaterialBalanceType.componentTotal

    @pytest.mark.unit
    def test_default_energy_balance_type(self, frame):
        assert frame.props[1].default_energy_balance_type() == \
            EnergyBalanceType.enthalpyTotal

    @pytest.mark.unit
    def test_get_material_flow_basis(self, frame):
        assert frame.props[1].get_material_flow_basis() == \
            MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_define_state_vars(self, frame):
        assert frame.props[1].define_state_vars() == \
            {"flow_mol": frame.props[1].flow_mol,
             "mole_frac_comp": frame.props[1].mole_frac_comp,
             "enth_mol": frame.props[1].enth_mol,
             "pressure": frame.props[1].pressure}

    @pytest.mark.unit
    def test_define_display_vars(self, frame):
        assert frame.props[1].define_display_vars() == \
            {"Total Molar Flowrate": frame.props[1].flow_mol,
             "Total Mole Fraction": frame.props[1].mole_frac_comp,
             "Molar Enthalpy": frame.props[1].enth_mol,
             "Pressure": frame.props[1].pressure}
