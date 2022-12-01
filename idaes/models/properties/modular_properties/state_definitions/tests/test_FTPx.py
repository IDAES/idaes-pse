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
Tests for FTP state formulation

Authors: Andrew Lee
"""

import pytest
import numpy as np
from pytest import approx
from sys import modules

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    value,
    Var,
    units as pyunits,
)
from pyomo.util.check_units import check_units_equivalent, assert_units_consistent

# Need define_default_scaling_factors, even though it is not used directly
from idaes.models.properties.modular_properties.state_definitions.FTPx import (
    FTPx,
    define_state,
    set_metadata,
    state_initialization,
    _modified_rachford_rice,
)
from idaes.core import (
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
    declare_process_block_class,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    pass


@pytest.mark.unit
def test_set_metadata():
    assert set_metadata(None) is None


# Dummy methods for dummied submodules
class dummy_pe:
    def return_expression(b, *args):
        # Return a dummy expression for the constraint
        return b.temperature == 100


def phase_equil(b, *args):
    pass


class TestInvalidBounds(object):
    @pytest.mark.unit
    def test_bad_name(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={"p1": {"equation_of_state": DummyEoS}},
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            state_bounds={"foo": (None, None, None)},
        )

        with pytest.raises(
            ConfigurationError,
            match="props\[1\] - found unexpected state_bounds key foo. "
            "Please ensure bounds are provided only for expected state "
            "variables and that you have typed the variable names "
            "correctly.",
        ):
            m.props = m.params.build_state_block([1], defined_state=True)

    @pytest.mark.unit
    def test_mole_frac(self, caplog):
        m = ConcreteModel()

        caplog.set_level(
            idaeslog.WARNING, logger=("idaes.models.properties.modular_properties.")
        )

        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={"p1": {"equation_of_state": DummyEoS}},
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            state_bounds={"mole_frac_comp": (None, None, None)},
        )

        # Create a dummy state block
        m.props = m.params.build_state_block([1], defined_state=True)

        assert (
            "props[1] - found state_bounds argument for mole_frac_comp."
            " Mole fraction bounds are set automatically and "
            "this argument will be ignored." in caplog.text
        )


class Test1PhaseDefinedStateFalseNoBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={"p1": {"equation_of_state": DummyEoS}},
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Create state block
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 1
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value is None
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 1
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1
        assert check_units_equivalent(frame.props[1].phase_frac, None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 3
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
            frame.props[1].flow_mol
            - frame.props[1].flow_mol_phase[frame.params.phase_list.at(1)]
        )

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].mole_frac_comp[i]
                - frame.props[1].mole_frac_phase_comp[frame.params.phase_list.at(1), i]
            )

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)
        assert len(frame.props[1].sum_mole_frac_out) == 1
        assert str(frame.props[1].sum_mole_frac_out.body) == str(
            sum(
                frame.props[1].mole_frac_comp[i]
                for i in frame.props[1].params.component_list
            )
        )

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 1
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == str(
                frame.props[1].phase_frac[i]
            )

        assert_units_consistent(frame.props[1])


class Test1PhaseDefinedStateTrueWithBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={"p1": {"equation_of_state": DummyEoS}},
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            state_bounds={
                "flow_mol": (0, 100, 200),
                "temperature": (290, 345, 400),
                "pressure": (100000.0, 300000.0, 500000.0),
            },
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Create state block
        m.props = m.params.build_state_block([1], defined_state=True)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 1
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value == 100
            assert frame.props[1].flow_mol_phase[i].lb == 0
            assert frame.props[1].flow_mol_phase[i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 1
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1
        assert check_units_equivalent(frame.props[1].phase_frac, None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 3
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [("p1", "c1"), ("p1", "c2"), ("p1", "c3")]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
            frame.props[1].flow_mol
            - frame.props[1].flow_mol_phase[frame.params.phase_list.at(1)]
        )

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].mole_frac_comp[i]
                - frame.props[1].mole_frac_phase_comp[frame.params.phase_list.at(1), i]
            )

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 1
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == str(
                frame.props[1].phase_frac[i]
            )

        assert_units_consistent(frame.props[1])

    @pytest.mark.unit
    def test_initialization(self, frame):
        state_initialization(frame.props[1])
        assert isinstance(frame.props[1].temperature, Var)
        assert isinstance(frame.props[1].pressure, Var)
        assert isinstance(frame.props[1].flow_mol, Var)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert isinstance(frame.props[1].flow_mol_phase_comp, Expression)
        assert isinstance(frame.props[1].phase_frac, Var)
        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)

        assert frame.props[1].temperature.value == approx(345)
        assert frame.props[1].pressure.value == approx(3e5)
        assert frame.props[1].flow_mol.value == approx(100)
        assert frame.props[1].phase_frac["p1"].value == approx(1)
        assert frame.props[1].flow_mol_phase["p1"].value == approx(100)
        for j in frame.props[1].component_list:
            assert frame.props[1].mole_frac_comp[j].value == approx(1 / 3)
            assert frame.props[1].mole_frac_phase_comp["p1", j].value == approx(1 / 3)
            assert approx(100 / 3) == value(frame.props[1].flow_mol_phase_comp["p1", j])


class Test2PhaseDefinedStateFalseNoBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Create state block
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 2
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value is None
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 2
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1 / 2
        assert check_units_equivalent(frame.props[1].phase_frac, None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [
                ("p1", "c1"),
                ("p1", "c2"),
                ("p1", "c3"),
                ("p2", "c1"),
                ("p2", "c2"),
                ("p2", "c3"),
            ]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
            sum(
                frame.props[1].flow_mol_phase[p]
                for p in frame.props[1].params.phase_list
            )
            - frame.props[1].flow_mol
        )

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].flow_mol * frame.props[1].mole_frac_comp[i]
                - sum(
                    frame.props[1].flow_mol_phase[p]
                    * frame.props[1].mole_frac_phase_comp[p, i]
                    for p in frame.props[1].params.phase_list
                )
            )

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 1
        assert str(frame.props[1].sum_mole_frac.body) == str(
            sum(
                frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list.at(1), i
                ]
                for i in frame.props[1].params.component_list
            )
            - sum(
                frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list.at(2), i
                ]
                for i in frame.props[1].params.component_list
            )
        )

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)
        assert len(frame.props[1].sum_mole_frac_out) == 1
        assert str(frame.props[1].sum_mole_frac_out.body) == str(
            sum(
                frame.props[1].mole_frac_comp[i]
                for i in frame.props[1].params.component_list
            )
        )

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 2
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == str(
                frame.props[1].phase_frac[i] * frame.props[1].flow_mol
                - frame.props[1].flow_mol_phase[i]
            )

        assert_units_consistent(frame.props[1])


class Test2PhaseDefinedStateTrueWithBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            state_bounds={
                "flow_mol": (0, 100, 200),
                "temperature": (290, 345, 400),
                "pressure": (100000.0, 300000.0, 500000.0),
            },
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Create state block
        m.props = m.params.build_state_block([1], defined_state=True)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 2
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value == 100 / 2
            assert frame.props[1].flow_mol_phase[i].lb == 0
            assert frame.props[1].flow_mol_phase[i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 2
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1 / 2
        assert check_units_equivalent(frame.props[1].phase_frac, None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [
                ("p1", "c1"),
                ("p1", "c2"),
                ("p1", "c3"),
                ("p2", "c1"),
                ("p2", "c2"),
                ("p2", "c3"),
            ]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].total_flow_balance, Constraint)
        assert len(frame.props[1].total_flow_balance) == 1
        assert str(frame.props[1].total_flow_balance.body) == str(
            sum(
                frame.props[1].flow_mol_phase[p]
                for p in frame.props[1].params.phase_list
            )
            - frame.props[1].flow_mol
        )

        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for i in frame.props[1].component_flow_balances:
            assert i in frame.props[1].params.component_list
            assert str(frame.props[1].component_flow_balances[i].body) == str(
                frame.props[1].flow_mol * frame.props[1].mole_frac_comp[i]
                - sum(
                    frame.props[1].flow_mol_phase[p]
                    * frame.props[1].mole_frac_phase_comp[p, i]
                    for p in frame.props[1].params.phase_list
                )
            )

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 1
        assert str(frame.props[1].sum_mole_frac.body) == str(
            sum(
                frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list.at(1), i
                ]
                for i in frame.props[1].params.component_list
            )
            - sum(
                frame.props[1].mole_frac_phase_comp[
                    frame.props[1].params.phase_list.at(2), i
                ]
                for i in frame.props[1].params.component_list
            )
        )

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 2
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == str(
                frame.props[1].phase_frac[i] * frame.props[1].flow_mol
                - frame.props[1].flow_mol_phase[i]
            )

        assert_units_consistent(frame.props[1])

    @pytest.mark.unit
    def test_initialization(self, frame):
        state_initialization(frame.props[1])
        assert isinstance(frame.props[1].temperature, Var)
        assert isinstance(frame.props[1].pressure, Var)
        assert isinstance(frame.props[1].flow_mol, Var)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert isinstance(frame.props[1].flow_mol_phase_comp, Expression)
        assert isinstance(frame.props[1].phase_frac, Var)
        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)

        assert frame.props[1].temperature.value == approx(345)
        assert frame.props[1].pressure.value == approx(3e5)
        assert frame.props[1].flow_mol.value == approx(100)

        for p in frame.props[1].phase_list:
            assert frame.props[1].phase_frac[p].value == approx(0.5)
            assert frame.props[1].flow_mol_phase[p].value == approx(50)

            for j in frame.props[1].component_list:
                assert frame.props[1].mole_frac_comp[j].value == approx(1 / 3)
                assert frame.props[1].mole_frac_phase_comp[p, j].value == approx(1 / 3)
                assert approx(50 / 3) == value(frame.props[1].flow_mol_phase_comp[p, j])

        frame.props[1].phase_frac["p1"].value = 0.4
        state_initialization(frame.props[1])
        assert frame.props[1].phase_frac["p1"].value == approx(0.4)
        assert frame.props[1].flow_mol_phase["p1"].value == approx(40)
        for j in frame.props[1].component_list:
            assert frame.props[1].mole_frac_comp[j].value == approx(1 / 3)
            assert frame.props[1].mole_frac_phase_comp["p1", j].value == approx(1 / 3)
            assert approx(40 / 3) == value(frame.props[1].flow_mol_phase_comp["p1", j])

        assert frame.props[1].phase_frac["p2"].value == approx(0.5)
        assert frame.props[1].flow_mol_phase["p2"].value == approx(50)
        for j in frame.props[1].component_list:
            assert frame.props[1].mole_frac_phase_comp["p2", j].value == approx(1 / 3)
            assert approx(50 / 3) == value(frame.props[1].flow_mol_phase_comp["p2", j])
        # To avoid side effects
        frame.props[1].phase_frac["p1"].value = 0.5
        state_initialization(frame.props[1])


class Test3PhaseDefinedStateFalseNoBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
                "p3": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Create state block
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 3
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value is None
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 3
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].phase_frac, None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 9
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [
                ("p1", "c1"),
                ("p1", "c2"),
                ("p1", "c3"),
                ("p2", "c1"),
                ("p2", "c2"),
                ("p2", "c3"),
                ("p3", "c1"),
                ("p3", "c2"),
                ("p3", "c3"),
            ]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for j in frame.props[1].component_flow_balances:
            assert j in frame.params.component_list
            assert str(frame.props[1].component_flow_balances[j].body) == str(
                frame.props[1].flow_mol * frame.props[1].mole_frac_comp[j]
                - sum(
                    frame.props[1].flow_mol_phase[p]
                    * frame.props[1].mole_frac_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                )
            )

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 3
        for p in frame.props[1].sum_mole_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].sum_mole_frac[p].body) == str(
                sum(
                    frame.props[1].mole_frac_phase_comp[p, i]
                    for i in frame.props[1].params.component_list
                )
            )

        assert isinstance(frame.props[1].sum_mole_frac_out, Constraint)
        assert len(frame.props[1].sum_mole_frac_out) == 1
        assert str(frame.props[1].sum_mole_frac_out.body) == str(
            sum(
                frame.props[1].mole_frac_comp[i]
                for i in frame.props[1].params.component_list
            )
        )

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 3
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == str(
                frame.props[1].phase_frac[i] * frame.props[1].flow_mol
                - frame.props[1].flow_mol_phase[i]
            )

        assert_units_consistent(frame.props[1])


class Test3PhaseDefinedStateTrueWithBounds(object):
    # Test define_state method with no bounds and defined_State = False
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
                "p3": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            state_bounds={
                "flow_mol": (0, 100, 200),
                "temperature": (290, 345, 400),
                "pressure": (100000.0, 300000.0, 500000.0),
            },
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Create state block
        m.props = m.params.build_state_block([1], defined_state=True)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert len(frame.props[1].mole_frac_comp) == 3
        for i in frame.props[1].mole_frac_comp:
            assert i in frame.props[1].params.component_list
            assert frame.props[1].mole_frac_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_comp, None)

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert len(frame.props[1].flow_mol_phase) == 3
        for i in frame.props[1].flow_mol_phase:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].flow_mol_phase[i].value == 100 / 3
            assert frame.props[1].flow_mol_phase[i].lb == 0
            assert frame.props[1].flow_mol_phase[i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].phase_frac, Var)
        assert len(frame.props[1].phase_frac) == 3
        for i in frame.props[1].phase_frac:
            assert i in frame.props[1].params.phase_list
            assert frame.props[1].phase_frac[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].phase_frac, None)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 9
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in [
                ("p1", "c1"),
                ("p1", "c2"),
                ("p1", "c3"),
                ("p2", "c1"),
                ("p2", "c2"),
                ("p2", "c3"),
                ("p3", "c1"),
                ("p3", "c2"),
                ("p3", "c3"),
            ]
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].component_flow_balances, Constraint)
        assert len(frame.props[1].component_flow_balances) == 3
        for j in frame.props[1].component_flow_balances:
            assert j in frame.params.component_list
            assert str(frame.props[1].component_flow_balances[j].body) == str(
                frame.props[1].flow_mol * frame.props[1].mole_frac_comp[j]
                - sum(
                    frame.props[1].flow_mol_phase[p]
                    * frame.props[1].mole_frac_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                )
            )

        assert isinstance(frame.props[1].sum_mole_frac, Constraint)
        assert len(frame.props[1].sum_mole_frac) == 3
        for p in frame.props[1].sum_mole_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].sum_mole_frac[p].body) == str(
                sum(
                    frame.props[1].mole_frac_phase_comp[p, i]
                    for i in frame.props[1].params.component_list
                )
            )

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 3
        for i in frame.props[1].phase_fraction_constraint:
            assert i in frame.props[1].params.phase_list
            assert str(frame.props[1].phase_fraction_constraint[i].body) == str(
                frame.props[1].phase_frac[i] * frame.props[1].flow_mol
                - frame.props[1].flow_mol_phase[i]
            )

        assert_units_consistent(frame.props[1])

    @pytest.mark.unit
    def test_initialization(self, frame):
        state_initialization(frame.props[1])
        assert isinstance(frame.props[1].temperature, Var)
        assert isinstance(frame.props[1].pressure, Var)
        assert isinstance(frame.props[1].flow_mol, Var)

        assert isinstance(frame.props[1].mole_frac_comp, Var)
        assert isinstance(frame.props[1].flow_mol_phase, Var)
        assert isinstance(frame.props[1].flow_mol_phase_comp, Expression)
        assert isinstance(frame.props[1].phase_frac, Var)
        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)

        assert frame.props[1].temperature.value == approx(345)
        assert frame.props[1].pressure.value == approx(3e5)
        assert frame.props[1].flow_mol.value == approx(100)

        for p in frame.props[1].phase_list:
            assert frame.props[1].phase_frac[p].value == approx(1 / 3)
            assert frame.props[1].flow_mol_phase[p].value == approx(100 / 3)

            for j in frame.props[1].component_list:
                assert frame.props[1].mole_frac_comp[j].value == approx(1 / 3)
                assert frame.props[1].mole_frac_phase_comp[p, j].value == approx(1 / 3)
                assert approx(100 / 9) == value(
                    frame.props[1].flow_mol_phase_comp[p, j]
                )

        frame.props[1].phase_frac["p1"].value = 0.2
        frame.props[1].phase_frac["p2"].value = 0.5
        frame.props[1].phase_frac["p3"].value = 0.3
        state_initialization(frame.props[1])
        assert frame.props[1].phase_frac["p1"].value == approx(0.2)
        assert frame.props[1].flow_mol_phase["p1"].value == approx(20)
        for j in frame.props[1].component_list:
            assert frame.props[1].mole_frac_comp[j].value == approx(1 / 3)
            assert frame.props[1].mole_frac_phase_comp["p1", j].value == approx(1 / 3)
            assert approx(20 / 3) == value(frame.props[1].flow_mol_phase_comp["p1", j])

        assert frame.props[1].phase_frac["p2"].value == approx(0.5)
        assert frame.props[1].flow_mol_phase["p2"].value == approx(50)
        for j in frame.props[1].component_list:
            assert frame.props[1].mole_frac_phase_comp["p2", j].value == approx(1 / 3)
            assert approx(50 / 3) == value(frame.props[1].flow_mol_phase_comp["p2", j])

        assert frame.props[1].phase_frac["p3"].value == approx(0.3)
        assert frame.props[1].flow_mol_phase["p3"].value == approx(30)
        for j in frame.props[1].component_list:
            assert frame.props[1].mole_frac_comp[j].value == approx(1 / 3)
            assert frame.props[1].mole_frac_phase_comp["p3", j].value == approx(1 / 3)
            assert approx(30 / 3) == value(frame.props[1].flow_mol_phase_comp["p3", j])

        # To avoid side effects
        for p in frame.props[1].phase_list:
            frame.props[1].phase_frac[p].value = 1 / 3
        state_initialization(frame.props[1])


class TestCommon(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={"c1": {}, "c2": {}, "c3": {}},
            phases={
                "a": {"equation_of_state": DummyEoS},
                "b": {"equation_of_state": DummyEoS},
            },
            state_definition=FTPx,
            pressure_ref=100000.0,
            temperature_ref=300,
            state_bounds={
                "flow_mol": (0, 0.1, 0.2, pyunits.kmol / pyunits.s),
                "temperature": (522, 621, 720, pyunits.degR),
                "pressure": (1, 3, 5, pyunits.bar),
            },
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=True)

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
        assert check_units_equivalent(frame.props[1].flow_mol, pyunits.mol / pyunits.s)

        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

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

        assert not hasattr(frame.props[1], "sum_mole_frac_out")

        assert isinstance(frame.props[1].phase_fraction_constraint, Constraint)
        assert len(frame.props[1].phase_fraction_constraint) == 2

        assert_units_consistent(frame)

    @pytest.mark.unit
    def test_calculate_scaling_factors(self, frame):
        frame.props[1].calculate_scaling_factors()

        assert len(frame.props[1].scaling_factor) == 22
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_phase["a"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_phase["b"]] == 1e-2
        assert (
            frame.props[1].scaling_factor[frame.props[1].flow_mol_phase_comp["a", "c1"]]
            == 1e-2
        )
        assert (
            frame.props[1].scaling_factor[frame.props[1].flow_mol_phase_comp["a", "c2"]]
            == 1e-2
        )
        assert (
            frame.props[1].scaling_factor[frame.props[1].flow_mol_phase_comp["a", "c3"]]
            == 1e-2
        )
        assert (
            frame.props[1].scaling_factor[frame.props[1].flow_mol_phase_comp["b", "c1"]]
            == 1e-2
        )
        assert (
            frame.props[1].scaling_factor[frame.props[1].flow_mol_phase_comp["b", "c2"]]
            == 1e-2
        )
        assert (
            frame.props[1].scaling_factor[frame.props[1].flow_mol_phase_comp["b", "c3"]]
            == 1e-2
        )
        assert frame.props[1].scaling_factor[frame.props[1].dens_mol_phase["a"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].dens_mol_phase["b"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].mole_frac_comp["c1"]] == 1e3
        assert frame.props[1].scaling_factor[frame.props[1].mole_frac_comp["c2"]] == 1e3
        assert frame.props[1].scaling_factor[frame.props[1].mole_frac_comp["c3"]] == 1e3
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["a", "c1"]
            ]
            == 1e3
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["a", "c2"]
            ]
            == 1e3
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["a", "c3"]
            ]
            == 1e3
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["b", "c1"]
            ]
            == 1e3
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["b", "c2"]
            ]
            == 1e3
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["b", "c3"]
            ]
            == 1e3
        )
        assert frame.props[1].scaling_factor[frame.props[1].pressure] == 1e-5
        assert frame.props[1].scaling_factor[frame.props[1].temperature] == 1e-2

    # Test General Methods
    @pytest.mark.unit
    def test_get_material_flow_terms(self, frame):
        for (p, j) in frame.params._phase_component_set:
            assert str(frame.props[1].get_material_flow_terms(p, j)) == str(
                frame.props[1].flow_mol_phase_comp[p, j]
            )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, frame):
        for p in frame.params.phase_list:
            assert str(frame.props[1].get_enthalpy_flow_terms(p)) == str(
                frame.props[1]._enthalpy_flow_term[p]
            )
            assert str(frame.props[1]._enthalpy_flow_term[p].expr) == str(
                frame.props[1].flow_mol_phase[p] * frame.props[1].enth_mol_phase[p]
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, frame):
        for (p, j) in frame.params._phase_component_set:
            assert str(frame.props[1].get_material_density_terms(p, j)) == str(
                frame.props[1]._material_density_term[p, j]
            )
            assert str(frame.props[1]._material_density_term[p, j].expr) == str(
                frame.props[1].dens_mol_phase[p]
                * frame.props[1].mole_frac_phase_comp[p, j]
            )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, frame):
        for p in frame.params.phase_list:
            assert str(frame.props[1].get_energy_density_terms(p)) == str(
                frame.props[1]._energy_density_term[p]
            )
            assert str(frame.props[1]._energy_density_term[p].expr) == str(
                frame.props[1].dens_mol_phase[p]
                * frame.props[1].energy_internal_mol_phase[p]
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, frame):
        assert (
            frame.props[1].default_material_balance_type()
            == MaterialBalanceType.componentTotal
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, frame):
        assert (
            frame.props[1].default_energy_balance_type()
            == EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_get_material_flow_basis(self, frame):
        assert frame.props[1].get_material_flow_basis() == MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_define_state_vars(self, frame):
        assert frame.props[1].define_state_vars() == {
            "flow_mol": frame.props[1].flow_mol,
            "mole_frac_comp": frame.props[1].mole_frac_comp,
            "temperature": frame.props[1].temperature,
            "pressure": frame.props[1].pressure,
        }

    @pytest.mark.unit
    def test_define_display_vars(self, frame):
        assert frame.props[1].define_display_vars() == {
            "Total Molar Flowrate": frame.props[1].flow_mol,
            "Total Mole Fraction": frame.props[1].mole_frac_comp,
            "Temperature": frame.props[1].temperature,
            "Pressure": frame.props[1].pressure,
        }

    @pytest.mark.unit
    def test_conc_mol(self, frame):
        assert isinstance(frame.props[1].conc_mol_comp, Expression)
        assert len(frame.props[1].conc_mol_comp) == 3
        assert isinstance(frame.props[1].conc_mol_phase_comp, Expression)
        assert len(frame.props[1].conc_mol_phase_comp) == 6

    @pytest.mark.unit
    def test_unphysical_mol_fraction_fail(self, frame):
        frame.props[1].mole_frac_comp["c1"].value = -0.1
        with pytest.raises(
            ValueError,
            match="Component c1 has a negative mole fraction "
            "in block props\[1\]. Check your initialization.",
        ):
            frame.props[1].params.config.state_definition.state_initialization(
                frame.props[1]
            )


class TestModifiedRachfordRice(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel(name="George")
        m.component_list = ["a", "b", "c"]
        m.mole_frac_comp = Var(m.component_list, initialize=1 / len(m.component_list))
        m.K = {"a": 0.5, "b": 1, "c": 3}
        return m

    @pytest.mark.unit
    def test_flash(self, model):
        m = model
        vl_comps_list = [["a", "b", "c"], ["b", "c"], ["a", "b"], ["b"], ["c"], [], []]
        l_only_comps_list = [[], ["a"], [], ["a"], ["a", "b"], ["a", "b", "c"], []]
        v_only_comps_list = [[], [], ["c"], ["c"], [], [], ["a", "b", "c"]]
        eps_list = [None, 0, 0.01]
        # Who will validate the validation?
        assert len(vl_comps_list) == len(l_only_comps_list)
        assert len(vl_comps_list) == len(v_only_comps_list)

        expected_output = np.array(
            [
                [0.75, 0.25, 1 - 1e-5, 0.5, 1e-5, 1e-5, 1 - 1e-5],
                [0.75, 0.25, 1, 0.5, 0, 0, 1],
                [0.75, 0.25, 0.99, 0.5, 0.01, 0.01, 0.99],
            ]
        )
        for i in range(len(vl_comps_list)):
            for j in range(len(eps_list)):
                if eps_list[j] is not None:
                    vap_frac = _modified_rachford_rice(
                        m,
                        m.K,
                        vl_comps_list[i],
                        l_only_comps_list[i],
                        v_only_comps_list[i],
                        eps=eps_list[j],
                    )
                else:
                    vap_frac = _modified_rachford_rice(
                        m,
                        m.K,
                        vl_comps_list[i],
                        l_only_comps_list[i],
                        v_only_comps_list[i],
                    )
                # Convergence criterion for Newton's method is 1e-6 (because
                # we expect to pass it of to IPOPT later). We cannot expect
                # machine precision here.
                assert expected_output[j, i] == approx(vap_frac, rel=5e-5)

    @pytest.mark.unit
    def test_negative_K(self, model, caplog):
        m = model
        m.K["a"] = -1
        vap_frac = _modified_rachford_rice(m, m.K, m.component_list, [], [])
        assert vap_frac is None
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert record.levelno == idaeslog.WARNING
        assert record.getMessage() == (
            "While initializing block George, the vapor/liquid split ratio "
            "of Component a was calculated to be negative. Check the "
            "implementation of the saturation pressure, Henry's law method, "
            "or liquid density."
        )
        m.K["a"] = 0.5

    @pytest.mark.unit
    def test_unphysical_mole_fracs(self, model, caplog):
        m = model
        m.mole_frac_comp["a"] = -20
        m.mole_frac_comp["b"] = -20
        m.mole_frac_comp["c"] = -20
        vap_frac = _modified_rachford_rice(m, m.K, ["a"], ["b"], ["c"])
        assert vap_frac is None
        assert len(caplog.records) == 1
        record = caplog.records[0]
        assert record.levelno == idaeslog.WARNING
        assert record.getMessage() == (
            "Block George - phase faction initialization using "
            "modified Rachford-Rice failed. This could be "
            "because a component is essentially "
            "nonvolatile or noncondensible, or "
            "because mole fractions sum to more than "
            "one."
        )
        m.mole_frac_comp["a"] = 1 / 3
        m.mole_frac_comp["b"] = 1 / 3
        m.mole_frac_comp["c"] = 1 / 3
