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
Tests for FpcTP state formulation

Authors: Andrew Lee
"""

import pytest
from sys import modules

from pyomo.environ import ConcreteModel, Constraint, Expression, Var, units as pyunits
from pyomo.util.check_units import check_units_equivalent, assert_units_consistent

# Need define_default_scaling_factors, even though it is not used directly
from idaes.models.properties.modular_properties.state_definitions.FpcTP import (
    FpcTP,
    define_state,
    set_metadata,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
    declare_process_block_class,
    Solvent,
    Solute,
    AqueousPhase,
    VaporPhase,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS
from idaes.models.properties.modular_properties.pure.Perrys import Perrys
from idaes.models.properties.modular_properties.pure.NIST import NIST
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import degrees_of_freedom, large_residuals_set


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    pass


@pytest.mark.unit
def test_set_metadata():
    assert set_metadata(None) is None


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
            # Build state block
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

        with pytest.raises(
            ConfigurationError,
            match="props\[1\] - found unexpected state_bounds key "
            "mole_frac_comp. Please ensure bounds are provided only for "
            "expected state variables and that you have typed the "
            "variable names correctly.",
        ):
            # Build state block
            m.props = m.params.build_state_block([1], defined_state=False)


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

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert not frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol_phase_comp, Var)
        assert len(frame.props[1].flow_mol_phase_comp) == 3
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert (p, i) in frame.props[1].params._phase_component_set
            assert frame.props[1].flow_mol_phase_comp[p, i].value is None
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 3
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in frame.params._phase_component_set
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_expressions(self, frame):
        assert isinstance(frame.props[1].flow_mol, Expression)
        assert len(frame.props[1].flow_mol) == 1
        assert str(frame.props[1].flow_mol.expr) == str(
            sum(
                frame.props[1].flow_mol_phase_comp[i]
                for i in frame.props[1].params._phase_component_set
            )
        )

        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 1
        for p in frame.props[1].flow_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].flow_mol_phase[p].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for j in frame.props[1].params.component_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].flow_mol_comp, Expression)
        assert len(frame.props[1].flow_mol_comp) == 3
        for j in frame.props[1].flow_mol_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].flow_mol_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].mole_frac_comp, Expression)
        assert len(frame.props[1].mole_frac_comp) == 3
        for j in frame.props[1].mole_frac_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].mole_frac_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
                / frame.props[1].flow_mol
            )

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 1
        for p in frame.props[1].phase_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].phase_frac[p].expr) == str(1.0)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 3
        for (p, i) in frame.props[1].mole_frac_phase_comp_eq:
            assert str(frame.props[1].mole_frac_phase_comp_eq[p, i].body) == str(
                frame.props[1].mole_frac_phase_comp[p, i]
                * frame.props[1].flow_mol_phase[p]
                - frame.props[1].flow_mol_phase_comp[p, i]
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
                "flow_mol_phase_comp": (0, 100, 200),
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

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=True)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert not frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol_phase_comp, Var)
        assert len(frame.props[1].flow_mol_phase_comp) == 3
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert (p, i) in frame.props[1].params._phase_component_set
            assert frame.props[1].flow_mol_phase_comp[p, i].value == 100
            assert frame.props[1].flow_mol_phase_comp[p, i].lb == 0
            assert frame.props[1].flow_mol_phase_comp[p, i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

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

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 3
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in frame.params._phase_component_set
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_expressions(self, frame):
        assert isinstance(frame.props[1].flow_mol, Expression)
        assert len(frame.props[1].flow_mol) == 1
        assert str(frame.props[1].flow_mol.expr) == str(
            sum(
                frame.props[1].flow_mol_phase_comp[i]
                for i in frame.props[1].params._phase_component_set
            )
        )

        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 1
        for p in frame.props[1].flow_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].flow_mol_phase[p].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for j in frame.props[1].params.component_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].flow_mol_comp, Expression)
        assert len(frame.props[1].flow_mol_comp) == 3
        for j in frame.props[1].flow_mol_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].flow_mol_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].mole_frac_comp, Expression)
        assert len(frame.props[1].mole_frac_comp) == 3
        for j in frame.props[1].mole_frac_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].mole_frac_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
                / frame.props[1].flow_mol
            )

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 1
        for p in frame.props[1].phase_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].phase_frac[p].expr) == str(1.0)

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 3
        for (p, i) in frame.props[1].mole_frac_phase_comp_eq:
            assert str(frame.props[1].mole_frac_phase_comp_eq[p, i].body) == str(
                frame.props[1].mole_frac_phase_comp[p, i]
                * frame.props[1].flow_mol_phase[p]
                - frame.props[1].flow_mol_phase_comp[p, i]
            )

        assert_units_consistent(frame.props[1])


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

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert not frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol_phase_comp, Var)
        assert len(frame.props[1].flow_mol_phase_comp) == 6
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert (p, i) in frame.props[1].params._phase_component_set
            assert frame.props[1].flow_mol_phase_comp[p, i].value is None
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in frame.params._phase_component_set
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_expressions(self, frame):
        assert isinstance(frame.props[1].flow_mol, Expression)
        assert len(frame.props[1].flow_mol) == 1
        assert str(frame.props[1].flow_mol.expr) == str(
            sum(
                frame.props[1].flow_mol_phase_comp[i]
                for i in frame.props[1].params._phase_component_set
            )
        )

        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 2
        for p in frame.props[1].flow_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].flow_mol_phase[p].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for j in frame.props[1].params.component_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].flow_mol_comp, Expression)
        assert len(frame.props[1].flow_mol_comp) == 3
        for j in frame.props[1].flow_mol_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].flow_mol_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].mole_frac_comp, Expression)
        assert len(frame.props[1].mole_frac_comp) == 3
        for j in frame.props[1].mole_frac_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].mole_frac_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
                / frame.props[1].flow_mol
            )

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 2
        for p in frame.props[1].phase_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].phase_frac[p].expr) == str(
                frame.props[1].flow_mol_phase[p] / frame.props[1].flow_mol
            )

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 6
        for (p, i) in frame.props[1].mole_frac_phase_comp_eq:
            assert str(frame.props[1].mole_frac_phase_comp_eq[p, i].body) == str(
                frame.props[1].mole_frac_phase_comp[p, i]
                * frame.props[1].flow_mol_phase[p]
                - frame.props[1].flow_mol_phase_comp[p, i]
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
                "flow_mol_phase_comp": (0, 100, 200),
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

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=True)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert not frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol_phase_comp, Var)
        assert len(frame.props[1].flow_mol_phase_comp) == 6
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert (p, i) in frame.props[1].params._phase_component_set
            assert frame.props[1].flow_mol_phase_comp[p, i].value == 100
            assert frame.props[1].flow_mol_phase_comp[p, i].lb == 0
            assert frame.props[1].flow_mol_phase_comp[p, i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

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

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in frame.params._phase_component_set
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_expressions(self, frame):
        assert isinstance(frame.props[1].flow_mol, Expression)
        assert len(frame.props[1].flow_mol) == 1
        assert str(frame.props[1].flow_mol.expr) == str(
            sum(
                frame.props[1].flow_mol_phase_comp[i]
                for i in frame.props[1].params._phase_component_set
            )
        )

        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 2
        for p in frame.props[1].flow_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].flow_mol_phase[p].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for j in frame.props[1].params.component_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].flow_mol_comp, Expression)
        assert len(frame.props[1].flow_mol_comp) == 3
        for j in frame.props[1].flow_mol_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].flow_mol_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].mole_frac_comp, Expression)
        assert len(frame.props[1].mole_frac_comp) == 3
        for j in frame.props[1].mole_frac_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].mole_frac_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
                / frame.props[1].flow_mol
            )

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 2
        for p in frame.props[1].phase_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].phase_frac[p].expr) == str(
                frame.props[1].flow_mol_phase[p] / frame.props[1].flow_mol
            )

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 6
        for (p, i) in frame.props[1].mole_frac_phase_comp_eq:
            assert str(frame.props[1].mole_frac_phase_comp_eq[p, i].body) == str(
                frame.props[1].mole_frac_phase_comp[p, i]
                * frame.props[1].flow_mol_phase[p]
                - frame.props[1].flow_mol_phase_comp[p, i]
            )

        assert_units_consistent(frame.props[1])


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

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert not frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol_phase_comp, Var)
        assert len(frame.props[1].flow_mol_phase_comp) == 9
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert (p, i) in frame.props[1].params._phase_component_set
            assert frame.props[1].flow_mol_phase_comp[p, i].value is None
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

        assert isinstance(frame.props[1].pressure, Var)
        assert frame.props[1].pressure.value is None
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert isinstance(frame.props[1].temperature, Var)
        assert frame.props[1].temperature.value is None
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 9
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in frame.params._phase_component_set
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_expressions(self, frame):
        assert isinstance(frame.props[1].flow_mol, Expression)
        assert len(frame.props[1].flow_mol) == 1
        assert str(frame.props[1].flow_mol.expr) == str(
            sum(
                frame.props[1].flow_mol_phase_comp[i]
                for i in frame.props[1].params._phase_component_set
            )
        )

        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 3
        for p in frame.props[1].flow_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].flow_mol_phase[p].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for j in frame.props[1].params.component_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].flow_mol_comp, Expression)
        assert len(frame.props[1].flow_mol_comp) == 3
        for j in frame.props[1].flow_mol_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].flow_mol_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].mole_frac_comp, Expression)
        assert len(frame.props[1].mole_frac_comp) == 3
        for j in frame.props[1].mole_frac_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].mole_frac_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
                / frame.props[1].flow_mol
            )

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 3
        for p in frame.props[1].phase_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].phase_frac[p].expr) == str(
                frame.props[1].flow_mol_phase[p] / frame.props[1].flow_mol
            )

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 9
        for (p, i) in frame.props[1].mole_frac_phase_comp_eq:
            assert str(frame.props[1].mole_frac_phase_comp_eq[p, i].body) == str(
                frame.props[1].mole_frac_phase_comp[p, i]
                * frame.props[1].flow_mol_phase[p]
                - frame.props[1].flow_mol_phase_comp[p, i]
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
                "flow_mol_phase_comp": (0, 100, 200),
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

        # Build state block
        m.props = m.params.build_state_block([1], defined_state=True)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_always_flash(self, frame):
        define_state(frame.props[1])

        assert not frame.props[1].always_flash

    @pytest.mark.unit
    def test_vars(self, frame):
        # Check that all necessary variables have been constructed and have
        # the correct values
        assert isinstance(frame.props[1].flow_mol_phase_comp, Var)
        assert len(frame.props[1].flow_mol_phase_comp) == 9
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert (p, i) in frame.props[1].params._phase_component_set
            assert frame.props[1].flow_mol_phase_comp[p, i].value == 100
            assert frame.props[1].flow_mol_phase_comp[p, i].lb == 0
            assert frame.props[1].flow_mol_phase_comp[p, i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

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

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 9
        for i in frame.props[1].mole_frac_phase_comp:
            assert i in frame.params._phase_component_set
            assert frame.props[1].mole_frac_phase_comp[i].value == 1 / 3
        assert check_units_equivalent(frame.props[1].mole_frac_phase_comp, None)

    @pytest.mark.unit
    def test_expressions(self, frame):
        assert isinstance(frame.props[1].flow_mol, Expression)
        assert len(frame.props[1].flow_mol) == 1
        assert str(frame.props[1].flow_mol.expr) == str(
            sum(
                frame.props[1].flow_mol_phase_comp[i]
                for i in frame.props[1].params._phase_component_set
            )
        )

        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 3
        for p in frame.props[1].flow_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].flow_mol_phase[p].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for j in frame.props[1].params.component_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].flow_mol_comp, Expression)
        assert len(frame.props[1].flow_mol_comp) == 3
        for j in frame.props[1].flow_mol_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].flow_mol_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
            )

        assert isinstance(frame.props[1].mole_frac_comp, Expression)
        assert len(frame.props[1].mole_frac_comp) == 3
        for j in frame.props[1].mole_frac_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].mole_frac_comp[j].expr) == str(
                sum(
                    frame.props[1].flow_mol_phase_comp[p, j]
                    for p in frame.props[1].params.phase_list
                    if (p, j) in frame.props[1].params._phase_component_set
                )
                / frame.props[1].flow_mol
            )

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 3
        for p in frame.props[1].phase_frac:
            assert p in frame.params.phase_list
            assert str(frame.props[1].phase_frac[p].expr) == str(
                frame.props[1].flow_mol_phase[p] / frame.props[1].flow_mol
            )

    @pytest.mark.unit
    def test_constraints(self, frame):
        # Check that the correct constraints are present
        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 9
        for (p, i) in frame.props[1].mole_frac_phase_comp_eq:
            assert str(frame.props[1].mole_frac_phase_comp_eq[p, i].body) == str(
                frame.props[1].mole_frac_phase_comp[p, i]
                * frame.props[1].flow_mol_phase[p]
                - frame.props[1].flow_mol_phase_comp[p, i]
            )

        assert_units_consistent(frame.props[1])


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
            state_definition=FpcTP,
            pressure_ref=100000.0,
            temperature_ref=300,
            state_bounds={
                "flow_mol_phase_comp": (0, 0.1, 0.2, pyunits.kmol / pyunits.s),
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
        m.props = m.params.build_state_block([1], defined_state=False)

        # Add necessary variables that would be built by other methods
        m.props[1].dens_mol_phase = Var(m.params.phase_list, initialize=1)
        m.props[1].enth_mol_phase = Var(m.params.phase_list, initialize=1)

        return m

    @pytest.mark.unit
    def test_convert_vars(self, frame):
        # Check that all state var values and bounds were converted correctly
        for (p, i) in frame.props[1].flow_mol_phase_comp:
            assert frame.props[1].flow_mol_phase_comp[p, i].value == 100
            assert frame.props[1].flow_mol_phase_comp[p, i].lb == 0
            assert frame.props[1].flow_mol_phase_comp[p, i].ub == 200
        assert check_units_equivalent(
            frame.props[1].flow_mol_phase_comp, pyunits.mol / pyunits.s
        )

        assert frame.props[1].pressure.value == 3e5
        assert frame.props[1].pressure.lb == 1e5
        assert frame.props[1].pressure.ub == 5e5
        assert check_units_equivalent(frame.props[1].pressure, pyunits.Pa)

        assert frame.props[1].temperature.value == 345
        assert frame.props[1].temperature.lb == 290
        assert frame.props[1].temperature.ub == 400
        assert check_units_equivalent(frame.props[1].temperature, pyunits.K)

        # Check supporting variables
        assert isinstance(frame.props[1].flow_mol_phase, Expression)
        assert len(frame.props[1].flow_mol_phase) == 2

        assert isinstance(frame.props[1].mole_frac_phase_comp, Var)
        assert len(frame.props[1].mole_frac_phase_comp) == 6

        assert isinstance(frame.props[1].phase_frac, Expression)
        assert len(frame.props[1].phase_frac) == 2

        assert isinstance(frame.props[1].mole_frac_phase_comp_eq, Constraint)
        assert len(frame.props[1].mole_frac_phase_comp_eq) == 6

    @pytest.mark.unit
    def test_calculate_scaling_factors(self, frame):
        frame.props[1].calculate_scaling_factors()

        assert len(frame.props[1].scaling_factor) == 25
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_phase["a"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_phase["b"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_comp["c1"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_comp["c2"]] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].flow_mol_comp["c3"]] == 1e-2
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
            "flow_mol_phase_comp": frame.props[1].flow_mol_phase_comp,
            "temperature": frame.props[1].temperature,
            "pressure": frame.props[1].pressure,
        }

    @pytest.mark.unit
    def test_define_display_vars(self, frame):
        assert frame.props[1].define_display_vars() == {
            "Molar Flowrate": frame.props[1].flow_mol_phase_comp,
            "Temperature": frame.props[1].temperature,
            "Pressure": frame.props[1].pressure,
        }


# -----------------------------------------------------------------------------
# Test example from Austin Ladshaw via WaterTAP project
# This example revealed some issues with initialization of FpcTP states with
# phase equilibrium
thermo_config_no_rxn = {
    "components": {
        "H2O": {
            "type": Solvent,
            # Define the methods used to calculate the following properties
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": NIST,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            # Parameter data is always associated with the methods defined above
            "parameter_data": {
                "mw": (18.0153, pyunits.g / pyunits.mol),
                "pressure_crit": (220.64e5, pyunits.Pa),
                "temperature_crit": (647, pyunits.K),
                # Comes from Perry's Handbook:  p. 2-98
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (5.459, pyunits.kmol * pyunits.m**-3),
                    "2": (0.30542, pyunits.dimensionless),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, pyunits.dimensionless),
                },
                "enth_mol_form_liq_comp_ref": (-285.830, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (0, pyunits.kJ / pyunits.mol),
                # Comes from Perry's Handbook:  p. 2-174
                "cp_mol_liq_comp_coeff": {
                    "1": (2.7637e5, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (30.09200, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        6.832514,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-1,
                    ),
                    "C": (
                        6.793435,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-2,
                    ),
                    "D": (
                        -2.534480,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-3,
                    ),
                    "E": (
                        0.082139,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**2,
                    ),
                    "F": (-250.8810, pyunits.kJ / pyunits.mol),
                    "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "entr_mol_form_liq_comp_ref": (
                    69.95,
                    pyunits.J / pyunits.K / pyunits.mol,
                ),
                "pressure_sat_comp_coeff": {
                    "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                    "B": (1435.264, pyunits.K),
                    "C": (-64.848, pyunits.K),
                },
            },
        },
        "CO2": {
            "type": Solute,
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "cp_mol_liq_comp": Perrys,
            "entr_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": NIST,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (44.0095, pyunits.g / pyunits.mol),
                "pressure_crit": (73.825e5, pyunits.Pa),
                "temperature_crit": (304.23, pyunits.K),
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (0.000789, pyunits.kmol * pyunits.m**-3),
                    "2": (0.000956, pyunits.dimensionless),
                    "3": (500.78, pyunits.K),
                    "4": (0.94599, pyunits.dimensionless),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (24.99735, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        55.18696,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-1,
                    ),
                    "C": (
                        -33.69137,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-2,
                    ),
                    "D": (
                        7.948387,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-3,
                    ),
                    "E": (
                        -0.136638,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**2,
                    ),
                    "F": (-403.6075, pyunits.kJ / pyunits.mol),
                    "G": (228.2431, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (-8.3043e6, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (1.0437e5, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (4.333e2, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (6.0052e-1, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (-285.83, pyunits.kJ / pyunits.mol),
                "enth_mol_form_vap_comp_ref": (-393.52, pyunits.kJ / pyunits.mol),
                "entr_mol_form_liq_comp_ref": (0, pyunits.J / pyunits.K / pyunits.mol),
                "entr_mol_form_vap_comp_ref": (213.6, pyunits.J / pyunits.mol),
                "pressure_sat_comp_coeff": {
                    "A": (6.81228, None),
                    "B": (1301.679, pyunits.K),
                    "C": (-3.494, pyunits.K),
                },
            },
        },
    },
    "phases": {
        "Liq": {"type": AqueousPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
    },
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FpcTP,
    "state_bounds": {
        "temperature": (273.15, 300, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}


@pytest.mark.component
def test_phase_equilibrium_initialization():
    # Create a pyomo model object
    model = ConcreteModel()
    model.fs = FlowsheetBlock(dynamic=False)

    model.fs.thermo_params = GenericParameterBlock(**thermo_config_no_rxn)

    model.fs.state = model.fs.thermo_params.build_state_block(
        model.fs.time, defined_state=False
    )

    model.fs.state[0].pressure.set_value(101325.0)
    model.fs.state[0].temperature.set_value(298.0)

    model.fs.state[0].flow_mol_phase_comp["Vap", "CO2"].set_value(0.0005 * 10)
    model.fs.state[0].flow_mol_phase_comp["Liq", "CO2"].set_value(1e-8)
    model.fs.state[0].flow_mol_phase_comp["Vap", "H2O"].set_value(1e-8)
    model.fs.state[0].flow_mol_phase_comp["Liq", "H2O"].set_value((1 - 0.0005) * 10)

    assert_units_consistent(model)
    # We expect 6 state variables, but two additional constraints for phase equilibrium
    assert degrees_of_freedom(model) == 6 - 2

    model.fs.state.initialize()

    # Check that degrees of freedom are still the same
    assert degrees_of_freedom(model) == 6 - 2

    # As the phase equilibrium constraints were not solved, we expect these to have a large residual
    large_res = large_residuals_set(model.fs.state[0])
    assert len(large_res) == 2
    for i in large_res:
        assert i.name in [
            "fs.state[0.0].equilibrium_constraint[Vap,Liq,H2O]",
            "fs.state[0.0].equilibrium_constraint[Vap,Liq,CO2]",
        ]
