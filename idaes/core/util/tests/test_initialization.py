##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for math util methods.
"""

import pytest
from pyomo.environ import (Block, ConcreteModel,  Constraint, Expression,
                           Set, SolverFactory, Var, value, 
                           TransformationFactory, TerminationCondition)
from pyomo.network import Arc, Port

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
        AqueousEnzymeParameterBlock, EnzymeReactionParameterBlock,
        EnzymeReactionBlock)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.unit_models.cstr import CSTR
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            propagate_state,
                                            solve_indexed_blocks,
                                            initialize_by_time_element)

__author__ = "Andrew Lee"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = m.fs.pp.state_block_class(default={'parameters': m.fs.pp})

    for i in m.fs.sb.flow_mol_phase_comp:
        assert not m.fs.sb.flow_mol_phase_comp[i].fixed
        assert m.fs.sb.flow_mol_phase_comp[i].value == 2
    assert not m.fs.sb.pressure.fixed
    assert m.fs.sb.pressure.value == 1e5
    assert not m.fs.sb.temperature.fixed
    assert m.fs.sb.temperature.value == 300

    return m


def test_fix_state_vars_basic(model):
    flags = fix_state_vars(model.fs.sb)

    for i in model.fs.sb.flow_mol_phase_comp:
        assert model.fs.sb.flow_mol_phase_comp[i].fixed
        assert model.fs.sb.flow_mol_phase_comp[i].value == 2
    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1e5
    assert model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 300

    assert not flags[None, "component_flow_phase", ("p1", "c1")]
    assert not flags[None, "component_flow_phase", ("p1", "c2")]
    assert not flags[None, "component_flow_phase", ("p2", "c1")]
    assert not flags[None, "component_flow_phase", ("p2", "c2")]
    assert not flags[None, "pressure", None]
    assert not flags[None, "temperature", None]


def test_fix_state_vars_None_value(model):
    model.fs.sb.pressure.value = None

    with pytest.raises(ConfigurationError):
        fix_state_vars(model.fs.sb)


def test_fix_state_vars_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {"component_flow_phase": {("p1", "c1"): 1,
                                     ("p1", "c2"): 2,
                                     ("p2", "c1"): 3,
                                     ("p2", "c2"): 4},
                  "pressure": 2e5,
                  "temperature": 500}

    flags = fix_state_vars(model.fs.sb, state_args)

    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 1
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 3
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 4

    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 2e5
    assert model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 500

    assert not flags[None, "component_flow_phase", ("p1", "c1")]
    assert not flags[None, "component_flow_phase", ("p1", "c2")]
    assert not flags[None, "component_flow_phase", ("p2", "c1")]
    assert not flags[None, "component_flow_phase", ("p2", "c2")]
    assert not flags[None, "pressure", None]
    assert not flags[None, "temperature", None]


def test_fix_state_vars_partial_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {"component_flow_phase": {("p1", "c1"): 1,
                                     ("p1", "c2"): 2,
                                     ("p2", "c1"): 3,
                                     ("p2", "c2"): 4}}

    flags = fix_state_vars(model.fs.sb, state_args)

    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 1
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 3
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 4

    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1e5
    assert model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 300

    assert not flags[None, "component_flow_phase", ("p1", "c1")]
    assert not flags[None, "component_flow_phase", ("p1", "c2")]
    assert not flags[None, "component_flow_phase", ("p2", "c1")]
    assert not flags[None, "component_flow_phase", ("p2", "c2")]
    assert not flags[None, "pressure", None]
    assert not flags[None, "temperature", None]


def test_fix_state_vars_guesses_mismatch_index(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {"component_flow_phase": {("p1", "c1"): 1,
                                     ("p1", "c2"): 2,
                                     ("p2", "c1"): 3},
                  "pressure": 2e5,
                  "temperature": 500}

    with pytest.raises(ConfigurationError):
        fix_state_vars(model.fs.sb, state_args)


def test_fix_state_vars_fixed_no_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    model.fs.sb.flow_mol_phase_comp["p1", "c1"].fix(10)
    model.fs.sb.pressure.fix(1.5e5)

    flags = fix_state_vars(model.fs.sb)

    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 10
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 2
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 2

    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1.5e5
    assert model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 300

    # Pressure and component_flow_phase[p1, c1] should be True
    assert flags[None, "component_flow_phase", ("p1", "c1")]
    assert not flags[None, "component_flow_phase", ("p1", "c2")]
    assert not flags[None, "component_flow_phase", ("p2", "c1")]
    assert not flags[None, "component_flow_phase", ("p2", "c2")]
    assert flags[None, "pressure", None]
    assert not flags[None, "temperature", None]


def test_fix_state_vars_fixed_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    model.fs.sb.flow_mol_phase_comp["p1", "c1"].fix(10)
    model.fs.sb.pressure.fix(1.5e5)

    state_args = {"component_flow_phase": {("p1", "c1"): 1,
                                     ("p1", "c2"): 2,
                                     ("p2", "c1"): 3,
                                     ("p2", "c2"): 4},
                  "pressure": 2e5,
                  "temperature": 500}

    flags = fix_state_vars(model.fs.sb, state_args)

    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 10
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 3
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 4

    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1.5e5
    assert model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 500

    # Pressure and component_flow_phase[p1, c1] should be True
    assert flags[None, "component_flow_phase", ("p1", "c1")]
    assert not flags[None, "component_flow_phase", ("p1", "c2")]
    assert not flags[None, "component_flow_phase", ("p2", "c1")]
    assert not flags[None, "component_flow_phase", ("p2", "c2")]
    assert flags[None, "pressure", None]
    assert not flags[None, "temperature", None]


def test_revert_state_vars_basic(model):
    flags = fix_state_vars(model.fs.sb)

    revert_state_vars(model.fs.sb, flags)

    for i in model.fs.sb.flow_mol_phase_comp:
        assert not model.fs.sb.flow_mol_phase_comp[i].fixed
        assert model.fs.sb.flow_mol_phase_comp[i].value == 2
    assert not model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1e5
    assert not model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 300


def test_revert_state_vars_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {"component_flow_phase": {("p1", "c1"): 1,
                                     ("p1", "c2"): 2,
                                     ("p2", "c1"): 3,
                                     ("p2", "c2"): 4},
                  "pressure": 2e5,
                  "temperature": 500}

    flags = fix_state_vars(model.fs.sb, state_args)

    revert_state_vars(model.fs.sb, flags)

    assert not model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 1
    assert not model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert not model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 3
    assert not model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 4

    assert not model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 2e5
    assert not model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 500


def test_revert_state_vars_fixed_no_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    model.fs.sb.flow_mol_phase_comp["p1", "c1"].fix(10)
    model.fs.sb.pressure.fix(1.5e5)

    flags = fix_state_vars(model.fs.sb)
    revert_state_vars(model.fs.sb, flags)

    # Pressure and componet_flow[p1, c1] should still be fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 10
    assert not model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert not model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 2
    assert not model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 2

    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1.5e5
    assert not model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 300


def test_revert_state_vars_fixed_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    model.fs.sb.flow_mol_phase_comp["p1", "c1"].fix(10)
    model.fs.sb.pressure.fix(1.5e5)

    state_args = {"component_flow_phase": {("p1", "c1"): 1,
                                     ("p1", "c2"): 2,
                                     ("p2", "c1"): 3,
                                     ("p2", "c2"): 4},
                  "pressure": 2e5,
                  "temperature": 500}

    flags = fix_state_vars(model.fs.sb, state_args)
    revert_state_vars(model.fs.sb, flags)

    # Pressure and componet_flow[p1, c1] should still be fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c1")].value == 10
    assert not model.fs.sb.flow_mol_phase_comp[("p1", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p1", "c2")].value == 2
    assert not model.fs.sb.flow_mol_phase_comp[("p2", "c1")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c1")].value == 3
    assert not model.fs.sb.flow_mol_phase_comp[("p2", "c2")].fixed
    assert model.fs.sb.flow_mol_phase_comp[("p2", "c2")].value == 4

    assert model.fs.sb.pressure.fixed
    assert model.fs.sb.pressure.value == 1.5e5
    assert not model.fs.sb.temperature.fixed
    assert model.fs.sb.temperature.value == 500


def test_revert_state_vars_flag_mismatch(model):
    flags = {(None, 'component_flow_phase', ('p1', 'c1')): False,
             (None, 'component_flow_phase', ('p1', 'c2')): False,
             (None, 'component_flow_phase', ('p2', 'c1')): False,
             (None, 'pressure', None): False}

    with pytest.raises(ConfigurationError):
        revert_state_vars(model.fs.sb, flags)


def test_propagate_state():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    # Set values on first block
    m.b1.v1.value = 10
    m.b1.v2[1].value = 20
    m.b1.v2[2].value = 30

    # Make sure vars in block 2 haven't been changed accidentally
    assert m.b2.v1.value is None
    assert m.b2.v2[1].value is None
    assert m.b2.v2[2].value is None

    propagate_state(m.s1)

    # Check that values were propagated correctly
    assert m.b2.v1.value == m.b1.v1.value
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value

    assert m.b1.v1.fixed is False
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v1.fixed is False
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


def test_propagate_state_reverse():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    # Test reverse propogation - set values on second block
    m.b2.v1.value = 100
    m.b2.v2[1].value = 200
    m.b2.v2[2].value = 300

    # Make sure vars in block 1 haven't been changed accidentally
    assert m.b1.v1.value is None
    assert m.b1.v2[1].value is None
    assert m.b1.v2[2].value is None

    propagate_state(m.s1, direction="backward")

    # Check that values were propagated correctly
    assert m.b2.v1.value == m.b1.v1.value
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value

    assert m.b1.v1.fixed is False
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v1.fixed is False
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


def test_propagate_state_fixed():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    # Set values on first block
    m.b1.v1.value = 10
    m.b1.v2[1].value = 20
    m.b1.v2[2].value = 30

    # Make sure vars in block 2 haven't been changed accidentally
    assert m.b2.v1.value is None
    assert m.b2.v2[1].value is None
    assert m.b2.v2[2].value is None

    # Fix v1 in block 2
    m.b2.v1.fix(500)

    propagate_state(m.s1)

    # Check that values were propagated correctly
    assert m.b2.v1.value == 500
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value

    assert m.b1.v1.fixed is False
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v1.fixed is True
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


def test_propagate_state_Expression():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.e = Expression(expr=b.v1)

        b.p = Port()
        b.p.add(b.e, "E")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    with pytest.raises(TypeError):
        propagate_state(m.s1)


def test_propagate_state_invalid_stream():
    m = ConcreteModel()

    with pytest.raises(TypeError):
        propagate_state(m)


def test_propagate_state_invalid_direction():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    with pytest.raises(ValueError):
        propagate_state(m.s1, direction="foo")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve_indexed_block_list():
    # Create an indexed block and try to solve it
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3])

    def block_rule(b, x):
        b.v = Var(initialize=1.0)
        b.c = Constraint(expr=b.v == 2.0)
    m.b = Block(m.s, rule=block_rule)

    solve_indexed_blocks(solver=solver, blocks=[m.b])

    for i in m.s:
        assert value(m.b[i].v == 2.0)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve_indexed_block_IndexedBlock():
    # Create an indexed block and try to solve it
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3])

    def block_rule(b, x):
        b.v = Var(initialize=1.0)
        b.c = Constraint(expr=b.v == 2.0)
    m.b = Block(m.s, rule=block_rule)

    solve_indexed_blocks(solver=solver, blocks=m.b)

    for i in m.s:
        assert value(m.b[i].v == 2.0)


def test_solve_indexed_block_error():
    # Try solve_indexed_block on non-block object
    with pytest.raises(TypeError):
        solve_indexed_blocks(solver=None, blocks=[1, 2, 3])


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_by_time_element():
    horizon = 6
    time_set = [0, horizon]
    ntfe = 60 # For a finite element every six seconds
    ntcp = 2
    m = ConcreteModel(name='CSTR model for testing')
    m.fs = FlowsheetBlock(default={'dynamic': True,
                                   'time_set': time_set})

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(
            default={'property_package': m.fs.properties})
    m.fs.cstr = CSTR(default={"property_package": m.fs.properties,
                              "reaction_package": m.fs.reactions,
                              "material_balance_type": MaterialBalanceType.componentTotal,
                              "energy_balance_type": EnergyBalanceType.none,
                              "momentum_balance_type": MomentumBalanceType.none})

    # Time discretization
    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')

    # Fix geometry variables
    m.fs.cstr.volume.fix(1.0)

    # Fix initial conditions:
    for p, j in m.fs.properties.phase_list*m.fs.properties.component_list:
        m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0)

    # Fix inlet conditions
    for t, j in m.fs.time*m.fs.properties.component_list:
        if t <= 2:
            if j == 'E':
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(11.91*0.1)
            elif j == 'S':
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(12.92*2.1)
            else:
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(0)
        elif t <= 4:
            if j == 'E':
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(5.95*0.1)
            elif j == 'S':
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(12.92*2.1)
            else:
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(0)
        else:
            if j == 'E':
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(8.95*0.1)
            elif j == 'S':
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(16.75*2.1)
            else:
                m.fs.cstr.inlet.flow_mol_comp[t, j].fix(0)

    m.fs.cstr.inlet.flow_rate.fix(2.2)
    m.fs.cstr.inlet.temperature.fix(300)

    # Fix outlet conditions
    m.fs.cstr.outlet.flow_rate.fix(2.2)
    m.fs.cstr.outlet.temperature.fix(300)

    assert degrees_of_freedom(m) == 0

    initialize_by_time_element(m.fs, m.fs.time, solver=solver)

    assert degrees_of_freedom(m) == 0

    # Assert that the result looks how we expect
    assert m.fs.cstr.outlet.conc_mol[0, 'S'].value == 0
    assert abs(m.fs.cstr.outlet.conc_mol[2, 'S'].value - 11.076) < 1e-2
    assert abs(m.fs.cstr.outlet.conc_mol[4, 'P'].value - 0.3577) < 1e-3
    assert abs(m.fs.cstr.outlet.conc_mol[6, 'E'].value - 0.0286) < 1e-3

    # Assert that model is still fixed and deactivated as expected
    assert (
    m.fs.cstr.control_volume.material_holdup[m.fs.time.first(), 'aq', 'S'].fixed)
    for t in m.fs.time:
        if t != m.fs.time.first():
            assert (not 
    m.fs.cstr.control_volume.material_holdup[t, 'aq', 'S'].fixed)
        assert (
    m.fs.cstr.control_volume.material_holdup_calculation[t, 'aq', 'C'].active)
        assert m.fs.cstr.control_volume.properties_out[t].active
        assert m.fs.cstr.outlet.temperature[t].fixed
        assert not m.fs.cstr.outlet.flow_mol_comp[t, 'S'].fixed
        assert m.fs.cstr.inlet.flow_mol_comp[t, 'S'].fixed

    # Assert that constraints are feasible after initialization
    for con in m.fs.component_data_objects(Constraint, active=True):
        assert value(con.body) - value(con.upper) < 1e-5
        assert value(con.lower) - value(con.body) < 1e-5

    results = solver.solve(m.fs)
    assert results.solver.termination_condition == TerminationCondition.optimal

if __name__ == '__main__':
    test_initialize_by_time_element()
