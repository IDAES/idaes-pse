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
Tests for math util methods.
"""

import pytest
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    exp,
    Set,
    Var,
    value,
    Param,
    Reals,
    units as pyunits,
    TransformationFactory,
    check_optimal_termination,
)
from pyomo.network import Arc, Port

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    ReactionParameterBlock,
    ReactionBlockBase,
    ReactionBlockDataBase,
    MaterialFlowBasis,
    Component,
    Phase,
)
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.unit_models import CSTR
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    propagate_state,
    solve_indexed_blocks,
    initialize_by_time_element,
)
from idaes.core.solvers import get_solver

__author__ = "Andrew Lee"


# Set up solver
solver = get_solver()


@declare_process_block_class("AqueousEnzymeParameterBlock")
class ParameterData(PhysicalParameterBlock):
    """
    Parameter block for the aqueous enzyme reaction in the biochemical CSTR
    used by Christofides and Daoutidis, 1996, presented by Heineken et al, 1967
    """

    def build(self):
        super(ParameterData, self).build()

        # all components are in the aqueous phase
        self.aq = Phase()
        self.S = Component()
        self.E = Component()
        self.C = Component()
        self.P = Component()
        self.Solvent = Component()

        self._state_block_class = AqueousEnzymeStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.minute,
                "length": pyunits.m,
                "amount": pyunits.kmol,
                "temperature": pyunits.K,
                "mass": pyunits.kg,
            }
        )


class _AqueousEnzymeStateBlock(StateBlock):
    def initialize(blk):
        pass


@declare_process_block_class(
    "AqueousEnzymeStateBlock", block_class=_AqueousEnzymeStateBlock
)
class AqueousEnzymeStateBlockData(StateBlockData):
    def build(self):
        super(AqueousEnzymeStateBlockData, self).build()

        self.conc_mol = Var(
            self._params.component_list,
            domain=Reals,
            doc="Component molar concentration [kmol/m^3]",
        )

        self.flow_mol_comp = Var(
            self._params.component_list,
            domain=Reals,
            doc="Molar component flow rate [kmol/min]",
        )

        self.flow_vol = Var(
            domain=Reals, doc="Volumetric flow rate out of reactor [m^3/min]"
        )

        self.temperature = Var(
            initialize=303, domain=Reals, doc="Temperature within reactor [K]"
        )

        if not self.config.defined_state:
            self.conc_mol["Solvent"].fix(1.0)

        def flow_mol_comp_rule(b, j):
            return b.flow_mol_comp[j] == b.flow_vol * b.conc_mol[j]

        self.flow_mol_comp_eqn = Constraint(
            self._params.component_list,
            rule=flow_mol_comp_rule,
            doc="Outlet component molar flow rate equation",
        )

    def get_material_density_terms(b, p, j):
        return b.conc_mol[j]

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_comp[j]

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def get_enthalpy_flow_terms(b, p):
        return b.flow_vol * b.temperature

    def get_energy_density_terms(b, p):
        return b.temperature

    def define_state_vars(b):
        return {
            "conc_mol": b.conc_mol,
            "temperature": b.temperature,
            "flow_vol": b.flow_vol,
        }


@declare_process_block_class("EnzymeReactionParameterBlock")
class EnzymeReactionParameterData(ReactionParameterBlock):
    """
    Enzyme reaction:
    S + E <-> C -> P + E
    """

    def build(self):
        super(EnzymeReactionParameterData, self).build()

        self._reaction_block_class = EnzymeReactionBlock

        self.rate_reaction_idx = Set(initialize=["R1", "R2", "R3"])
        self.rate_reaction_stoichiometry = {
            ("R1", "aq", "S"): -1,
            ("R1", "aq", "E"): -1,
            ("R1", "aq", "C"): 1,
            ("R1", "aq", "P"): 0,
            ("R1", "aq", "Solvent"): 0,
            ("R2", "aq", "S"): 1,
            ("R2", "aq", "E"): 1,
            ("R2", "aq", "C"): -1,
            ("R2", "aq", "P"): 0,
            ("R2", "aq", "Solvent"): 0,
            ("R3", "aq", "S"): 0,
            ("R3", "aq", "E"): 1,
            ("R3", "aq", "C"): -1,
            ("R3", "aq", "P"): 1,
            ("R3", "aq", "Solvent"): 0,
        }

        self.act_energy = Param(
            self.rate_reaction_idx,
            initialize={"R1": 8.0e3, "R2": 9.0e3, "R3": 1.0e4},
            doc="Activation energy [kcal/kmol]",
        )

        self.gas_const = Param(initialize=1.987, doc="Gas constant R [kcal/kmol/K]")

        self.temperature_ref = Param(initialize=300.0, doc="Reference temperature")

        self.k_rxn = Param(
            self.rate_reaction_idx,
            initialize={"R1": 3.36e6, "R2": 1.80e6, "R3": 5.79e7},
            doc="Pre-exponential rate constant in Arrhenius expression",
        )

    #    self.reaction_block_class = EnzymeReactionBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": pyunits.minute,
                "length": pyunits.m,
                "amount": pyunits.kmol,
                "temperature": pyunits.K,
                "mass": pyunits.kg,
            }
        )


class _EnzymeReactionBlock(ReactionBlockBase):
    def initialize(blk):
        # initialize for reaction rates for each data object
        pass


@declare_process_block_class("EnzymeReactionBlock", block_class=_EnzymeReactionBlock)
class EnzymeReactionBlockData(ReactionBlockDataBase):
    def build(self):
        super(EnzymeReactionBlockData, self).build()

        self.reaction_coef = Var(
            self._params.rate_reaction_idx,
            domain=Reals,
            doc="Reaction rate coefficient",
        )

        self.reaction_rate = Var(
            self._params.rate_reaction_idx,
            domain=Reals,
            doc="Reaction rate [kmol/m^3/min]",
        )

        self.dh_rxn = Param(
            self._params.rate_reaction_idx,
            domain=Reals,
            doc="Heat of reaction",
            initialize={
                "R1": 1e3 / 900 / 0.231,
                "R2": 1e3 / 900 / 0.231,
                "R3": 5e3 / 900 / 0.231,
            },
        )

        def reaction_rate_rule(b, r):
            if r == "R1":
                return (
                    b.reaction_rate[r]
                    == b.reaction_coef[r]
                    * b.state_ref.conc_mol["S"]
                    * b.state_ref.conc_mol["E"]
                )
            elif r == "R2":
                return (
                    b.reaction_rate[r] == b.reaction_coef[r] * b.state_ref.conc_mol["C"]
                )
            elif r == "R3":
                return (
                    b.reaction_rate[r] == b.reaction_coef[r] * b.state_ref.conc_mol["C"]
                )

        self.reaction_rate_eqn = Constraint(
            self._params.rate_reaction_idx, rule=reaction_rate_rule
        )

        def arrhenius_rule(b, r):
            return b.reaction_coef[r] == b._params.k_rxn[r] * exp(
                -b._params.act_energy[r] / b._params.gas_const / b.state_ref.temperature
            )

        self.arrhenius_eqn = Constraint(
            self._params.rate_reaction_idx, rule=arrhenius_rule
        )

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar


@pytest.fixture
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.sb = m.fs.pp.state_block_class(parameters=m.fs.pp)

    for i in m.fs.sb.flow_mol_phase_comp:
        assert not m.fs.sb.flow_mol_phase_comp[i].fixed
        assert m.fs.sb.flow_mol_phase_comp[i].value == 2
    assert not m.fs.sb.pressure.fixed
    assert m.fs.sb.pressure.value == 1e5
    assert not m.fs.sb.temperature.fixed
    assert m.fs.sb.temperature.value == 300

    return m


@pytest.mark.unit
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


@pytest.mark.unit
def test_fix_state_vars_None_value(model):
    model.fs.sb.pressure.value = None

    with pytest.raises(ConfigurationError):
        fix_state_vars(model.fs.sb)


@pytest.mark.unit
def test_fix_state_vars_guesses(model):
    # Note that flow_mol_phase_comp is labled as component_flow
    # in define_state_vars
    state_args = {
        "component_flow_phase": {
            ("p1", "c1"): 1,
            ("p1", "c2"): 2,
            ("p2", "c1"): 3,
            ("p2", "c2"): 4,
        },
        "pressure": 2e5,
        "temperature": 500,
    }

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


@pytest.mark.unit
def test_fix_state_vars_partial_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {
        "component_flow_phase": {
            ("p1", "c1"): 1,
            ("p1", "c2"): 2,
            ("p2", "c1"): 3,
            ("p2", "c2"): 4,
        }
    }

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


@pytest.mark.unit
def test_fix_state_vars_guesses_mismatch_index(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {
        "component_flow_phase": {("p1", "c1"): 1, ("p1", "c2"): 2, ("p2", "c1"): 3},
        "pressure": 2e5,
        "temperature": 500,
    }

    with pytest.raises(ConfigurationError):
        fix_state_vars(model.fs.sb, state_args)


@pytest.mark.unit
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


@pytest.mark.unit
def test_fix_state_vars_fixed_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    model.fs.sb.flow_mol_phase_comp["p1", "c1"].fix(10)
    model.fs.sb.pressure.fix(1.5e5)

    state_args = {
        "component_flow_phase": {
            ("p1", "c1"): 1,
            ("p1", "c2"): 2,
            ("p2", "c1"): 3,
            ("p2", "c2"): 4,
        },
        "pressure": 2e5,
        "temperature": 500,
    }

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


@pytest.mark.unit
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


@pytest.mark.unit
def test_revert_state_vars_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    state_args = {
        "component_flow_phase": {
            ("p1", "c1"): 1,
            ("p1", "c2"): 2,
            ("p2", "c1"): 3,
            ("p2", "c2"): 4,
        },
        "pressure": 2e5,
        "temperature": 500,
    }

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


@pytest.mark.unit
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


@pytest.mark.unit
def test_revert_state_vars_fixed_guesses(model):
    # Note that flow_mol_phase_comp is labled as compoennt_flow
    # in define_state_vars
    model.fs.sb.flow_mol_phase_comp["p1", "c1"].fix(10)
    model.fs.sb.pressure.fix(1.5e5)

    state_args = {
        "component_flow_phase": {
            ("p1", "c1"): 1,
            ("p1", "c2"): 2,
            ("p2", "c1"): 3,
            ("p2", "c2"): 4,
        },
        "pressure": 2e5,
        "temperature": 500,
    }

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


@pytest.mark.unit
def test_revert_state_vars_flag_mismatch(model):
    flags = {
        (None, "component_flow_phase", ("p1", "c1")): False,
        (None, "component_flow_phase", ("p1", "c2")): False,
        (None, "component_flow_phase", ("p2", "c1")): False,
        (None, "pressure", None): False,
    }

    with pytest.raises(ConfigurationError):
        revert_state_vars(model.fs.sb, flags)


@pytest.mark.unit
def test_propagate_state():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)
        b.e1 = Expression(expr=b.v1)

        @b.Expression(b.s)
        def e2(blk, i):
            return b.v2[i] * b.v1

        b.p1 = Param(mutable=True, initialize=5)
        b.p2 = Param(b.s, mutable=True, initialize=6)

        b.port1 = Port()
        b.port1.add(b.v1, "V1")
        b.port1.add(b.v2, "V2")

        b.port2 = Port()
        b.port2.add(b.v1, "V1")
        b.port2.add(b.e2, "V2")

        b.port3 = Port()
        b.port3.add(b.e1, "V1")
        b.port3.add(b.v2, "V2")

        b.port4 = Port()
        b.port4.add(b.p1, "V1")
        b.port4.add(b.v2, "V2")

        b.port5 = Port()
        b.port5.add(b.v1, "V1")
        b.port5.add(b.p2, "V2")

        b.port6 = Port()
        b.port6.add(b.v1, "V1")
        b.port6.add(b.v1, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.port1, destination=m.b2.port1)
    m.s2 = Arc(source=m.b1.port1, destination=m.b2.port2)
    m.s3 = Arc(source=m.b1.port1, destination=m.b2.port3)
    m.s4 = Arc(source=m.b1.port1, destination=m.b2.port4)
    m.s5 = Arc(source=m.b1.port1, destination=m.b2.port5)
    m.s6 = Arc(source=m.b1.port2, destination=m.b2.port1)
    m.s7 = Arc(source=m.b1.port3, destination=m.b2.port1)
    m.s8 = Arc(source=m.b1.port4, destination=m.b2.port1)
    m.s9 = Arc(source=m.b1.port5, destination=m.b2.port1)
    m.s10 = Arc(source=m.b1.port6, destination=m.b2.port1)
    m.s11 = Arc(source=m.b2.port6, destination=m.b1.port1)

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

    with pytest.raises(TypeError):
        propagate_state(m.s2)

    with pytest.raises(TypeError):
        propagate_state(m.s3)

    with pytest.raises(TypeError):
        propagate_state(m.s4)

    with pytest.raises(TypeError):
        propagate_state(m.s5)

    propagate_state(m.s6)
    assert value(m.b1.v1) == value(m.b2.v1)
    assert value(m.b1.e2[1]) == value(m.b2.v2[1])
    assert value(m.b1.e2[2]) == value(m.b2.v2[2])

    propagate_state(m.s7)
    assert value(m.b1.e1) == value(m.b2.v1)
    assert value(m.b1.v2[1]) == value(m.b2.v2[1])
    assert value(m.b1.v2[2]) == value(m.b2.v2[2])

    propagate_state(m.s8)
    assert value(m.b1.p1) == value(m.b2.v1)
    assert value(m.b1.v2[1]) == value(m.b2.v2[1])
    assert value(m.b1.v2[2]) == value(m.b2.v2[2])

    propagate_state(m.s9)
    assert value(m.b1.v1) == value(m.b2.v1)
    assert value(m.b1.p2[1]) == value(m.b2.v2[1])
    assert value(m.b1.p2[2]) == value(m.b2.v2[2])

    with pytest.raises(KeyError):
        propagate_state(m.s10)

    with pytest.raises(KeyError):
        propagate_state(m.s11)


@pytest.mark.unit
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


@pytest.mark.unit
def test_propagate_state_indexed_fixed():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port(b.s)
        b.p[1].add(b.v1, "V1")
        b.p[2].add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    def arc_rule(m, i):
        return {"source": m.b1.p[i], "destination": m.b2.p[i]}

    m.s1 = Arc([1, 2], rule=arc_rule)

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

    propagate_state(m.s1[1])

    # Check that values were propagated correctly
    assert m.b2.v1.value == 500
    assert m.b1.v1.fixed is False
    assert m.b2.v1.fixed is True

    propagate_state(m.s1[1], overwrite_fixed=True)

    # Check that values were propagated correctly
    assert m.b2.v1.value == 10
    assert m.b1.v1.fixed is False
    assert m.b2.v1.fixed is True

    propagate_state(m.s1[2])

    # Check that values were propagated correctly
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


@pytest.mark.unit
def test_propagate_state_indexed():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port(b.s)
        b.p[1].add(b.v1, "V1")
        b.p[2].add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    def arc_rule(m, i):
        return {"source": m.b1.p[i], "destination": m.b2.p[i]}

    m.s1 = Arc([1, 2], rule=arc_rule)

    with pytest.raises(AttributeError):
        propagate_state(m.s1)


@pytest.mark.unit
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


@pytest.mark.unit
def test_propagate_state_invalid_stream():
    m = ConcreteModel()

    with pytest.raises(RuntimeError):
        propagate_state(m)


@pytest.mark.unit
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
@pytest.mark.unit
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
@pytest.mark.unit
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


@pytest.mark.unit
def test_solve_indexed_block_error():
    # Try solve_indexed_block on non-block object
    with pytest.raises(TypeError):
        solve_indexed_blocks(solver=None, blocks=[1, 2, 3])


@pytest.mark.integration
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize_by_time_element():
    horizon = 6
    time_set = [0, horizon]
    ntfe = 60  # For a finite element every six seconds
    ntcp = 2
    m = ConcreteModel(name="CSTR model for testing")
    m.fs = FlowsheetBlock(dynamic=True, time_set=time_set, time_units=pyunits.minute)

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(property_package=m.fs.properties)
    m.fs.cstr = CSTR(
        property_package=m.fs.properties,
        reaction_package=m.fs.reactions,
        material_balance_type=MaterialBalanceType.componentTotal,
        energy_balance_type=EnergyBalanceType.enthalpyTotal,
        momentum_balance_type=MomentumBalanceType.none,
        has_heat_of_reaction=True,
    )

    # Time discretization
    disc = TransformationFactory("dae.collocation")
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme="LAGRANGE-RADAU")

    # Fix geometry variables
    m.fs.cstr.volume[0].fix(1.0)

    # Fix initial conditions:
    for p, j in m.fs.properties.phase_list * m.fs.properties.component_list:
        if j == "Solvent":
            continue
        m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0)

    # Fix inlet conditions
    # This is a huge hack because I didn't know that the proper way to
    # have multiple inlets to a CSTR was to use a mixer.
    # I 'combine' both my inlet streams before sending them to the CSTR.
    for t, j in m.fs.time * m.fs.properties.component_list:
        if t <= 2:
            if j == "E":
                m.fs.cstr.inlet.conc_mol[t, j].fix(11.91 * 0.1 / 2.2)
            elif j == "S":
                m.fs.cstr.inlet.conc_mol[t, j].fix(12.92 * 2.1 / 2.2)
            elif j != "Solvent":
                m.fs.cstr.inlet.conc_mol[t, j].fix(0)
        elif t <= 4:
            if j == "E":
                m.fs.cstr.inlet.conc_mol[t, j].fix(5.95 * 0.1 / 2.2)
            elif j == "S":
                m.fs.cstr.inlet.conc_mol[t, j].fix(12.92 * 2.1 / 2.2)
            elif j != "Solvent":
                m.fs.cstr.inlet.conc_mol[t, j].fix(0)
        else:
            if j == "E":
                m.fs.cstr.inlet.conc_mol[t, j].fix(8.95 * 0.1 / 2.2)
            elif j == "S":
                m.fs.cstr.inlet.conc_mol[t, j].fix(16.75 * 2.1 / 2.2)
            elif j != "Solvent":
                m.fs.cstr.inlet.conc_mol[t, j].fix(0)

    m.fs.cstr.inlet.conc_mol[:, "Solvent"].fix(1.0)

    m.fs.cstr.inlet.flow_vol.fix(2.2)
    m.fs.cstr.inlet.temperature.fix(300)

    # Fix outlet conditions
    m.fs.cstr.outlet.flow_vol.fix(2.2)
    m.fs.cstr.outlet.temperature[m.fs.time.first()].fix(300)

    assert degrees_of_freedom(m) == 0

    initialize_by_time_element(m.fs, m.fs.time, solver=solver)

    assert degrees_of_freedom(m) == 0

    # Assert that the result looks how we expect
    assert m.fs.cstr.outlet.conc_mol[0, "S"].value == 0
    assert abs(m.fs.cstr.outlet.conc_mol[2, "S"].value - 11.389) < 1e-2
    assert abs(m.fs.cstr.outlet.conc_mol[4, "P"].value - 0.2191) < 1e-3
    assert abs(m.fs.cstr.outlet.conc_mol[6, "E"].value - 0.0327) < 1e-3
    assert abs(m.fs.cstr.outlet.temperature[6].value - 289.7) < 1

    # Assert that model is still fixed and deactivated as expected
    assert m.fs.cstr.control_volume.material_holdup[m.fs.time.first(), "aq", "S"].fixed

    for t in m.fs.time:
        if t != m.fs.time.first():
            assert not m.fs.cstr.control_volume.material_holdup[t, "aq", "S"].fixed

            assert not m.fs.cstr.outlet.temperature[t].fixed
        assert m.fs.cstr.control_volume.material_holdup_calculation[t, "aq", "C"].active

        assert m.fs.cstr.control_volume.properties_out[t].active
        assert not m.fs.cstr.outlet.conc_mol[t, "S"].fixed
        assert m.fs.cstr.inlet.conc_mol[t, "S"].fixed

    # Assert that constraints are feasible after initialization
    for con in m.fs.component_data_objects(Constraint, active=True):
        assert value(con.body) - value(con.upper) < 1e-5
        assert value(con.lower) - value(con.body) < 1e-5

    results = solver.solve(m.fs)
    assert check_optimal_termination(results)
