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
This module contains utility functions for use in testing IDAES models.
"""

__author__ = "Andrew Lee"


from pyomo.environ import Constraint, Set, units, Var
from pyomo.common.config import ConfigBlock

from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    ReactionParameterBlock,
    ReactionBlockBase,
    ReactionBlockDataBase,
    MaterialFlowBasis,
    MaterialBalanceType,
    EnergyBalanceType,
    Component,
    Phase,
)

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def initialization_tester(m, dof=0, unit=None, **init_kwargs):
    """
    A method to test initialization methods on IDAES models. This method is
    designed to be used as part of the tests for most models.

    This method checks that the initialization methods runs as expected
    and that the state of the model (active/inactive and fixed/unfixed)
    remains the same.

    This method also add some dummy constraints to the model and deactivates
    them to make sure that the initialization does not affect their status.

    Args:
        m: a Concrete model which contains a flowsheet and a model named unit
            (i.e. m.fs.unit) which will be initialized
        dof: expected degrees of freedom during initialization, default=0
        unit: unit object to test, if None assume m.fs.unit, default='None'
        init_kwargs: model specific arguments to pass to initialize method
                     (e.g. initial guesses for states)

    Returns:
        None

    Raises:
        AssertionErrors if an issue is found
    """
    if unit is None:
        unit = m.fs.unit
    # Add some extra constraints and deactivate them to make sure
    # they remain deactivated
    # Test both indexed and unindexed constraints
    unit.__dummy_var = Var()
    unit.__dummy_equality = Constraint(expr=unit.__dummy_var == 5)
    unit.__dummy_inequality = Constraint(expr=unit.__dummy_var <= 10)

    def deq_idx(b, i):
        return unit.__dummy_var == 5

    unit.__dummy_equality_idx = Constraint([1], rule=deq_idx)

    def dieq_idx(b, i):
        return unit.__dummy_var <= 10

    unit.__dummy_inequality_idx = Constraint([1], rule=dieq_idx)

    unit.__dummy_equality.deactivate()
    unit.__dummy_inequality.deactivate()
    unit.__dummy_equality_idx[1].deactivate()
    unit.__dummy_inequality_idx[1].deactivate()

    orig_fixed_vars = fixed_variables_set(m)
    orig_act_consts = activated_constraints_set(m)

    unit.initialize(**init_kwargs)

    assert degrees_of_freedom(m) == dof

    fin_fixed_vars = fixed_variables_set(m)
    fin_act_consts = activated_constraints_set(m)

    assert len(fin_act_consts) == len(orig_act_consts)
    assert len(fin_fixed_vars) == len(orig_fixed_vars)

    for c in fin_act_consts:
        assert c in orig_act_consts
    for v in fin_fixed_vars:
        assert v in orig_fixed_vars

    # Check dummy constraints and clean up
    assert not unit.__dummy_equality.active
    assert not unit.__dummy_inequality.active
    assert not unit.__dummy_equality_idx[1].active
    assert not unit.__dummy_inequality_idx[1].active

    unit.del_component(unit.__dummy_inequality)
    unit.del_component(unit.__dummy_equality)
    unit.del_component(unit.__dummy_inequality_idx)
    unit.del_component(unit.__dummy_equality_idx)
    unit.del_component(unit.__dummy_var)


# -----------------------------------------------------------------------------
# Define some generic PhysicalBlock and ReactionBlock classes for testing
@declare_process_block_class("PhysicalParameterTestBlock")
class _PhysicalParameterBlock(PhysicalParameterBlock):
    def build(self):
        super(_PhysicalParameterBlock, self).build()

        self.p1 = Phase()
        self.p2 = Phase()

        self.c1 = Component()
        self.c2 = Component()

        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])
        self.element_list = Set(initialize=["H", "He", "Li"])
        self.element_comp = {
            "c1": {"H": 1, "He": 2, "Li": 3},
            "c2": {"H": 4, "He": 5, "Li": 6},
        }

        self.phase_equilibrium_list = {
            "e1": ["c1", ("p1", "p2")],
            "e2": ["c2", ("p1", "p2")],
        }

        # Add inherent reactions for use when needed
        self.inherent_reaction_idx = Set(initialize=["i1", "i2"])
        self.inherent_reaction_stoichiometry = {
            ("i1", "p1", "c1"): 1,
            ("i1", "p1", "c2"): 1,
            ("i1", "p2", "c1"): 1,
            ("i1", "p2", "c2"): 1,
            ("i2", "p1", "c1"): 1,
            ("i2", "p1", "c2"): 1,
            ("i2", "p2", "c1"): 1,
            ("i2", "p2", "c2"): 1,
        }

        # Attribute to switch flow basis for testing
        self.basis_switch = 1
        self.default_balance_switch = 1

        self._state_block_class = TestStateBlock

        self.set_default_scaling("flow_vol", 100)
        self.set_default_scaling("flow_mol", 101)
        self.set_default_scaling("flow_mol_phase_comp", 102)
        self.set_default_scaling("test_var", 103)
        self.set_default_scaling("pressure", 104)
        self.set_default_scaling("temperature", 105)
        self.set_default_scaling("enth_mol", 106)
        self.set_default_scaling("gibbs_mol_phase_comp", 107)
        self.set_default_scaling("entr_mol", 108)
        self.set_default_scaling("mole_frac_phase_comp", 109)
        self.set_default_scaling("enthalpy_flow", 110)
        self.set_default_scaling("energy_dens", 111)
        self.set_default_scaling("material_flow_mol", 112)
        self.set_default_scaling("material_dens_mol", 113)
        self.set_default_scaling("material_flow_mass", 114)
        self.set_default_scaling("material_dens_mass", 115)

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


class SBlockBase(StateBlock):
    def initialize(
        blk,
        outlvl=idaeslog.NOTSET,
        optarg=None,
        solver=None,
        hold_state=False,
        state_args=None,
    ):
        for k in blk.keys():
            blk[k].init_test = True
            blk[k].hold_state = hold_state

    def release_state(blk, flags=None, outlvl=idaeslog.NOTSET):
        for k in blk.keys():
            blk[k].hold_state = not blk[k].hold_state


@declare_process_block_class("TestStateBlock", block_class=SBlockBase)
class StateTestBlockData(StateBlockData):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(StateTestBlockData, self).build()

        self.flow_vol = Var(initialize=20, units=units.m**3 / units.s)
        self.flow_mol = Var(
            initialize=14,
            units=units.mol / units.s,
        )
        self.flow_mol_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=2,
            units=units.mol / units.s,
        )
        self.test_var = Var(initialize=1)
        self.enthalpy_flow = Var(initialize=1, units=units.J / units.s)
        self.energy_dens = Var(initialize=1, units=units.J / units.m**3)
        self.material_flow_mol = Var(initialize=1, units=units.mol / units.s)
        self.material_dens_mol = Var(initialize=1, units=units.mol / units.m**3)
        self.material_flow_mass = Var(initialize=1, units=units.kg / units.s)
        self.material_dens_mass = Var(initialize=1, units=units.kg / units.m**3)
        self.pressure = Var(initialize=1e5, units=units.Pa)
        self.temperature = Var(initialize=300, units=units.K)

        self.cp_mol = Var(initialize=100, units=units.J / units.mol / units.K)
        self.enth_mol = Var(initialize=10000, units=units.J / units.mol)

        self.gibbs_mol_phase_comp = Var(
            self.params.phase_list,
            self.params.component_list,
            initialize=50,
            units=units.J / units.mol,
        )
        self.entr_mol = Var(initialize=1000, units=units.J / units.mol / units.K)

        self.mole_frac_phase_comp = Var(
            self.params.phase_list, self.params.component_list, initialize=0.5
        )

    def get_material_flow_terms(b, p, j):
        if b.config.parameters.basis_switch == 2:
            return b.material_flow_mass
        else:
            return b.material_flow_mol

    def get_material_density_terms(b, p, j):
        if b.config.parameters.basis_switch == 2:
            return b.material_dens_mass
        else:
            return b.material_dens_mol

    def get_enthalpy_flow_terms(b, p):
        return b.enthalpy_flow

    def get_energy_density_terms(b, p):
        return b.energy_dens

    def model_check(self):
        self.check = True

    def get_material_flow_basis(b):
        if b.config.parameters.basis_switch == 1:
            return MaterialFlowBasis.molar
        elif b.config.parameters.basis_switch == 2:
            return MaterialFlowBasis.mass
        else:
            return MaterialFlowBasis.other

    def default_material_balance_type(self):
        if self.params.default_balance_switch == 1:
            return MaterialBalanceType.componentPhase
        else:
            raise NotImplementedError

    def default_energy_balance_type(self):
        if self.params.default_balance_switch == 1:
            return EnergyBalanceType.enthalpyTotal
        else:
            raise NotImplementedError

    def define_state_vars(self):
        return {
            "component_flow_phase": self.flow_mol_phase_comp,
            "temperature": self.temperature,
            "pressure": self.pressure,
        }


@declare_process_block_class("ReactionParameterTestBlock")
class _ReactionParameterBlock(ReactionParameterBlock):
    def build(self):
        super(_ReactionParameterBlock, self).build()

        self.rate_reaction_idx = Set(initialize=["r1", "r2"])
        self.equilibrium_reaction_idx = Set(initialize=["e1", "e2"])

        self.rate_reaction_stoichiometry = {
            ("r1", "p1", "c1"): 1,
            ("r1", "p1", "c2"): 1,
            ("r1", "p2", "c1"): 1,
            ("r1", "p2", "c2"): 1,
            ("r2", "p1", "c1"): 1,
            ("r2", "p1", "c2"): 1,
            ("r2", "p2", "c1"): 1,
            ("r2", "p2", "c2"): 1,
        }
        self.equilibrium_reaction_stoichiometry = {
            ("e1", "p1", "c1"): 1,
            ("e1", "p1", "c2"): 1,
            ("e1", "p2", "c1"): 1,
            ("e1", "p2", "c2"): 1,
            ("e2", "p1", "c1"): 1,
            ("e2", "p1", "c2"): 1,
            ("e2", "p2", "c1"): 1,
            ("e2", "p2", "c2"): 1,
        }

        self._reaction_block_class = ReactionBlock

        # Attribute to switch flow basis for testing
        self.basis_switch = 1

        self.set_default_scaling("reaction_rate", 101, "r1")
        self.set_default_scaling("reaction_rate", 102, "r2")

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )

    @classmethod
    def get_required_properties(self):
        return {}


class RBlockBase(ReactionBlockBase):
    def initialize(
        blk, outlvl=idaeslog.NOTSET, optarg=None, solver=None, state_vars_fixed=False
    ):
        for k in blk.keys():
            blk[k].init_test = True


@declare_process_block_class("ReactionBlock", block_class=RBlockBase)
class ReactionBlockData(ReactionBlockDataBase):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(ReactionBlockData, self).build()

        self.reaction_rate = Var(["r1", "r2"], units=units.mol / units.m**3 / units.s)

        self.dh_rxn = {
            "r1": 10 * units.J / units.mol,
            "r2": 20 * units.J / units.mol,
            "e1": 30 * units.J / units.mol,
            "e2": 40 * units.J / units.mol,
        }

    def model_check(self):
        self.check = True

    def get_reaction_rate_basis(b):
        if b.config.parameters.basis_switch == 1:
            return MaterialFlowBasis.molar
        elif b.config.parameters.basis_switch == 2:
            return MaterialFlowBasis.mass
        else:
            return MaterialFlowBasis.other
