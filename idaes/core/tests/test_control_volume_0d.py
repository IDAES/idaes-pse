##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for ControlVolumeBase.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Constraint, Expression, Set, Var
from pyomo.common.config import ConfigBlock
from idaes.core import (ControlVolume0D, ControlVolumeBase, FlowsheetBlockData,
                        declare_process_block_class, FlowDirection,
                        PhysicalParameterBase, StateBlockBase,
                        StateBlockDataBase, ReactionParameterBase,
                        ReactionBlockBase, ReactionBlockDataBase)
from idaes.core.util.exceptions import (ConfigurationError,
                                        PropertyNotSupportedError)


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("PhysicalParameterBlock")
class _PhysicalParameterBlock(PhysicalParameterBase):
    def build(self):
        super(_PhysicalParameterBlock, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])
        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])
        self.element_list = Set(initialize=["H", "He", "Li"])

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


@declare_process_block_class("StateBlock", block_class=StateBlockBase)
class StateBlockData(StateBlockDataBase):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(StateBlockData, self).build()

        self.test_var = Var()
        self.phase_equilibrium_idx = Set(initialize=["e1", "e2"])
        self.phase_equilibrium_list = \
            {"e1": ["c1", ("p1", "p2")],
             "e2": ["c2", ("p1", "p2")]}
        self.element_comp = {"c1": {"H": 1, "He": 2, "Li": 3},
                             "c2": {"H": 4, "He": 5, "Li": 6}}

    def get_material_flow_terms(b, p, j):
        return b.test_var

    def get_material_density_terms(b, p, j):
        return b.test_var


@declare_process_block_class("ReactionParameterBlock")
class _ReactionParameterBlock(ReactionParameterBase):
    def build(self):
        super(_ReactionParameterBlock, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])
        self.rate_reaction_idx = Set(initialize=["r1", "r2"])
        self.equilibrium_reaction_idx = Set(initialize=["e1", "e2"])

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})

    @classmethod
    def get_required_properties(self):
        return {}


@declare_process_block_class("ReactionBlock", block_class=ReactionBlockBase)
class ReactionBlockData(ReactionBlockDataBase):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(ReactionBlockData, self).build()

        self.rate_reaction_idx = Set(initialize=["r1", "r2"])
        self.reaction_rate = Var(self.rate_reaction_idx)
        self.rate_reaction_stoichiometry = {("r1", "p1", "c1"): 1,
                                            ("r1", "p1", "c2"): 1,
                                            ("r1", "p2", "c1"): 1,
                                            ("r1", "p2", "c2"): 1,
                                            ("r2", "p1", "c1"): 1,
                                            ("r2", "p1", "c2"): 1,
                                            ("r2", "p2", "c1"): 1,
                                            ("r2", "p2", "c2"): 1}
        self.equilibrium_reaction_idx = Set(initialize=["e1", "e2"])
        self.equilibrium_reaction_stoichiometry = {
                                            ("e1", "p1", "c1"): 1,
                                            ("e1", "p1", "c2"): 1,
                                            ("e1", "p2", "c1"): 1,
                                            ("e1", "p2", "c2"): 1,
                                            ("e2", "p1", "c1"): 1,
                                            ("e2", "p1", "c2"): 1,
                                            ("e2", "p2", "c1"): 1,
                                            ("e2", "p2", "c2"): 1}


@declare_process_block_class("CVFrame")
class CVFrameData(ControlVolume0D):
    def build(self):
        super(ControlVolumeBase, self).build()


# -----------------------------------------------------------------------------
# Basic tests
def test_base_build():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    assert len(m.fs.cv.config) == 6
    assert m.fs.cv.config.dynamic is False
    assert m.fs.cv.config.property_package == m.fs.pp
    assert isinstance(m.fs.cv.config.property_package_args, ConfigBlock)
    assert len(m.fs.cv.config.property_package_args) == 0
    assert m.fs.cv.config.reaction_package is None
    assert isinstance(m.fs.cv.config.reaction_package_args, ConfigBlock)
    assert len(m.fs.cv.config.reaction_package_args) == 0
    assert m.fs.cv.config.auto_construct is False

    assert hasattr(m.fs.cv, "time")
    assert hasattr(m.fs.cv, "phase_list")
    assert hasattr(m.fs.cv, "component_list")


# -----------------------------------------------------------------------------
# Test add_geometry
def test_add_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_geometry()

    assert hasattr(m.fs.cv, "volume")
    assert len(m.fs.cv.volume) == 1.0
    assert m.fs.cv.volume[0].value == 1.0


# -----------------------------------------------------------------------------
# Test add_state_blocks
def test_add_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks()

    assert hasattr(m.fs.cv, "properties_in")
    assert len(m.fs.cv.properties_in[0].config) == 3
    assert m.fs.cv.properties_in[0].config.defined_state is True
    assert m.fs.cv.properties_in[0].config.has_phase_equilibrium is False
    assert m.fs.cv.properties_in[0].config.parameters == m.fs.pp

    assert hasattr(m.fs.cv, "properties_out")
    assert len(m.fs.cv.properties_out[0].config) == 3
    assert m.fs.cv.properties_out[0].config.defined_state is False
    assert m.fs.cv.properties_out[0].config.has_phase_equilibrium is False
    assert m.fs.cv.properties_out[0].config.parameters == m.fs.pp


def test_add_state_block_forward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(information_flow=FlowDirection.forward)

    assert m.fs.cv.properties_in[0].config.defined_state is True
    assert m.fs.cv.properties_out[0].config.defined_state is False


def test_add_state_block_backward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(information_flow=FlowDirection.backward)

    assert m.fs.cv.properties_in[0].config.defined_state is False
    assert m.fs.cv.properties_out[0].config.defined_state is True


def test_add_state_blocks_has_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    assert m.fs.cv.properties_in[0].config.has_phase_equilibrium is True
    assert m.fs.cv.properties_out[0].config.has_phase_equilibrium is True


def test_add_state_blocks_custom_args():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(package_arguments={"test": "test"})

    assert len(m.fs.cv.properties_in[0].config) == 4
    assert m.fs.cv.properties_in[0].config.test == "test"

    assert len(m.fs.cv.properties_out[0].config) == 4
    assert m.fs.cv.properties_out[0].config.test == "test"


# -----------------------------------------------------------------------------
# Test add_reaction_blocks
def test_add_reaction_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    assert hasattr(m.fs.cv, "reactions")
    assert len(m.fs.cv.reactions[0].config) == 3
    assert m.fs.cv.reactions[0].config.state_block == m.fs.cv.properties_out
    assert m.fs.cv.reactions[0]._state == m.fs.cv.properties_out[0]
    assert m.fs.cv.reactions[0].config.has_equilibrium is False
    assert m.fs.cv.reactions[0].config.parameters == m.fs.rp


def test_add_reaction_blocks_has_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    assert m.fs.cv.reactions[0].config.has_equilibrium is True


# -----------------------------------------------------------------------------
# Test _add_phase_fractions
def test_add_phase_fractions():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv._add_phase_fractions()

    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert len(m.fs.cv.phase_fraction) == 2
    assert hasattr(m.fs.cv, "sum_of_phase_fractions")


def test_add_phase_fractions_single_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.pp.del_component(m.fs.pp.phase_list)
    m.fs.pp.phase_list = Set(initialize=["p1"])

    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv._add_phase_fractions()

    assert isinstance(m.fs.cv.phase_fraction, Expression)
    assert len(m.fs.cv.phase_fraction) == 1
    assert not hasattr(m.fs.cv, "sum_of_phase_fractions")


# -----------------------------------------------------------------------------
# Test add_phase_component_balances
def test_add_phase_component_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_phase_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 4


def test_add_phase_component_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_phase_component_balances(dynamic=True, has_holdup=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 8
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.material_holdup, Var)
    assert isinstance(m.fs.cv.material_accumulation, Var)


def test_add_phase_component_balances_dynamic_no_holdup():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(dynamic=True,
                                             has_holdup=False)


def test_add_phase_component_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    # Do not add geometry
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(dynamic=True,
                                             has_holdup=True)


def test_add_phase_component_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_phase_component_balances(has_rate_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.rate_reaction_generation, Var)
    assert isinstance(m.fs.cv.rate_reaction_idx, Set)
    assert isinstance(m.fs.cv.rate_reaction_extent, Var)
    assert isinstance(m.fs.cv.rate_reaction_stoichiometry_constraint,
                      Constraint)
    assert isinstance(m.fs.cv.rate_reaction_extents_constraint,
                      Constraint)


def test_add_phase_component_balances_rate_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_rate_reactions=True)


def test_add_phase_component_balances_rate_rxns_no_rxn_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.rp.del_component(m.fs.rp.rate_reaction_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_rate_reactions=True)


def test_add_phase_component_balances_rate_rxns_no_volume():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    # Do not add geometry
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_rate_reactions=True)


def test_add_phase_component_balances_rate_rxns_no_reaction_rate():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    m.fs.cv.reactions[0].del_component(
            m.fs.cv.reactions[0].reaction_rate)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_rate_reactions=True)


def test_add_phase_component_balances_eq_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    mb = m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.equilibrium_reaction_generation, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_idx, Set)
    assert isinstance(m.fs.cv.equilibrium_reaction_extent, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_stoichiometry_constraint,
                      Constraint)


def test_add_phase_component_balances_eq_rxns_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


def test_add_phase_component_balances_eq_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.rp.del_component(m.fs.rp.equilibrium_reaction_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


def test_add_phase_component_balances_eq_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


def test_add_phase_component_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.phase_equilibrium_generation, Var)
    assert isinstance(m.fs.cv.phase_equilibrium_idx, Set)


def test_add_phase_component_balances_phase_eq_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)


def test_add_phase_component_balances_phase_eq_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)


def test_add_phase_component_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_phase_component_balances(has_mass_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.mass_transfer_term, Var)


# -----------------------------------------------------------------------------
# Test add_total_component_balances
def test_add_total_component_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 2


def test_add_total_component_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_component_balances(dynamic=True, has_holdup=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.material_holdup, Var)
    assert isinstance(m.fs.cv.material_accumulation, Var)


def test_add_total_component_balances_dynamic_no_holdup():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(dynamic=True,
                                             has_holdup=False)


def test_add_total_component_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    # Do not add geometry
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(dynamic=True,
                                             has_holdup=True)


def test_add_total_component_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_component_balances(has_rate_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.rate_reaction_generation, Var)
    assert isinstance(m.fs.cv.rate_reaction_idx, Set)
    assert isinstance(m.fs.cv.rate_reaction_extent, Var)
    assert isinstance(m.fs.cv.rate_reaction_stoichiometry_constraint,
                      Constraint)
    assert isinstance(m.fs.cv.rate_reaction_extents_constraint,
                      Constraint)


def test_add_total_component_balances_rate_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_rate_reactions=True)


def test_add_total_component_balances_rate_rxns_no_rxn_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.rp.del_component(m.fs.rp.rate_reaction_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(has_rate_reactions=True)


def test_add_total_component_balances_rate_rxns_no_volume():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    # Do not add geometry
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_rate_reactions=True)


def test_add_total_component_balances_rate_rxns_no_reaction_rate():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    m.fs.cv.reactions[0].del_component(
            m.fs.cv.reactions[0].reaction_rate)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(has_rate_reactions=True)


def test_add_total_component_balances_eq_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    mb = m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.equilibrium_reaction_generation, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_idx, Set)
    assert isinstance(m.fs.cv.equilibrium_reaction_extent, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_stoichiometry_constraint,
                      Constraint)


def test_add_total_component_balances_eq_rxns_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)


def test_add_total_component_balances_eq_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.rp.del_component(m.fs.rp.equilibrium_reaction_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)


def test_add_total_component_balances_eq_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)


def test_add_total_component_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_component_balances(has_phase_equilibrium=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.phase_equilibrium_generation, Var)
    assert isinstance(m.fs.cv.phase_equilibrium_idx, Set)


def test_add_total_component_balances_phase_eq_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_phase_equilibrium=True)


def test_add_total_component_balances_phase_eq_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(has_phase_equilibrium=True)


def test_add_total_component_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_component_balances(has_mass_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.mass_transfer_term, Var)


# -----------------------------------------------------------------------------
# Test add_total_component_balances
def test_add_total_element_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_element_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 3


def test_add_total_element_balances_properties_not_supported():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.pp.del_component(m.fs.pp.element_list)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_element_balances()


def test_add_total_element_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_element_balances(dynamic=True, has_holdup=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 6
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.element_holdup, Var)
    assert isinstance(m.fs.cv.element_accumulation, Var)


def test_add_total_element_balances_dynamic_no_holdup():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(dynamic=True,
                                           has_holdup=False)


def test_add_total_element_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    # Do not add geometry
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(dynamic=True,
                                           has_holdup=True)


def test_add_total_element_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_rate_reactions=True)


def test_add_total_element_balances_eq_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    mb = m.fs.cv.add_total_element_balances(has_equilibrium_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 3
    assert isinstance(m.fs.cv.equilibrium_reaction_idx, Set)


def test_add_total_element_balances_eq_rxns_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_equilibrium_reactions=True)


def test_add_total_element_balances_eq_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.rp.del_component(m.fs.rp.equilibrium_reaction_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_element_balances(has_equilibrium_reactions=True)


def test_add_total_element_balances_eq_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_equilibrium_reactions=True)


def test_add_total_element_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_element_balances(has_phase_equilibrium=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 3
    assert isinstance(m.fs.cv.phase_equilibrium_idx, Set)


def test_add_total_element_balances_phase_eq_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_phase_equilibrium=True)


def test_add_total_element_balances_phase_eq_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks()

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_element_balances(has_phase_equilibrium=True)


def test_add_total_element_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PhysicalParameterBlock()
    m.fs.rp = ReactionParameterBlock(default={"property_package": m.fs.pp})

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp,
                                       "reaction_package": m.fs.rp})

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks()
    m.fs.cv.add_reaction_blocks()

    mb = m.fs.cv.add_total_element_balances(has_mass_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 3
    assert isinstance(m.fs.cv.elemental_mass_transfer_term, Var)
