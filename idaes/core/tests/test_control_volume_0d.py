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
from pyomo.environ import ConcreteModel, Set
from pyomo.common.config import ConfigBlock
from idaes.core import (ControlVolume0D, ControlVolumeBase, FlowsheetBlockData,
                        declare_process_block_class, useDefault, FlowDirection,
                        PhysicalParameterBase, StateBlockDataBase)


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("PropertyParameterBlock")
class _PropertyParameterBlock(PhysicalParameterBase):
    def build(self):
        super(_PropertyParameterBlock, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])


@declare_process_block_class("StateBlock")
class StateBlockData(StateBlockDataBase):
    CONFIG = ConfigBlock(implicit=True)

    def build(self):
        super(StateBlockData, self).build()


@declare_process_block_class("CVFrame")
class CVFrameData(ControlVolume0D):
    def build(self):
        super(ControlVolumeBase, self).build()


# -----------------------------------------------------------------------------
# Basic tests
def test_base_build():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PropertyParameterBlock()

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
# Test add_state_blocks
def test_add_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PropertyParameterBlock()

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
    m.fs.pp = PropertyParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(information_flow=FlowDirection.forward)

    assert m.fs.cv.properties_in[0].config.defined_state is True
    assert m.fs.cv.properties_out[0].config.defined_state is False


def test_add_state_block_backward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PropertyParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(information_flow=FlowDirection.backward)

    assert m.fs.cv.properties_in[0].config.defined_state is False
    assert m.fs.cv.properties_out[0].config.defined_state is True


def test_add_state_blocks_has_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PropertyParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    assert m.fs.cv.properties_in[0].config.has_phase_equilibrium is True
    assert m.fs.cv.properties_out[0].config.has_phase_equilibrium is True


def test_add_state_blocks_custom_args():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs.pp = PropertyParameterBlock()

    m.fs.cv = ControlVolume0D(default={"property_package": m.fs.pp})

    m.fs.cv.add_state_blocks(package_arguments={"test": "test"})

    assert len(m.fs.cv.properties_in[0].config) == 4
    assert m.fs.cv.properties_in[0].config.test == "test"

    assert len(m.fs.cv.properties_out[0].config) == 4
    assert m.fs.cv.properties_out[0].config.test == "test"


# -----------------------------------------------------------------------------
# Test add_reaction_blocks
