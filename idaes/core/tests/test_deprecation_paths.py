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
Tests for deprecation paths for moved base classes
"""
import pytest


@pytest.mark.unit
def test_import_rxn_param_block(caplog):
    from idaes.core.reaction_base import ReactionParameterBlock
    assert (
        "DEPRECATED: the 'ReactionParameterBlock' class has been moved to"
        "'idaes.core.base.reaction_base.ReactionParameterBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_rxn_block_base(caplog):
    from idaes.core.reaction_base import ReactionBlockBase
    assert (
        "DEPRECATED: the 'ReactionBlockBase' class has been moved to"
        "'idaes.core.base.reaction_base.ReactionBlockBase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_rxn_block_data_base(caplog):
    from idaes.core.reaction_base import ReactionBlockDataBase
    assert (
        "DEPRECATED: the 'ReactionBlockDataBase' class has been moved to"
        "'idaes.core.base.reaction_base.ReactionBlockDataBase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_mat_flow_basis(caplog):
    from idaes.core.process_base import MaterialFlowBasis
    assert (
        "DEPRECATED: the 'MaterialFlowBasis' class has been moved to"
        "'idaes.core.base.process_base.MaterialFlowBasis'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_proc_block_data(caplog):
    from idaes.core.process_base import ProcessBlockData
    assert (
        "DEPRECATED: the 'ProcessBlockData' class has been moved to"
        "'idaes.core.base.process_base.ProcessBlockData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_phys_param_block(caplog):
    from idaes.core.property_base import PhysicalParameterBlock
    assert (
        "DEPRECATED: the 'PhysicalParameterBlock' class has been moved to"
        "'idaes.core.base.property_base.PhysicalParameterBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_state_block(caplog):
    from idaes.core.property_base import StateBlock
    assert (
        "DEPRECATED: the 'StateBlock' class has been moved to"
        "'idaes.core.base.property_base.StateBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_state_block_data(caplog):
    from idaes.core.property_base import StateBlockData
    assert (
        "DEPRECATED: the 'StateBlockData' class has been moved to"
        "'idaes.core.base.property_base.StateBlockData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_mat_bal_type(caplog):
    from idaes.core.control_volume_base import MaterialBalanceType
    assert (
        "DEPRECATED: the 'MaterialBalanceType' class has been moved to"
        "'idaes.core.base.control_volume_base.MaterialBalanceType'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_e_bal_type(caplog):
    from idaes.core.control_volume_base import EnergyBalanceType
    assert (
        "DEPRECATED: the 'EnergyBalanceType' class has been moved to"
        "'idaes.core.base.control_volume_base.EnergyBalanceType'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_mom_bal_type(caplog):
    from idaes.core.control_volume_base import MomentumBalanceType
    assert (
        "DEPRECATED: the 'MomentumBalanceType' class has been moved to"
        "'idaes.core.base.control_volume_base.MomentumBalanceType'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_config_template(caplog):
    from idaes.core.control_volume_base import CONFIG_Template
    assert (
        "DEPRECATED: the 'CONFIG_Template' attribute has been moved to"
        "'idaes.core.base.control_volume_base.CONFIG_Template'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_cv_block_data(caplog):
    from idaes.core.control_volume_base import ControlVolumeBlockData
    assert (
        "DEPRECATED: the 'ControlVolumeBlockData' class has been moved to"
        "'idaes.core.base.control_volume_base.ControlVolumeBlockData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_flowsheet_block(caplog):
    from idaes.core.flowsheet_model import FlowsheetBlock
    assert (
        "DEPRECATED: the 'FlowsheetBlock' class has been moved to"
        "'idaes.core.base.flowsheet_model.FlowsheetBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_flowsheet_block_data(caplog):
    from idaes.core.flowsheet_model import FlowsheetBlockData
    assert (
        "DEPRECATED: the 'FlowsheetBlockData' class has been moved to"
        "'idaes.core.base.flowsheet_model.FlowsheetBlockData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_cv1d_block(caplog):
    from idaes.core.control_volume1d import ControlVolume1DBlock
    assert (
        "DEPRECATED: the 'ControlVolume1DBlock' class has been moved to"
        "'idaes.core.base.control_volume1d.ControlVolume1DBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_cv1d_block_data(caplog):
    from idaes.core.control_volume1d import ControlVolume1DBlockData
    assert (
        "DEPRECATED: the 'ControlVolume1DBlockData' class has been moved to"
        "'idaes.core.base.control_volume1d.ControlVolume1DBlockData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_dist_vars(caplog):
    from idaes.core.control_volume1d import DistributedVars
    assert (
        "DEPRECATED: the 'DistributedVars' class has been moved to"
        "'idaes.core.base.control_volume1d.DistributedVars'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_component(caplog):
    from idaes.core.components import Component
    assert (
        "DEPRECATED: the 'Component' class has been moved to"
        "'idaes.core.base.components.Component'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_solute(caplog):
    from idaes.core.components import Solute
    assert (
        "DEPRECATED: the 'Solute' class has been moved to"
        "'idaes.core.base.components.Solute'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_solvent(caplog):
    from idaes.core.components import Solvent
    assert (
        "DEPRECATED: the 'Solvent' class has been moved to"
        "'idaes.core.base.components.Solvent'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_ion(caplog):
    from idaes.core.components import Ion
    assert (
        "DEPRECATED: the 'Ion' class has been moved to"
        "'idaes.core.base.components.Ion'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_anion(caplog):
    from idaes.core.components import Anion
    assert (
        "DEPRECATED: the 'Anion' class has been moved to"
        "'idaes.core.base.components.Anion'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_cation(caplog):
    from idaes.core.components import Cation
    assert (
        "DEPRECATED: the 'Cation' class has been moved to"
        "'idaes.core.base.components.Cation'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_apparent(caplog):
    from idaes.core.components import Apparent
    assert (
        "DEPRECATED: the 'Apparent' class has been moved to"
        "'idaes.core.base.components.Apparent'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_cv0d_block(caplog):
    from idaes.core.control_volume0d import ControlVolume0DBlock
    assert (
        "DEPRECATED: the 'ControlVolume0DBlock' class has been moved to"
        "'idaes.core.base.control_volume0d.ControlVolume0DBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_cv0d_block_data(caplog):
    from idaes.core.control_volume0d import ControlVolume0DBlockData
    assert (
        "DEPRECATED: the 'ControlVolume0DBlockData' class has been moved to"
        "'idaes.core.base.control_volume0d.ControlVolume0DBlockData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_phase_type(caplog):
    from idaes.core.phases import PhaseType
    assert (
        "DEPRECATED: the 'PhaseType' class has been moved to"
        "'idaes.core.base.phases.PhaseType'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_phase(caplog):
    from idaes.core.phases import Phase
    assert (
        "DEPRECATED: the 'Phase' class has been moved to"
        "'idaes.core.base.phases.Phase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_phase_data(caplog):
    from idaes.core.phases import PhaseData
    assert (
        "DEPRECATED: the 'PhaseData' class has been moved to"
        "'idaes.core.base.phases.PhaseData'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_liquid_phase(caplog):
    from idaes.core.phases import LiquidPhase
    assert (
        "DEPRECATED: the 'LiquidPhase' class has been moved to"
        "'idaes.core.base.phases.LiquidPhase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_solid_phase(caplog):
    from idaes.core.phases import SolidPhase
    assert (
        "DEPRECATED: the 'SolidPhase' class has been moved to"
        "'idaes.core.base.phases.SolidPhase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_vapor_phase(caplog):
    from idaes.core.phases import VaporPhase
    assert (
        "DEPRECATED: the 'VaporPhase' class has been moved to"
        "'idaes.core.base.phases.VaporPhase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_aqueous_phase(caplog):
    from idaes.core.phases import AqueousPhase
    assert (
        "DEPRECATED: the 'AqueousPhase' class has been moved to"
        "'idaes.core.base.phases.AqueousPhase'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_unit_model_block(caplog):
    from idaes.core.unit_model import UnitModelBlock
    assert (
        "DEPRECATED: the 'UnitModelBlock' class has been moved to"
        "'idaes.core.base.unit_model.UnitModelBlock'"
        ) in caplog.text.replace('\n', '')


@pytest.mark.unit
def test_import_unit_model_block_data(caplog):
    from idaes.core.unit_model import UnitModelBlockData
    assert (
        "DEPRECATED: the 'UnitModelBlockData' class has been moved to"
        "'idaes.core.base.unit_model.UnitModelBlockData'"
        ) in caplog.text.replace('\n', '')
