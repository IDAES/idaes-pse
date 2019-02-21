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
Tests for Ideal + Liquid activity coefficient state block;
only tests for construction.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import ConcreteModel

from idaes.core import FlowsheetBlock
from idaes.property_models.NRTL.BTX_ideal_VLE import IdealParameterBlock
from idaes.ui.report import degrees_of_freedom

# -----------------------------------------------------------------------------
# Create a flowsheet for test
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

# vapor-liquid (NRTL)
m.fs.properties_NRTL = IdealParameterBlock(default={"valid_phase":
                                                    ('Liq', 'Vap'),
                                                    "activity_coeff_model":
                                                    'NRTL'})
m.fs.state_block_NRTL = m.fs.properties_NRTL.state_block_class(
    default={"parameters": m.fs.properties_NRTL,
             "defined_state": True})

# vapor-liquid (Wilson)
m.fs.properties_Wilson = IdealParameterBlock(default={"valid_phase":
                                                      ('Liq', 'Vap'),
                                                      "activity_coeff_model":
                                                      'Wilson'})
m.fs.state_block_Wilson = m.fs.properties_Wilson.state_block_class(
    default={"parameters": m.fs.properties_Wilson,
             "defined_state": True})


def test_build_vap_liq_inlet():
    assert len(m.fs.properties_NRTL.config) == 3

    # vapor-liquid (NRTL)
    assert m.fs.properties_NRTL.config.valid_phase == ('Vap', 'Liq') or \
        m.fs.properties_NRTL.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs.properties_NRTL.phase_list) == 2
    assert m.fs.properties_NRTL.phase_list == ["Liq", "Vap"]
    assert m.fs.state_block_NRTL.config.defined_state
    assert hasattr(m.fs.state_block_NRTL, "eq_Keq")
    assert hasattr(m.fs.state_block_NRTL, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_NRTL, "eq_mol_frac_out")

    # vapor-liquid (Wilson)
    assert len(m.fs.properties_Wilson.config) == 3

    assert m.fs.properties_Wilson.config.valid_phase == ('Vap', 'Liq') or \
        m.fs.properties_Wilson.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs.properties_Wilson.phase_list) == 2
    assert m.fs.properties_Wilson.phase_list == ["Liq", "Vap"]
    assert hasattr(m.fs.state_block_Wilson, "eq_Keq")
    assert hasattr(m.fs.state_block_Wilson, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_NRTL, "eq_mol_frac_out")


def test_setInputs_vap_liq_inlet():

    # vapor-liquid (NRTL)
    m.fs.state_block_NRTL.flow_mol.fix(1)
    m.fs.state_block_NRTL.temperature.fix(368)
    m.fs.state_block_NRTL.pressure.fix(101325)
    m.fs.state_block_NRTL.mole_frac["benzene"].fix(0.5)
    m.fs.state_block_NRTL.mole_frac["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_NRTL) == 6

    # vapor-liquid (Wilson)
    m.fs.state_block_Wilson.flow_mol.fix(1)
    m.fs.state_block_Wilson.temperature.fix(368)
    m.fs.state_block_Wilson.pressure.fix(101325)
    m.fs.state_block_Wilson.mole_frac["benzene"].fix(0.5)
    m.fs.state_block_Wilson.mole_frac["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_Wilson) == 4


m.fs1 = FlowsheetBlock(default={"dynamic": False})

# vapor-liquid (NRTL)
m.fs1.properties_NRTL = IdealParameterBlock(default={"valid_phase":
                                                     ('Liq', 'Vap'),
                                                     "activity_coeff_model":
                                                     'NRTL'})
m.fs1.state_block_NRTL = m.fs1.properties_NRTL.state_block_class(
    default={"parameters": m.fs1.properties_NRTL,
             "defined_state": False})

# vapor-liquid (Wilson)
m.fs1.properties_Wilson = IdealParameterBlock(default={"valid_phase":
                                                       ('Liq', 'Vap'),
                                                       "activity_coeff_model":
                                                       'Wilson'})
m.fs1.state_block_Wilson = m.fs1.properties_Wilson.state_block_class(
    default={"parameters": m.fs1.properties_Wilson,
             "defined_state": False})


def test_build_vap_liq_outlet():
    assert len(m.fs.properties_NRTL.config) == 3

    # vapor-liquid (NRTL)
    assert m.fs1.properties_NRTL.config.valid_phase == ('Vap', 'Liq') or \
        m.fs1.properties_NRTL.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs1.properties_NRTL.phase_list) == 2
    assert m.fs1.properties_NRTL.phase_list == ["Liq", "Vap"]
    assert not m.fs1.state_block_NRTL.config.defined_state
    assert hasattr(m.fs1.state_block_NRTL, "eq_Keq")
    assert hasattr(m.fs1.state_block_NRTL, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_NRTL, "eq_mol_frac_out")

    # vapor-liquid (Wilson)
    assert len(m.fs1.properties_Wilson.config) == 3

    assert m.fs1.properties_Wilson.config.valid_phase == ('Vap', 'Liq') or \
        m.fs1.properties_Wilson.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs1.properties_Wilson.phase_list) == 2
    assert m.fs1.properties_Wilson.phase_list == ["Liq", "Vap"]
    assert not m.fs1.state_block_NRTL.config.defined_state
    assert hasattr(m.fs1.state_block_Wilson, "eq_Keq")
    assert hasattr(m.fs1.state_block_Wilson, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_NRTL, "eq_mol_frac_out")


def test_setInputs_vap_liq_outlet():

    # vapor-liquid (NRTL)
    m.fs1.state_block_NRTL.flow_mol.fix(1)
    m.fs1.state_block_NRTL.temperature.fix(368)
    m.fs1.state_block_NRTL.pressure.fix(101325)
    m.fs1.state_block_NRTL.mole_frac["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_NRTL) == 6

    # vapor-liquid (Wilson)
    m.fs1.state_block_Wilson.flow_mol.fix(1)
    m.fs1.state_block_Wilson.temperature.fix(368)
    m.fs1.state_block_Wilson.pressure.fix(101325)
    m.fs1.state_block_Wilson.mole_frac["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_Wilson) == 4
