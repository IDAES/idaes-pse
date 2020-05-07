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
Tests for ControlVolumeBlockData, and for initializing the moving bed module

Author: Chinedu Okoli
"""

import sys
import os
import pytest

from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.gas_solid_contactors.unit_models.moving_bed import MovingBed
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import Gas_Phase_Thermo_ParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import Solid_Phase_Thermo_ParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    hetero_reactions import HeteroReactionParameterBlock


# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
else:
    solver = None


# -----------------------------------------------------------------------------
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
        default={"solid_property_package": m.fs.solid_properties,
                 "gas_property_package": m.fs.gas_properties})

    m.fs.MB = MovingBed(default={
        "finite_elements": 10,
        "transformation_method":
        "dae.finite_difference",
        "transformation_scheme": "BACKWARD",
        "energy_balance_type": EnergyBalanceType.enthalpyTotal,
        "gas_phase_config": {
            "property_package": m.fs.gas_properties,
            "has_pressure_change": True,
            "pressure_drop_type": "ergun_correlation"},
        "solid_phase_config": {
            "property_package": m.fs.solid_properties,
            "reaction_package": m.fs.hetero_reactions
            }})

    assert m.fs.MB.config.flow_type == "counter_current"
    assert m.fs.MB.config.collocation_points == 3
    assert m.fs.MB.config.finite_elements == 10

    assert m.fs.MB.config.transformation_scheme == "BACKWARD"

    # There should be 13 DOFs in this model:
    # Geometry - 2 (bed length and height)
    # Gas feed - 6 (inlet flow, temperature, pressure and mole fractions(3))
    # Solid feed - 5 (inlet flow, temperature, and mass fractions (3))
    assert degrees_of_freedom(m) == 13


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
        default={"solid_property_package": m.fs.solid_properties,
                 "gas_property_package": m.fs.gas_properties})

    m.fs.MB = MovingBed(default={
        "finite_elements": 10,
        "transformation_method":
            "dae.finite_difference",
        "transformation_scheme": "BACKWARD",
        "energy_balance_type": EnergyBalanceType.enthalpyTotal,
        "gas_phase_config": {
            "property_package": m.fs.gas_properties,
            "has_pressure_change": True,
            "pressure_drop_type": "ergun_correlation"},
        "solid_phase_config": {
            "property_package": m.fs.solid_properties,
            "reaction_package": m.fs.hetero_reactions
            }})

    # Fix bed geometry variables
    m.fs.MB.bed_diameter.fix(6.5)  # m
    m.fs.MB.bed_height.fix(5)  # m

    # Fix inlet port variables for gas and solid
    m.fs.MB.gas_inlet.flow_mol[0].fix(128.20513)  # mol/s
    m.fs.MB.gas_inlet.temperature[0].fix(298.15)  # K
    m.fs.MB.gas_inlet.pressure[0].fix(2.00)  # bar
    m.fs.MB.gas_inlet.mole_frac[0, "CO2"].fix(0.02499)
    m.fs.MB.gas_inlet.mole_frac[0, "H2O"].fix(0.00001)
    m.fs.MB.gas_inlet.mole_frac[0, "CH4"].fix(0.975)

    m.fs.MB.solid_inlet.flow_mass[0].fix(591.4)  # kg/s
    m.fs.MB.solid_inlet.temperature[0].fix(1183.15)  # K
    m.fs.MB.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.45)
    m.fs.MB.solid_inlet.mass_frac[0, "Fe3O4"].fix(1e-9)
    m.fs.MB.solid_inlet.mass_frac[0, "Al2O3"].fix(0.55)

    assert degrees_of_freedom(m.fs) == 0

    # State arguments for initializing property state blocks
    # Gas phase temperature is initialized at solid
    # temperature because thermal mass of solid >> thermal mass of gas
    blk = m.fs.MB
    gas_phase_state_args = {
        'flow_mol': blk.gas_inlet.flow_mol[0].value,
        'temperature': blk.solid_inlet.temperature[0].value,
        'pressure': blk.gas_inlet.pressure[0].value,
        'mole_frac': {
            'CH4': blk.gas_inlet.mole_frac[0, 'CH4'].value,
            'CO2': blk.gas_inlet.mole_frac[0, 'CO2'].value,
            'H2O': blk.gas_inlet.mole_frac[0, 'H2O'].value}}
    solid_phase_state_args = {
        'flow_mass': blk.solid_inlet.flow_mass[0].value,
        'temperature': blk.solid_inlet.temperature[0].value,
        'mass_frac': {
            'Fe2O3': blk.solid_inlet.mass_frac[0, 'Fe2O3'].value,
            'Fe3O4': blk.solid_inlet.mass_frac[0, 'Fe3O4'].value,
            'Al2O3': blk.solid_inlet.mass_frac[0, 'Al2O3'].value}}

    m.fs.MB.initialize(outlvl=5,
                       gas_phase_state_args=gas_phase_state_args,
                       solid_phase_state_args=solid_phase_state_args)

    assert degrees_of_freedom(m.fs) == 0
