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
Test for ControlVolumeBlockData, and for initializing the
bubbling fluidized bed module

Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           SolverFactory,
                           value,
                           Var)
from pyomo.common.config import ConfigBlock
from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_total_constraints,
                                              fixed_variables_set,
                                              activated_constraints_set,
                                              number_unused_variables,
                                              number_derivative_variables,
                                              unused_variables_set)
from idaes.core.util.testing import (get_default_solver,
                                     PhysicalParameterTestBlock,
                                     ReactionParameterTestBlock,
                                     initialization_tester)

from idaes.gas_solid_contactors.unit_models. \
    bubbling_fluidized_bed import BubblingFluidizedBed
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import GasPhaseThermoParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import SolidPhaseThermoParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    hetero_reactions import HeteroReactionParameterBlock

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.gas_properties = GasPhaseThermoParameterBlock()
    m.fs.solid_properties = SolidPhaseThermoParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m.fs.solid_properties,
                     "gas_property_package": m.fs.gas_properties})

    m.fs.unit = BubblingFluidizedBed(
            default={
                    "gas_phase_config":
                    {"property_package": m.fs.gas_properties},
                    "solid_phase_config":
                    {"property_package": m.fs.solid_properties,
                     "reaction_package": m.fs.hetero_reactions
                     }})

    # Check unit config arguments
    assert len(m.fs.unit.config) == 14
    assert isinstance(m.fs.unit.config.gas_phase_config, ConfigBlock)
    assert isinstance(m.fs.unit.config.solid_phase_config, ConfigBlock)

    assert m.fs.unit.config.finite_elements == 10
    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == 'BACKWARD'
    assert m.fs.unit.config.collocation_points == 3
    assert m.fs.unit.config.flow_type == "co_current"
    assert m.fs.unit.config.material_balance_type == \
        MaterialBalanceType.componentTotal
    assert m.fs.unit.config.energy_balance_type == \
        EnergyBalanceType.enthalpyTotal
    assert m.fs.unit.config.momentum_balance_type == \
        MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_pressure_change is True

    # Check gas phase config arguments
    assert len(m.fs.unit.config.gas_phase_config) == 7
    assert m.fs.unit.config.gas_phase_config.has_equilibrium_reactions is False
    assert m.fs.unit.config.gas_phase_config.property_package is \
        m.fs.gas_properties
    assert m.fs.unit.config.gas_phase_config.reaction_package is None

    # Check solid phase config arguments
    assert len(m.fs.unit.config.solid_phase_config) == 7
    assert m.fs.unit.config.solid_phase_config.has_equilibrium_reactions is \
        False
    assert m.fs.unit.config.solid_phase_config.property_package is \
        m.fs.solid_properties
    assert m.fs.unit.config.solid_phase_config.reaction_package is \
        m.fs.hetero_reactions


# -----------------------------------------------------------------------------
class TestIronOC(object):
    @pytest.fixture(scope="class")
    def iron_oc(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        # Set up thermo props and reaction props
        m.fs.gas_properties = GasPhaseThermoParameterBlock()
        m.fs.solid_properties = SolidPhaseThermoParameterBlock()
        m.fs.hetero_reactions = HeteroReactionParameterBlock(
                default={"solid_property_package": m.fs.solid_properties,
                         "gas_property_package": m.fs.gas_properties})

        m.fs.unit = BubblingFluidizedBed(
                default={
                        "gas_phase_config":
                        {"property_package": m.fs.gas_properties},
                        "solid_phase_config":
                        {"property_package": m.fs.solid_properties,
                         "reaction_package": m.fs.hetero_reactions
                         }})

        return m

    @pytest.mark.build
    def test_build(self, iron_oc):
        assert hasattr(iron_oc.fs.unit, "gas_inlet")
        assert len(iron_oc.fs.unit.gas_inlet.vars) == 4
        assert hasattr(iron_oc.fs.unit.gas_inlet, "flow_mol")
        assert hasattr(iron_oc.fs.unit.gas_inlet, "mole_frac")
        assert hasattr(iron_oc.fs.unit.gas_inlet, "temperature")
        assert hasattr(iron_oc.fs.unit.gas_inlet, "pressure")

        assert hasattr(iron_oc.fs.unit, "solid_inlet")
        assert len(iron_oc.fs.unit.solid_inlet.vars) == 3
        assert hasattr(iron_oc.fs.unit.solid_inlet, "flow_mass")
        assert hasattr(iron_oc.fs.unit.solid_inlet, "mass_frac")
        assert hasattr(iron_oc.fs.unit.solid_inlet, "temperature")

        assert hasattr(iron_oc.fs.unit, "gas_outlet")
        assert len(iron_oc.fs.unit.gas_outlet.vars) == 4
        assert hasattr(iron_oc.fs.unit.gas_outlet, "flow_mol")
        assert hasattr(iron_oc.fs.unit.gas_outlet, "mole_frac")
        assert hasattr(iron_oc.fs.unit.gas_outlet, "temperature")
        assert hasattr(iron_oc.fs.unit.gas_outlet, "pressure")

        assert hasattr(iron_oc.fs.unit, "solid_outlet")
        assert len(iron_oc.fs.unit.solid_outlet.vars) == 3
        assert hasattr(iron_oc.fs.unit.solid_outlet, "flow_mass")
        assert hasattr(iron_oc.fs.unit.solid_outlet, "mass_frac")
        assert hasattr(iron_oc.fs.unit.solid_outlet, "temperature")

        assert number_variables(iron_oc) == 1412
        assert number_total_constraints(iron_oc) == 1371
        assert number_unused_variables(iron_oc) == 15
#        print(unused_variables_set(iron_oc))

    def test_dof(self, iron_oc):
        # Fix geometry variables
        iron_oc.fs.unit.number_orifice.fix(2500)  # [-]
        iron_oc.fs.unit.bed_diameter.fix(6.5)  # m
        iron_oc.fs.unit.bed_height.fix(5)  # m

        # Fix inlet port variables for gas and solid
        iron_oc.fs.unit.gas_inlet.flow_mol[0].fix(272.81)  # mol/s
        iron_oc.fs.unit.gas_inlet.temperature[0].fix(1186)  # K
        iron_oc.fs.unit.gas_inlet.pressure[0].fix(1.86)  # bar
        iron_oc.fs.unit.gas_inlet.mole_frac[0, "CO2"].fix(0.4772)
        iron_oc.fs.unit.gas_inlet.mole_frac[0, "H2O"].fix(0.0646)
        iron_oc.fs.unit.gas_inlet.mole_frac[0, "CH4"].fix(0.4582)

        iron_oc.fs.unit.solid_inlet.flow_mass[0].fix(1422)  # kg/s
        iron_oc.fs.unit.solid_inlet.temperature[0].fix(1186)  # K
        iron_oc.fs.unit.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.45)
        iron_oc.fs.unit.solid_inlet.mass_frac[0, "Fe3O4"].fix(1e-9)
        iron_oc.fs.unit.solid_inlet.mass_frac[0, "Al2O3"].fix(0.55)

        assert degrees_of_freedom(iron_oc) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, iron_oc):
        initialization_tester(
                iron_oc,
                optarg={'tol': 1e-6},
                gas_phase_state_args={"flow_mol": 272.81,
                                      "temperature": 1186,
                                      "pressure": 1.86},
                solid_phase_state_args={"flow_mass": 1422,
                                        "temperature": 1186})

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, iron_oc):
        results = solver.solve(iron_oc)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.ui
    def test_report(self, iron_oc):
        iron_oc.fs.unit.report()
        
## -----------------------------------------------------------------------------
#def test_build():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#
#    # Set up thermo props and reaction props
#    m.fs.gas_properties = GasPhaseThermoParameterBlock()
#    m.fs.solid_properties = SolidPhaseThermoParameterBlock()
#    m.fs.hetero_reactions = HeteroReactionParameterBlock(
#            default={"solid_property_package": m.fs.solid_properties,
#                     "gas_property_package": m.fs.gas_properties})
#
#    m.fs.BFB = BubblingFluidizedBed(
#            default={
#                    "flow_type": "co_current",
#                    "finite_elements": 5,
#                    "transformation_method": "dae.collocation",
#                    "gas_phase_config":
#                    {"property_package": m.fs.gas_properties},
#                    "solid_phase_config":
#                    {"property_package": m.fs.solid_properties,
#                     "reaction_package": m.fs.hetero_reactions
#                     }})
#
#    assert m.fs.BFB.config.flow_type == "co_current"
#    assert m.fs.BFB.config.finite_elements == 5
#    assert m.fs.BFB.config.transformation_method == "dae.collocation"
#
#    # There should be 14 DOFs in this model:
#    # Geometry - 3 (bed length, bed height, and number of orifices)
#    # Gas feed - 6 (inlet flow, temperature, pressure and mole fractions(3))
#    # Solid feed - 5 (inlet flow, temperature, and mass fractions (3))
#
#    assert degrees_of_freedom(m) == 14
#
#
#@pytest.mark.skipif(solver is None, reason="Solver not available")
#def test_initialize():
#    m = ConcreteModel()
#    m.fs = FlowsheetBlock(default={"dynamic": False})
#
#    # Set up thermo props and reaction props
#    m.fs.gas_properties = GasPhaseThermoParameterBlock()
#    m.fs.solid_properties = SolidPhaseThermoParameterBlock()
#    m.fs.hetero_reactions = HeteroReactionParameterBlock(
#            default={"solid_property_package": m.fs.solid_properties,
#                     "gas_property_package": m.fs.gas_properties})
#
#    m.fs.BFB = BubblingFluidizedBed(
#            default={
#                    "flow_type": "co_current",
#                    "finite_elements": 5,
#                    "transformation_method": "dae.collocation",
#                    "gas_phase_config":
#                    {"property_package": m.fs.gas_properties},
#                    "solid_phase_config":
#                    {"property_package": m.fs.solid_properties,
#                     "reaction_package": m.fs.hetero_reactions
#                     }})
#    # Fix geometry variables
#    m.fs.BFB.number_orifice.fix(2500)  # [-]
#    m.fs.BFB.bed_diameter.fix(6.5)  # m
#    m.fs.BFB.bed_height.fix(5)  # m
#
#    # Fix inlet port variables for gas and solid
#    m.fs.BFB.gas_inlet.flow_mol[0].fix(272.81)  # mol/s
#    m.fs.BFB.gas_inlet.temperature[0].fix(373)  # K
#    m.fs.BFB.gas_inlet.pressure[0].fix(1.86)  # bar
#    m.fs.BFB.gas_inlet.mole_frac[0, "CO2"].fix(0.4772)
#    m.fs.BFB.gas_inlet.mole_frac[0, "H2O"].fix(0.0646)
#    m.fs.BFB.gas_inlet.mole_frac[0, "CH4"].fix(0.4582)
#
#    m.fs.BFB.solid_inlet.flow_mass[0].fix(1422)  # kg/s
#    m.fs.BFB.solid_inlet.temperature[0].fix(1186)  # K
#    m.fs.BFB.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.45)
#    m.fs.BFB.solid_inlet.mass_frac[0, "Fe3O4"].fix(1e-9)
#    m.fs.BFB.solid_inlet.mass_frac[0, "Al2O3"].fix(0.55)
#
#    assert degrees_of_freedom(m.fs) == 0
#
#    # State arguments for initializing property state blocks
#    # Bubble and gas_emulsion temperatures are initialized at solid
#    # temperature because thermal mass of solid >> thermal mass of gas
#    blk = m.fs.BFB
#    gas_phase_state_args = {
#            'flow_mol': blk.gas_inlet.flow_mol[0].value,
#            'temperature': blk.solid_inlet.temperature[0].value,
#            'pressure': blk.gas_inlet.pressure[0].value,
#            'mole_frac': {
#                'CH4': blk.gas_inlet.mole_frac[0, 'CH4'].value,
#                'CO2': blk.gas_inlet.mole_frac[0, 'CO2'].value,
#                'H2O': blk.gas_inlet.mole_frac[0, 'H2O'].value}}
#    solid_phase_state_args = {
#            'flow_mass': blk.solid_inlet.flow_mass[0].value,
#            'temperature': blk.solid_inlet.temperature[0].value,
#            'mass_frac': {
#                    'Fe2O3': blk.solid_inlet.mass_frac[0, 'Fe2O3'].value,
#                    'Fe3O4': blk.solid_inlet.mass_frac[0, 'Fe3O4'].value,
#                    'Al2O3': blk.solid_inlet.mass_frac[0, 'Al2O3'].value}}
#
#    m.fs.BFB.initialize(outlvl=5,
#                        gas_phase_state_args=gas_phase_state_args,
#                        solid_phase_state_args=solid_phase_state_args)
#
#    assert degrees_of_freedom(m.fs) == 0
