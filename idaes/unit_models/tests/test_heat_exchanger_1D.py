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
Tests for Heat Exchanger unit model.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import ConcreteModel, SolverFactory

from idaes.core import FlowsheetBlock
from idaes.unit_models.heat_exchanger_1D import HeatExchanger1D as HX1D

from idaes.property_models.BFW_properties import PhysicalParameterBlock
from idaes.ui.report import degrees_of_freedom

# -----------------------------------------------------------------------------
# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


# -----------------------------------------------------------------------------
def test_build():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterBlock()
    m.fs.HX1D = HX1D(default={"shell_property_package": m.fs.properties,
                              "tube_property_package": m.fs.properties})

    assert len(m.fs.HX1D.config) == 29

    # Check for inlets/outlets construction
    assert hasattr(m.fs.HX1D, "shell_inlet")
    assert hasattr(m.fs.HX1D, "shell_outlet")
    assert hasattr(m.fs.HX1D, "tube_inlet")
    assert hasattr(m.fs.HX1D, "tube_outlet")


# # -----------------------------------------------------------------------------
# # Test class methods
# def test_default_config_block(m):

#     assert m.fs.hx.config.dynamic == 'use_parent_value'
#     assert not m.fs.hx.config.include_holdup
#     assert m.fs.hx.config.material_balance_type == 'component_phase'
#     assert m.fs.hx.config.energy_balance_type == 'enthalpy_total'
#     assert m.fs.hx.config.momentum_balance_type == 'pressure'
#     assert not m.fs.hx.config.has_rate_reactions
#     assert not m.fs.hx.config.has_equilibrium_reactions
#     assert not m.fs.hx.config.has_mass_transfer
#     assert m.fs.hx.config.has_heat_transfer
#     assert not m.fs.hx.config.has_work_transfer
#     assert not m.fs.hx.config.has_pressure_change
#     assert m.fs.hx.config.velocity_type == 'none'
#
#     assert not hasattr(m.fs.hx.config, "has_phase_equilibrium")
#
#     assert m.fs.hx.config.shell_has_phase_equilibrium
#     m.fs.hx.config.shell_has_phase_equilibrium = False
#     with pytest.raises(ValueError):
#         m.fs.hx.config.shell_has_phase_equilibrium = 'foo'
#         m.fs.hx.config.shell_has_phase_equilibrium = 10
#
#     assert m.fs.hx.config.tube_has_phase_equilibrium
#     m.fs.hx.config.tube_has_phase_equilibrium = False
#     with pytest.raises(ValueError):
#         m.fs.hx.config.tube_has_phase_equilibrium = 'foo'
#         m.fs.hx.config.tube_has_phase_equilibrium = 10
#
#     assert m.fs.hx.config.shell_property_package == m.fs.props
#     assert m.fs.hx.config.shell_property_package_args == {}
#     assert m.fs.hx.config.shell_inlet_list is None
#     assert m.fs.hx.config.shell_num_inlets is None
#     assert m.fs.hx.config.shell_outlet_list is None
#     assert m.fs.hx.config.shell_num_outlets is None
#
#     assert m.fs.hx.config.tube_property_package == m.fs.props
#     assert m.fs.hx.config.tube_property_package_args == {}
#     assert m.fs.hx.config.tube_inlet_list is None
#     assert m.fs.hx.config.tube_num_inlets is None
#     assert m.fs.hx.config.tube_outlet_list is None
#     assert m.fs.hx.config.tube_num_outlets is None
#
#     assert m.fs.hx.config.shell_discretization_method == 'OCLR'
#     assert m.fs.hx.config.tube_discretization_method == 'OCLR'
#     assert m.fs.hx.config.finite_elements == 20
#     assert m.fs.hx.config.collocation_points == 5
#
#     assert not m.fs.hx.config.shell_has_mass_diffusion
#     m.fs.hx.config.shell_has_mass_diffusion = True
#     with pytest.raises(ValueError):
#         m.fs.hx.config.shell_has_mass_diffusion = 'foo'
#     with pytest.raises(ValueError):
#         m.fs.hx.config.shell_has_mass_diffusion = 10
#
#     assert not m.fs.hx.config.shell_has_energy_diffusion
#     m.fs.hx.config.shell_has_energy_diffusion = True
#     with pytest.raises(ValueError):
#         m.fs.hx.config.shell_has_energy_diffusion = 'foo'
#     with pytest.raises(ValueError):
#         m.fs.hx.config.shell_has_energy_diffusion = 10
#
#     assert not m.fs.hx.config.tube_has_mass_diffusion
#     m.fs.hx.config.tube_has_mass_diffusion = True
#     with pytest.raises(ValueError):
#         m.fs.hx.config.tube_has_mass_diffusion = 'foo'
#     with pytest.raises(ValueError):
#         m.fs.hx.config.tube_has_mass_diffusion = 10
#
#     assert not m.fs.hx.config.tube_has_energy_diffusion
#     m.fs.hx.config.tube_has_energy_diffusion = True
#     with pytest.raises(ValueError):
#         m.fs.hx.config.tube_has_energy_diffusion = 'foo'
#     with pytest.raises(ValueError):
#         m.fs.hx.config.tube_has_energy_diffusion = 10
#
#     assert m.fs.hx.config.flow_type == 'co_current'
#
#     assert m.fs.hx.config.has_wall_conduction == 'none'
