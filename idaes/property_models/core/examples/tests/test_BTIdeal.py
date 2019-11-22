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
Tests for Ideal + Ideal Liquid (i.e. no activity coefficient) state block;
only tests for construction as parameters need to be provided or estimated
from VLE data to compute the activity coefficients.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import ConcreteModel, Set, SolverFactory, TerminationCondition, \
    SolverStatus, value

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import get_default_solver

from idaes.property_models.core.state_definitions import FPTx
import idaes.property_models.core.eos.ideal as ideal
from idaes.property_models.core.phase_equil import smooth_VLE
from idaes.property_models.core.generic.bubble_dew import (bubble_temp_ideal,
                                                           dew_temp_ideal,
                                                           bubble_press_ideal,
                                                           dew_press_ideal)

import idaes.property_models.core.pure.Perrys as Perrys
import idaes.property_models.core.pure.RPP as RPP

from idaes.property_models.core.examples.BT_ideal \
    import BTIdealParameterBlock


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = BTIdealParameterBlock()

        return model

    def test_config(self, model):
        assert len(model.params.config) == 15

        assert model.params.config.state_definition == FPTx

        assert model.params.config.state_bounds == {
                "flow_mol": (0, 1000),
                "temperature": (273.15, 450),
                "pressure": (5e4, 1e6)}

        assert model.params.config.equation_of_state == {
                "Vap": ideal,
                "Liq": ideal}

        assert model.params.config.phase_equilibrium_formulation == smooth_VLE

        assert model.params.config.bubble_temperature == bubble_temp_ideal
        assert model.params.config.dew_temperature == dew_temp_ideal
        assert model.params.config.bubble_pressure == bubble_press_ideal
        assert model.params.config.dew_pressure == dew_press_ideal

        assert model.params.config.dens_mol_comp_liq == Perrys
        assert model.params.config.enth_mol_comp_liq == Perrys
        assert model.params.config.enth_mol_comp_ig == RPP
        assert model.params.config.entr_mol_comp_liq == Perrys
        assert model.params.config.entr_mol_comp_ig == RPP
        assert model.params.config.pressure_sat_comp == RPP

    def test_build(self, model):
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ['benzene',
                         'toluene']

        for i in model.params.phase_comp.values():
            assert i == model.params.component_list

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2

        assert model.params.phase_equilibrium_list == {
                1: ["benzene", ("Vap", "Liq")],
                2: ["toluene", ("Vap", "Liq")]}

#        # Thermodynamic reference state
#        self.pressure_ref = Param(mutable=True,
#                                  default=101325,
#                                  doc='Reference pressure [Pa]')
#        self.temperature_ref = Param(mutable=True,
#                                     default=298.15,
#                                     doc='Reference temperature [K]')
#
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        pressure_crit_data = {'benzene': 48.9e5,
#                              'toluene': 41e5}
#
#        self.pressure_crit = Param(
#            self.component_list,
#            within=NonNegativeReals,
#            mutable=False,
#            initialize=pressure_crit_data,
#            doc='Critical pressure [Pa]')
#
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        temperature_crit_data = {'benzene': 562.2,
#                                 'toluene': 591.8}
#
#        self.temperature_crit = Param(
#            self.component_list,
#            within=NonNegativeReals,
#            mutable=False,
#            initialize=temperature_crit_data,
#            doc='Critical temperature [K]')
#
#        # Gas Constant
#        self.gas_const = Param(within=NonNegativeReals,
#                               mutable=False,
#                               default=8.314,
#                               doc='Gas constant [J/mol.K]')
#
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        mw_comp_data = {'benzene': 78.1136E-3,
#                        'toluene': 92.1405E-3}
#
#        self.mw_comp = Param(self.component_list,
#                             mutable=False,
#                             initialize=mw_comp_data,
#                             doc="Molecular weight [kg/mol]")
#
#        # Constants for ideal gas specific enthalpy
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        cp_ig_data = {('benzene', 'A'): -3.392E1,
#                      ('benzene', 'B'): 4.739E-1,
#                      ('benzene', 'C'): -3.017E-4,
#                      ('benzene', 'D'): 7.130E-8,
#                      ('toluene', 'A'): -2.435E1,
#                      ('toluene', 'B'): 5.125E-1,
#                      ('toluene', 'C'): -2.765E-4,
#                      ('toluene', 'D'): 4.911E-8}
#
#        self.cp_ig = Param(self.component_list,
#                           ['A', 'B', 'C', 'D'],
#                           mutable=False,
#                           initialize=cp_ig_data,
#                           doc="Parameters for ideal gas heat capacity")
#
#        # Constants for liquid phase specific enthalpy
#        # Source: Perry's Chemical Engineers Handbook 7th Ed.
#        # Units converted to J/mol.K
#        cp_liq_data = {('benzene', '1'): 1.29E2,
#                       ('benzene', '2'): -1.7E-1,
#                       ('benzene', '3'): 6.48E-4,
#                       ('benzene', '4'): 0,
#                       ('benzene', '5'): 0,
#                       ('toluene', '1'): 1.40E2,
#                       ('toluene', '2'): -1.52E-1,
#                       ('toluene', '3'): 6.95E-4,
#                       ('toluene', '4'): 0,
#                       ('toluene', '5'): 0}
#
#        self.cp_liq = Param(self.component_list,
#                            ['1', '2', '3', '4', '5'],
#                            mutable=False,
#                            initialize=cp_liq_data,
#                            doc="Parameters for liquid cp")
#
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        pressure_sat_coeff_data = {('benzene', 'A'): -6.98273,
#                                   ('benzene', 'B'): 1.33213,
#                                   ('benzene', 'C'): -2.62863,
#                                   ('benzene', 'D'): -3.33399,
#                                   ('toluene', 'A'): -7.28607,
#                                   ('toluene', 'B'): 1.38091,
#                                   ('toluene', 'C'): -2.83433,
#                                   ('toluene', 'D'): -2.79168}
#
#        self.pressure_sat_coeff = Param(
#            self.component_list,
#            ['A', 'B', 'C', 'D'],
#            mutable=False,
#            initialize=pressure_sat_coeff_data,
#            doc="parameters to compute Cp_comp")
#
#        # Source: "Perry's Chemical Engineers Handbook by Robert H. Perry"
#        # 7th Edition, pg. 2-98
#        # Units converted to mol/m^3
#        dens_mol_liq_coeff_data = {('benzene', '1'): 1.0162*1e3,
#                                   ('benzene', '2'): 0.2655,
#                                   ('benzene', '3'): 562.16,
#                                   ('benzene', '4'): 0.28212,
#                                   ('toluene', '1'): 0.8488*1e3,
#                                   ('toluene', '2'): 0.26655,
#                                   ('toluene', '3'): 591.8,
#                                   ('toluene', '4'): 0.2878}
#        self.dens_mol_liq_coeff = Param(
#            self.component_list,
#            ['1', '2', '3', '4'],
#            mutable=False,
#            initialize=dens_mol_liq_coeff_data,
#            doc="Coefficients for calculating liquid molar densities")
#
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        dh_vap_data = {'benzene': 3.377e4,
#                       'toluene': 3.8262e4}
#
#        self.dh_vap_ref = Param(
#                self.component_list,
#                mutable=False,
#                initialize=dh_vap_data,
#                doc="Molar heat of vaporization @ Tref [J/mol]")
#
#        # Source: The Properties of Gases and Liquids (1987)
#        # 4th edition, Chemical Engineering Series - Robert C. Reid
#        ds_vap_data = {'benzene': 3.377e4/298.15,
#                       'toluene': 3.8262e4/298.15}
#
#        self.ds_vap_ref = Param(
#                self.component_list,
#                mutable=False,
#                initialize=ds_vap_data,
#                doc="Molar entropy of vaporization @ Tref [J/mol.K]")
