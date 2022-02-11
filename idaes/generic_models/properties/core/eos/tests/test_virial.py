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
Author: Paul Akula
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import VaporPhase, Component
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util import get_solver
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock, StateIndex)
from idaes.generic_models.properties.core.pure.RPP4 import (cp_mol_ig_comp,
                                                            enth_mol_ig_comp)
from idaes.generic_models.properties.core.eos.virial import (
    Virial, entr_mol_ig_comp_G_H_ref)
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent import PressureSatSolvent


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# ------------------------------------------------------------------------------
# Test Case For Pure Species
# components: CO2 and H2O
# temperature : 317 K
# pressure : 107650 Pa
# Expected Results
S_ref_CO2 = 2.95153
S_ref_H2O = -44.37364
Z_CO2 = 0.99562
Z_H2O = 0.97633
V_CO2 = 0.02438
V_H2O = 0.02390
log_fug_coeff_CO2 = -0.00438
log_fug_coeff_H2O = -0.02367
fug_coeff_CO2 = 0.99563
fug_coeff_H2O = 0.97660
fug_CO2 = 107179.09671
fug_H2O = 105131.49559
log_fug_coeff_sat_H2O = -0.00199
fug_coeff_sat_H2O = 0.99801
fug_sat_H2O = 9025.032021116094
H_CO2 = -392831.48003
H_H2O = -241433.33063
S_CO2 = 4.66236
S_H2O = -43.39847
U_CO2 = -395455.60987
U_H2O = -244006.61988
G_CO2 = -394309.44909
G_H2O = -227676.01610

# -----------------------------------------------------------------------------
# Test Case For Gas Mixture
# Limiting condition as Fluegas mixture (CO2,H2O,N2,O2) tends to CO2 gas
# @ Temperature : 317 K, @pressure : 107650 Pa and
# Composition: CO2= 0.99997, H2O= 0.00001, N2= 0.00001, O2= 0.00001
# Expected Results
# flue gas mixture properties becomes approximately equal to pure CO2 properties

# ------------------------------------------------------------------------------
# Test case for Partial molar properties
# Summability relationships
# Expected results
# solution property = sum over all species(mole fraction * partial molar
# property)

#------------------------------------------------------------------------------
# Configuration dictionary for flue gas (gas mixtue)

flue_gas = {
    # Specifying components
    "components": {
        'CO2': {"type": Component,
                "cp_mol_ig_comp": cp_mol_ig_comp,
                "enth_mol_ig_comp": enth_mol_ig_comp,
                "entr_mol_ig_comp": entr_mol_ig_comp_G_H_ref,
                "parameter_data": {
                    "mw": (0.04401, pyunits.kg / pyunits.mol),
                    "omega": 0.2236,
                    "pressure_crit": (73.773, pyunits.bar),
                    "temperature_crit": (304.128, pyunits.K),
                    "volume_crit": (94.117, pyunits.cm**3 / pyunits.mol),
                    "enth_mol_form_vap_comp_ref": (-393500,
                                                   pyunits.J / pyunits.mol),
                    "gibbs_mol_form_vap_comp_ref": (-394380,
                                                    pyunits.J / pyunits.mol),
                    "cp_mol_ig_comp_coeff": {
                        'A': 1.9795e1,
                        'B': 7.3437e-2,
                        'C': -5.6019e-5,
                        'D': 1.7153e-8}
                }},
        'H2O': {"type": Component,
                "cp_mol_ig_comp": cp_mol_ig_comp,
                "enth_mol_ig_comp": enth_mol_ig_comp,
                "entr_mol_ig_comp": entr_mol_ig_comp_G_H_ref,
                "pressure_sat_comp": PressureSatSolvent,
                "parameter_data": {
                    "mw": (0.01802, pyunits.kg / pyunits.mol),
                    "omega": 0.3443,
                    "pressure_crit": (220.64, pyunits.bar),
                    "volume_crit": (55.947, pyunits.cm**3 / pyunits.mol),
                    "temperature_crit": (647.096, pyunits.K),
                    "enth_mol_form_vap_comp_ref": (-241820,
                                                   pyunits.J / pyunits.mol),
                    "gibbs_mol_form_vap_comp_ref": (-228590,
                                                    pyunits.J / pyunits.mol),
                    "cp_mol_ig_comp_coeff": {
                        'A': 32.22,
                        'B': 1.923E-3,
                        'C': 10.548E-6,
                        'D': -3.594E-9},
                    "pressure_sat_comp_coeff": {
                        '1': 73.649,
                        '2': -7258.2,
                        '3': -7.3037,
                        '4': 4.1653e-6},
                }},
        'N2': {"type": Component,
               "cp_mol_ig_comp": cp_mol_ig_comp,
               "enth_mol_ig_comp": enth_mol_ig_comp,
               "entr_mol_ig_comp": entr_mol_ig_comp_G_H_ref,
               "parameter_data": {
                   "mw": (0.02801, pyunits.kg / pyunits.mol),
                   "omega": 0.0372,
                   "pressure_crit": (33.958, pyunits.bar),
                   "temperature_crit": (126.192, pyunits.K),
                   "volume_crit": (89.416, pyunits.cm**3 / pyunits.mol),
                   "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
                   "gibbs_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
                   "cp_mol_ig_comp_coeff": {
                       'A': 31.128,
                       'B': -13.556E-3,
                       'C': 26.777E-6,
                       'D': -11.673E-9}
               }},
        'O2': {"type": Component,
               "cp_mol_ig_comp": cp_mol_ig_comp,
               "enth_mol_ig_comp": enth_mol_ig_comp,
               "entr_mol_ig_comp": entr_mol_ig_comp_G_H_ref,
               "parameter_data": {
                   "mw": (0.032, pyunits.kg / pyunits.mol),
                   "omega": 0.0221,
                   "pressure_crit": (50.464, pyunits.bar),
                   "volume_crit": (74.948, pyunits.cm**3 / pyunits.mol),
                   "temperature_crit": (154.599, pyunits.K),
                   "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
                   "gibbs_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),
                   "cp_mol_ig_comp_coeff": {
                       'A': 28.087,
                       'B': -0.004E-3,
                       'C': 17.447E-6,
                       'D': -10.644E-9}
               }}},

    # Specifying phases
    "phases": {'Vap': {"type": VaporPhase,
                       "equation_of_state": Virial,
                       "equation_of_state_options": {
                           "use_pseudocritical_rules": False}
                       }
               },
    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 1, 1000, pyunits.mol / pyunits.s),
                     "temperature": (273.15, 298.15, 450, pyunits.K),
                     "pressure": (5e4, 101325, 1e6, pyunits.Pa)},
    "state_components": StateIndex.apparent,
    "include_enthalpy_of_formation": True,
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K)}


class TestVirialEOS(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=flue_gas)

        model.props = model.params.build_state_block(
            [1],
            default={"defined_state": True})

        model.props[1].calculate_scaling_factors()
        # Fix state
        model.props[1].flow_mol.fix(21)
        model.props[1].temperature.fix(317)
        model.props[1].pressure.fix(107650)
        model.props[1].mole_frac_comp["CO2"].fix(0.11453)
        model.props[1].mole_frac_comp["H2O"].fix(0.08526)
        model.props[1].mole_frac_comp["N2"].fix(0.73821)
        model.props[1].mole_frac_comp["O2"].fix(0.06200)

        return model

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):

        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_pure_properties(self, model):
        solver.solve(model)

        assert pytest.approx(S_ref_CO2, abs=1e-5) == value(
            model.params.CO2.entr_mol_form_vap_comp_ref)

        assert pytest.approx(S_ref_H2O, abs=1e-5) == value(
            model.params.H2O.entr_mol_form_vap_comp_ref)

        assert pytest.approx(Z_CO2, abs=1e-5) == value(
            Virial.compress_fact_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(Z_H2O, abs=1e-5) == value(
            Virial.compress_fact_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(V_CO2, abs=1e-5) == value(
            Virial.vol_mol_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(V_H2O, abs=1e-5) == value(
            Virial.vol_mol_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(log_fug_coeff_CO2, abs=1e-5) == value(
            Virial.log_fug_coeff_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(log_fug_coeff_H2O, abs=1e-5) == value(
            Virial.log_fug_coeff_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(fug_coeff_CO2, abs=1e-5) == value(
            Virial.fug_coeff_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(fug_coeff_H2O, abs=1e-5) == value(
            Virial.fug_coeff_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(fug_CO2, abs=1e-5) == value(
            Virial.fug_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(fug_H2O, abs=1e-5) == value(
            Virial.fug_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(log_fug_coeff_sat_H2O, abs=1e-5) == value(
            Virial.log_fug_coeff_vap_sat_comp(model.props[1], 'H2O'))

        assert pytest.approx(fug_coeff_sat_H2O, abs=1e-5) == value(
            Virial.fug_coeff_vap_sat_comp(model.props[1], 'H2O'))

        assert pytest.approx(fug_sat_H2O, abs=1e-5) == value(
            Virial.fug_vap_sat_comp(model.props[1], 'H2O'))

        assert pytest.approx(H_CO2, abs=1e-5) == value(
            Virial.enth_mol_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(H_H2O, abs=1e-5) == value(
            Virial.enth_mol_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(S_CO2, abs=1e-5) == value(
            Virial.entr_mol_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(S_H2O, abs=1e-5) == value(
            Virial.entr_mol_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(U_CO2, abs=1e-5) == value(
            Virial.energy_internal_mol_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(U_H2O, abs=1e-5) == value(
            Virial.energy_internal_mol_vap_comp_pure(model.props[1], 'H2O'))

        assert pytest.approx(G_CO2, abs=1e-5) == value(
            Virial.gibbs_mol_vap_comp_pure(model.props[1], 'CO2'))

        assert pytest.approx(G_H2O, abs=1e-5) == value(
            Virial.gibbs_mol_vap_comp_pure(model.props[1], 'H2O'))

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_summability_relationships(self, model):

        solver.solve(model)

        G_total = value(sum(model.props[1].mole_frac_comp[i] *
                            model.props[1].gibbs_mol_phase_comp['Vap', i]
                            for i in model.props.component_list))
        G = value(model.props[1].gibbs_mol_phase['Vap'])

        assert G / G_total == pytest.approx(1.0000, rel=1e-5)

        H_total = value(sum(model.props[1].mole_frac_comp[i] *
                            model.props[1].enth_mol_phase_comp['Vap', i]
                            for i in model.props.component_list))
        H = value(model.props[1].enth_mol_phase['Vap'])

        assert H / H_total == pytest.approx(1.0000, rel=1e-5)

        S_total = value(sum(model.props[1].mole_frac_comp[i] *
                            model.props[1].entr_mol_phase_comp['Vap', i]
                            for i in model.props.component_list))
        S = value(model.props[1].entr_mol_phase['Vap'])

        assert S / S_total == pytest.approx(1.0000, rel=1e-5)

        U_total = value(sum(model.props[1].mole_frac_comp[i] *
                            model.props[1].energy_internal_mol_phase_comp['Vap', i]
                            for i in model.props.component_list))
        U = value(model.props[1].energy_internal_mol_phase['Vap'])

        assert U / U_total == pytest.approx(1.0000, rel=1e-5)

        V_total = value(sum(model.props[1].mole_frac_comp[i] *
                            model.props[1].vol_mol_phase_comp['Vap', i]
                            for i in model.props.component_list))
        V = value(model.props[1].vol_mol_phase['Vap'])

        assert V / V_total == pytest.approx(1.0000, rel=1e-5)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_mixture_properties(self, model):
        # Fix state
        model.props[1].flow_mol.fix(21)
        model.props[1].temperature.fix(317)
        model.props[1].pressure.fix(107650)
        # Fix mixture composition to approximatley  pure CO2 gas
        model.props[1].mole_frac_comp["CO2"].fix(0.99997)
        model.props[1].mole_frac_comp["H2O"].fix(0.00001)
        model.props[1].mole_frac_comp["N2"].fix(0.00001)
        model.props[1].mole_frac_comp["O2"].fix(0.00001)
        solver.solve(model)

        # pure CO2 properties
        fug_CO2_from_pure_pp = value(Virial.fug_vap_comp_pure(
                                     model.props[1], 'CO2'))
        Z_CO2_from_pure_pp = value(Virial.compress_fact_vap_comp_pure(
                                   model.props[1], 'CO2'))
        H_CO2_from_pure_pp = value(Virial.enth_mol_vap_comp_pure(
                                   model.props[1], 'CO2'))
        S_CO2_from_pure_pp = value(Virial.entr_mol_vap_comp_pure(
                                   model.props[1], 'CO2'))
        U_CO2_from_pure_pp = value(Virial.energy_internal_mol_vap_comp_pure(
                                   model.props[1], 'CO2'))
        G_CO2_from_pure_pp = value(Virial.gibbs_mol_vap_comp_pure(
                                   model.props[1], 'CO2'))
        V_CO2_from_pure_pp = value(Virial.vol_mol_vap_comp_pure(
                                   model.props[1], 'CO2'))

        # Mixture properties with composition approximately equal to pure CO2
        fug_CO2_from_mixture_pp = value(Virial.fug_phase_comp(
                                        model.props[1], 'Vap', 'CO2'))
        Z_CO2_from_mixture_pp = value(
            model.props[1].compress_fact_phase['Vap'])
        H_CO2_from_mixture_pp = value(model.props[1].enth_mol_phase['Vap'])
        S_CO2_from_mixture_pp = value(model.props[1].entr_mol_phase['Vap'])
        G_CO2_from_mixture_pp = value(model.props[1].gibbs_mol_phase['Vap'])
        V_CO2_from_mixture_pp = value(model.props[1].vol_mol_phase['Vap'])
        U_CO2_from_mixture_pp = \
            value(model.props[1].energy_internal_mol_phase['Vap'])

        assert U_CO2_from_pure_pp / U_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)

        assert H_CO2_from_pure_pp / H_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)

        assert S_CO2_from_pure_pp / S_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)

        assert G_CO2_from_pure_pp / G_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)

        assert V_CO2_from_pure_pp / V_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)

        assert fug_CO2_from_pure_pp / fug_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)

        assert Z_CO2_from_pure_pp / Z_CO2_from_mixture_pp ==\
            pytest.approx(1.000, abs=1e-3)
