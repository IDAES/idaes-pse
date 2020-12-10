##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Authors: Anuja Deshpande, Andrew Lee
"""
import pytest
import numpy as np
import sys
import os

from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (MaterialBalanceType,
                        EnergyBalanceType,
                        MaterialFlowBasis,
                        Component)

from idaes.core.flowsheet_model import FlowsheetBlock

from idaes.generic_models.unit_models import Flash

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util.testing import get_default_solver

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.phase_equil import smooth_VLE

from idaes.generic_models.properties.core.examples.CO2_H2O_Ideal_VLE import (
    configuration)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ['H2O',
                         'CO2']
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 3
        for i in model.params._phase_component_set:
            assert i in [("Liq", "H2O"), ("Vap", "H2O"),
                         ("Vap", "CO2")]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
            "flow_mol": (0, 10, 20, pyunits.mol/pyunits.s),
            "temperature": (273.15, 323.15, 1000, pyunits.K),
            "pressure": (5e4, 108900, 1e7, pyunits.Pa),
            "mole_frac_comp": {"H2O":(0,0.5,1),"CO2":(0,0.5,1)}}

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): smooth_VLE}

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 1
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"H2O": ("Vap", "Liq")}}

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.H2O.mw.value == 18.0153E-3
        assert model.params.H2O.pressure_crit.value == 220.64E5
        assert model.params.H2O.temperature_crit.value == 647

        assert model.params.CO2.mw.value == 44.0095E-3
        assert model.params.CO2.pressure_crit.value == 73.825E5
        assert model.params.CO2.temperature_crit.value == 304.23

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 10
        assert model.props[1].flow_mol.ub == 20
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 108900
        assert model.props[1].pressure.ub == 1E7
        assert model.props[1].pressure.lb == 5E4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 323.15
        assert model.props[1].temperature.ub == 1000
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 2
        assert value(model.props[1].mole_frac_comp["H2O"]) == 0.5
        assert value(model.props[1].mole_frac_comp["CO2"]) == 0.5

        # Check supporting variables
        assert isinstance(model.props[1].flow_mol_phase, Var)
        assert len(model.props[1].flow_mol_phase) == 2

        assert isinstance(model.props[1].mole_frac_phase_comp, Var)
        assert len(model.props[1].mole_frac_phase_comp) == 3

        assert isinstance(model.props[1].phase_frac, Var)
        assert len(model.props[1].phase_frac) == 2

        assert isinstance(model.props[1].total_flow_balance, Constraint)
        assert len(model.props[1].total_flow_balance) == 1

        assert isinstance(model.props[1].component_flow_balances, Constraint)
        assert len(model.props[1].component_flow_balances) == 2

        assert isinstance(model.props[1].sum_mole_frac, Constraint)
        assert len(model.props[1].sum_mole_frac) == 1

        assert not hasattr(model.props[1], "sum_mole_frac_out")

        assert isinstance(model.props[1].phase_fraction_constraint, Constraint)
        assert len(model.props[1].phase_fraction_constraint) == 2

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                if (p,j)!=("Liq","CO2"):
                    assert model.props[1].get_material_flow_terms(p, j) == (
                        model.props[1].flow_mol_phase[p] *
                        model.props[1].mole_frac_phase_comp[p, j])

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_enthalpy_flow_terms(p) == (
                model.props[1].flow_mol_phase[p] *
                model.props[1].enth_mol_phase[p])

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                if (p,j)!=("Liq","CO2"):
                    assert model.props[1].get_material_density_terms(p, j) == (
                        model.props[1].dens_mol_phase[p] *
                        model.props[1].mole_frac_phase_comp[p, j])

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_energy_density_terms(p) == (
                model.props[1].dens_mol_phase[p] *
                model.props[1].enth_mol_phase[p])

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert model.props[1].default_material_balance_type() == \
            MaterialBalanceType.componentTotal

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert model.props[1].default_energy_balance_type() == \
            EnergyBalanceType.enthalpyTotal

    @pytest.mark.unit
    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == \
            MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["Total Molar Flowrate",
                         "Total Mole Fraction",
                         "Temperature",
                         "Pressure"]

    @pytest.mark.unit
    def test_dof(self, model):
        # Fix state
        model.props[1].flow_mol.fix(10)
        model.props[1].temperature.fix(323.15)
        model.props[1].pressure.fix(108900)
        model.props[1].mole_frac_comp["H2O"].fix(0.5)
        model.props[1].mole_frac_comp["CO2"].fix(0.5)

        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
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

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp["Liq", "H2O"].value == \
            pytest.approx(1, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp["Vap", "H2O"].value == \
            pytest.approx(0.11501, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp["Vap", "CO2"].value == \
            pytest.approx(0.88499, abs=1e-4)
        assert model.props[1].phase_frac["Vap"].value == \
            pytest.approx(0.56498, abs=1e-4)
            
    @pytest.mark.integration
    def test_temp_swing(self):
        # Create a flash model with the CO2-H2O property package
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = GenericParameterBlock(default=configuration)
        m.fs.flash = Flash(default={"property_package": m.fs.properties})
        
        # Fix inlet stream state variables
        m.fs.flash.inlet.flow_mol.fix(9.89433124673833) # mol/s
        m.fs.flash.inlet.mole_frac_comp[0, 'CO2'].fix(0.13805801934749645)
        m.fs.flash.inlet.mole_frac_comp[0, 'H2O'].fix(0.8619419806525035)
        m.fs.flash.inlet.pressure.fix(183430) # Pa
        m.fs.flash.inlet.temperature.fix(396.79057912844183) # K
        
        # Fix flash and its outlet conditions
        m.fs.flash.deltaP.fix(0)
        m.fs.flash.vap_outlet.temperature.fix(313.15)
        
        #Initialize the flash model
        initial = m.fs.flash.initialize(outlvl=0)
        
        # Create a dictionary of expected solution for flash outlet temperature 
        #sweep
        temp_range = list(np.linspace(313,396))
        
        expected_vapor_frac = [0.14388,0.14445,0.14508,0.14576,0.1465,0.14731,
                               0.14818,0.14913,0.15017,0.15129,0.15251,0.15384,
                               0.15528,0.15685,0.15856,0.16042,0.16245,0.16467,
                               0.16709,0.16974,0.17265,0.17584,0.17935,0.18323,
                               0.18751,0.19226,0.19755,0.20346,0.21008,0.21755,
                               0.22601,0.23565,0.24673,0.25956,0.27456,0.29229,
                               0.31354,0.33942,0.37157,0.4125,0.46628,0.53993,
                               0.64678,0.81547,1.0,1.0,1.0,1.0,1.0,1.0]
        
        expected_heat_duty = [-319572.91269,-318207.35249,-316825.31456,
                              -315425.46452,-314006.35514,-312566.41317,
                              -311103.92414,-309617.01485,-308103.633,
                              -306561.52359,-304988.20138,-303380.91855,
                              -301736.62677,-300051.9324,-298323.04327,
                              -296545.70535,-294715.12686,-292825.88683,
                              -290871.8244,-288845.90382,-286740.04885,
                              -284544.93807,-282249.74995,-279841.84278,
                              -277306.349,-274625.65625,-271778.73623,
                              -268740.26687,-265479.46928,-261958.54546,
                              -258130.54709,-253936.4179,-249300.81058,
                              -244126.04073,-238283.13381,-231598.19109,
                              -223830.94744,-214639.75318,-203521.76902,
                              -189705.16711,-171941.46564,-148070.35054,
                              -114001.22402,-60937.07508,-3219.88715,
                              -2631.66029,-2043.09203,-1454.17752,-864.91582,
                              -275.30623]
        
        outvals = zip(expected_vapor_frac,expected_heat_duty)
        expected_sol = dict(zip(temp_range,outvals))
        
        # Solve the model for a range of flash outlet temperatures
        # Perform flash outlet temperature sweep and test the solution
        for t in temp_range:
            m.fs.flash.vap_outlet.temperature.fix(t)
            res = solver.solve(m)
            assert res.solver.termination_condition == "optimal"
            frac = value(m.fs.flash.vap_outlet.flow_mol[0]) \
                /value(m.fs.flash.inlet.flow_mol[0])
            assert frac == pytest.approx(expected_sol[t][0],abs=1e-4)
            hduty = value(m.fs.flash.heat_duty[0])
            assert hduty == pytest.approx(expected_sol[t][1],abs=1e-4)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()