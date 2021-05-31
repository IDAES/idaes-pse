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
from pyomo.common.unittest import assertStructuredAlmostEqual

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

from idaes.generic_models.properties.core.examples.CO2_H2O_MEA_Ideal_VLE_Final import (
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
        assert len(model.params.component_list) == 3
        for i in model.params.component_list:
            assert i in ['H2O',
                          'CO2',
                          'MEA']
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 5
        for i in model.params._phase_component_set:
            assert i in [("Liq", "H2O"), ("Vap", "H2O"),
                          ("Liq", "CO2"), ("Vap", "CO2"),
                          ("Liq", "MEA")]

        assert model.params.config.state_definition == FTPx
        
        assertStructuredAlmostEqual(
            model.params.config.state_bounds,
            { "flow_mol": (0, 92.0817, 100, pyunits.mol/pyunits.s),
              "temperature": (350, 380, 500, pyunits.K),
              "pressure": (5e4, 183430, 1e7, pyunits.Pa)},
             item_callback=lambda x: value(x) * (
                pyunits.get_units(x) or pyunits.dimensionless)._get_pint_unit()
        )

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): smooth_VLE}

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1","PE2"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"H2O": ("Vap", "Liq")},
            "PE2": {"CO2": ("Vap", "Liq")}}

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
        assert value(model.props[1].flow_mol) == 92.0817
        assert model.props[1].flow_mol.ub == 100
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 183430
        assert model.props[1].pressure.ub == 1E7
        assert model.props[1].pressure.lb == 5E4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 380
        assert model.props[1].temperature.ub == 500
        assert model.props[1].temperature.lb == 350

        # Check supporting variables
        assert isinstance(model.props[1].flow_mol_phase, Var)
        assert len(model.props[1].flow_mol_phase) == 2

        assert isinstance(model.props[1].mole_frac_phase_comp, Var)
        assert len(model.props[1].mole_frac_phase_comp) == 5

        assert isinstance(model.props[1].phase_frac, Var)
        assert len(model.props[1].phase_frac) == 2

        assert isinstance(model.props[1].total_flow_balance, Constraint)
        assert len(model.props[1].total_flow_balance) == 1

        assert isinstance(model.props[1].component_flow_balances, Constraint)
        assert len(model.props[1].component_flow_balances) == 3

        assert isinstance(model.props[1].sum_mole_frac, Constraint)
        assert len(model.props[1].sum_mole_frac) == 1

        assert not hasattr(model.props[1], "sum_mole_frac_out")

        assert isinstance(model.props[1].phase_fraction_constraint, Constraint)
        assert len(model.props[1].phase_fraction_constraint) == 2

        assert_units_consistent(model)

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
        model.props[1].flow_mol.fix(92.0817)
        model.props[1].temperature.fix(380)
        model.props[1].pressure.fix(183430)
        model.props[1].mole_frac_comp["H2O"].fix(0.8787745)
        model.props[1].mole_frac_comp["CO2"].fix(0.018288)
        model.props[1].mole_frac_comp["MEA"].fix(0.1029375)

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