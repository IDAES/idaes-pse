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
Tests for methanol flowsheet.

"""

import pytest

from io import StringIO
import sys

from idaes.models.flowsheets.methanol_flowsheet_w_recycle import (
    build_model, set_inputs, scale_flowsheet, initialize_flowsheet,
    add_costing, report)

from pyomo.environ import (Constraint,
                           ConcreteModel,
                           value,
                           Var)
from pyomo.network import Arc
from pyomo.environ import TerminationCondition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.scaling import (extreme_jacobian_columns,
                                     extreme_jacobian_rows)

from idaes.models.properties.modular_properties.base.generic_property import \
    GenericParameterBlock
from idaes.models.properties.modular_properties.base.generic_reaction import \
    GenericReactionParameterBlock

from idaes.models.unit_models import (
    Mixer,
    Heater,
    Compressor,
    Turbine,
    StoichiometricReactor,
    Flash,
    Separator as Splitter)


@pytest.fixture(scope='module')
def model():
    m = ConcreteModel()
    m = build_model(m)

    return m


@pytest.mark.unit
def test_build_flowsheet(model):
    assert isinstance(model.fs, FlowsheetBlock)

    assert isinstance(model.fs.thermo_params_VLE, GenericParameterBlock)
    assert isinstance(model.fs.thermo_params_vapor, GenericParameterBlock)
    assert isinstance(model.fs.reaction_params, GenericReactionParameterBlock)

    assert isinstance(model.fs.M101, Mixer)
    assert isinstance(model.fs.M102, Mixer)
    assert isinstance(model.fs.C101, Compressor)
    assert isinstance(model.fs.H101, Heater)
    assert isinstance(model.fs.R101, StoichiometricReactor)
    assert isinstance(model.fs.T101, Turbine)
    assert isinstance(model.fs.H102, Heater)
    assert isinstance(model.fs.F101, Flash)
    assert isinstance(model.fs.S101, Splitter)

    assert isinstance(model.fs.s01, Arc)
    assert isinstance(model.fs.s02, Arc)
    assert isinstance(model.fs.s03, Arc)
    assert isinstance(model.fs.s04, Arc)
    assert isinstance(model.fs.s05, Arc)
    assert isinstance(model.fs.s06, Arc)
    assert isinstance(model.fs.s07, Arc)
    assert isinstance(model.fs.s08, Arc)
    assert isinstance(model.fs.s09, Arc)

    assert degrees_of_freedom(model) == 24


@pytest.mark.unit
def test_set_inputs(model):
    set_inputs(model)

    assert degrees_of_freedom(model) == 0


@pytest.mark.unit
def test_scaling(model):
    scale_flowsheet(model)

    # check that less than 10% of model variables are badly scaled pre-solve
    badly_scaled_vars = extreme_jacobian_columns(model)
    all_var = list(model.component_data_objects(Var, descend_into=True))
    assert len(badly_scaled_vars)/len(all_var) < 0.1

    # check that less than 10% of model constraints are badly scaled pre-solve
    badly_scaled_cons = extreme_jacobian_rows(model)
    all_con = list(model.component_data_objects(Constraint, descend_into=True))
    assert len(badly_scaled_cons)/len(all_con) < 0.1

    # check manual scaling
    for unit in ('M101', 'M102', 'C101', 'H101', 'R101',
                 'T101', 'H102', 'F101', 'S101'):
        block = getattr(model.fs, unit)
        if 'mixer' in unit:
            assert iscale.get_scaling_factor(
                block.inlet_1_state[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.inlet_2_state[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.mixed_state[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.inlet_1_state[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.inlet_2_state[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.mixed_state[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
        if 'splitter' in unit:
            assert iscale.get_scaling_factor(
                block.mixed_state[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.outlet_1_state[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.outlet_2_state[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.mixed_state[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.outlet_1_state[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.outlet_2_state[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
        if hasattr(block, "control_volume"):
            assert iscale.get_scaling_factor(
                block.control_volume.properties_in[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.control_volume.properties_out[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.control_volume.properties_in[0.0].enth_mol_phase) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.control_volume.properties_out[0.0].enth_mol_phase) == \
                pytest.approx(1e2, rel=1e-5)
            if hasattr(block.control_volume, "rate_reaction_extent"):
                assert iscale.get_scaling_factor(
                    block.control_volume.rate_reaction_extent) == \
                    pytest.approx(1e3, rel=1e-5)
            if hasattr(block.control_volume, "heat"):
                assert iscale.get_scaling_factor(
                    block.control_volume.heat) == \
                    pytest.approx(1e-3, rel=1e-5)
            if hasattr(block.control_volume, "work"):
                assert iscale.get_scaling_factor(
                    block.control_volume.work) == \
                    pytest.approx(1e-3, rel=1e-5)
        if hasattr(block, "properties_isentropic"):
            assert iscale.get_scaling_factor(
                block.properties_isentropic[0.0].mole_frac_comp) == \
                pytest.approx(1e2, rel=1e-5)
            assert iscale.get_scaling_factor(
                block.properties_isentropic[0.0].enth_mol_phase) == \
                pytest.approx(1e-3, rel=1e-5)
        if hasattr(block, "properties"):
            assert iscale.get_scaling_factor(
                block.properties[0.0].mole_frac_comp) == \
                    pytest.approx(1e2, rel=1e-5)

    # checking specific constraint scaling factors
    for name in ('M101', 'M102', 'C101', 'H101', 'R101',
                 'T101', 'H102', 'F101', 'S101'):
        unit = getattr(model.fs, name)
        # mixer constraints
        if hasattr(unit, 'material_mixing_equations'):
            for (t, j), c in unit.material_mixing_equations.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'enthalpy_mixing_equations'):
            for t, c in unit.enthalpy_mixing_equations.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'minimum_pressure_constraint'):
            for (t, i), c in unit.minimum_pressure_constraint.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'mixture_pressure'):
            for t, c in unit.mixture_pressure.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'pressure_equality_constraints'):
            for (t, i), c in unit.pressure_equality_constraints.items():
                assert pytest.approx(1e-5, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        # splitter constraints
        if hasattr(unit, 'material_splitting_eqn'):
            for (t, o, j), c in unit.material_splitting_eqn.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'temperature_equality_eqn'):
            for (t, o), c in unit.temperature_equality_eqn.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'molar_enthalpy_equality_eqn'):
            for (t, o), c in unit.molar_enthalpy_equality_eqn.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'molar_enthalpy_splitting_eqn'):
            for (t, o), c in unit.molar_enthalpy_splitting_eqn.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'pressure_equality_eqn'):
            for (t, o), c in unit.pressure_equality_eqn.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'sum_split_frac'):
            for t, c in unit.sum_split_frac.items():
                assert pytest.approx(1e2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        # flash adds same as splitter, plus one more
        if hasattr(unit, 'split_fraction_eq'):
            for (t, o), c in unit.split_fraction_eq.items():
                assert pytest.approx(1e2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        # pressurechanger constraints

        if hasattr(unit, "ratioP_calculation"):
            for t, c in unit.ratioP_calculation.items():
                assert pytest.approx(1e-5, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "fluid_work_calculation"):
            for t, c in unit.fluid_work_calculation.items():
                assert pytest.approx(1e-5, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "actual_work"):
            for t, c in unit.actual_work.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "isentropic_pressure"):
            for t, c in unit.isentropic_pressure.items():
                assert pytest.approx(1e-5, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "isothermal"):
            for t, c in unit.isothermal.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "isentropic"):
            for t, c in unit.isentropic.items():
                assert pytest.approx(1e-1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "isentropic_energy_balance"):
            for t, c in unit.isentropic_energy_balance.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "zero_work_equation"):
            for t, c in unit.zero_work_equation.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "state_material_balances"):
            for (t, j), c in unit.state_material_balances.items():
                assert pytest.approx(1e-2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        # heater and reactor only add 0D control volume constraints
        if hasattr(unit, 'material_holdup_calculation'):
            for (t, p, j), c in unit.material_holdup_calculation.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'rate_reaction_stoichiometry_constraint'):
            for (t, p, j), c in (
                    unit.rate_reaction_stoichiometry_constraint.items()):
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'equilibrium_reaction_stoichiometry_constraint'):
            for (t, p, j), c in (
                    unit.equilibrium_reaction_stoichiometry_constraint
                    .items()):
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'inherent_reaction_stoichiometry_constraint'):
            for (t, p, j), c in (
                    unit.inherent_reaction_stoichiometry_constraint.items()):
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'material_balances'):
            for (t, p, j), c in unit.material_balances.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'element_balances'):
            for (t, e), c in unit.element_balances.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'elemental_holdup_calculation'):
            for (t, e), c in unit.elemental_holdup_calculation.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'enthalpy_balances'):
            for t, c in unit.enthalpy_balances.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'energy_holdup_calculation'):
            for (t, p), c in unit.energy_holdup_calculation.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'pressure_balance'):
            for t, c in unit.pressure_balance.items():
                assert pytest.approx(1e-5, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, 'sum_of_phase_fractions'):
            for t, c in unit.sum_of_phase_fractions.items():
                assert pytest.approx(1e2, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)
        if hasattr(unit, "material_accumulation_disc_eq"):
            for (t, p, j), c in unit.material_accumulation_disc_eq.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "energy_accumulation_disc_eq"):
            for (t, p), c in unit.energy_accumulation_disc_eq.items():
                assert pytest.approx(1e-3, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

        if hasattr(unit, "element_accumulation_disc_eq"):
            for (t, e), c in unit.element_accumulation_disc_eq.items():
                assert pytest.approx(1, rel=1e-5) == \
                    iscale.get_constraint_transform_applied_scaling_factor(c)

    # equality constraints between ports at Arc sources and destinations
    for arc in model.fs.component_data_objects(Arc, descend_into=True):
        for c in arc.component_data_objects(Constraint, descend_into=True):
            if hasattr(unit, "enth_mol_equality"):
                for t, c in unit.enth_mol_equality.items():
                    assert pytest.approx(1e-3, rel=1e-5) == \
                        (iscale.
                         get_constraint_transform_applied_scaling_factor(c))
            if hasattr(unit, "flow_mol_equality"):
                for t, c in unit.flow_mol_equality.items():
                    assert pytest.approx(1, rel=1e-5) == \
                        (iscale.
                         get_constraint_transform_applied_scaling_factor(c))
            if hasattr(unit, "mole_frac_comp_equality"):
                for (t, j), c in unit.mole_frac_comp_equality.items():
                    assert pytest.approx(1e2, rel=1e-5) == \
                        (iscale.
                         get_constraint_transform_applied_scaling_factor(c))
            if hasattr(unit, "pressure_equality"):
                for t, c in unit.pressure_equality.items():
                    assert pytest.approx(1e-5, rel=1e-5) == \
                        (iscale.
                         get_constraint_transform_applied_scaling_factor(c))


@pytest.mark.integration
def test_initialize_flowsheet(model):
    initialize_flowsheet(model)

    assert degrees_of_freedom(model) == 0

    assert model.fs.M101.outlet.flow_mol[0].value == pytest.approx(
        954, 1e-3)
    assert model.fs.M102.outlet.flow_mol[0].value == pytest.approx(
        954.04, 1e-3)
    assert model.fs.C101.outlet.flow_mol[0].value == pytest.approx(
        954.04, 1e-3)
    assert model.fs.H101.outlet.flow_mol[0].value == pytest.approx(
        954.04, 1e-3)
    assert model.fs.R101.outlet.flow_mol[0].value == pytest.approx(
        478.82, 1e-3)
    assert model.fs.T101.outlet.flow_mol[0].value == pytest.approx(
        478.82, 1e-3)
    assert model.fs.H102.outlet.flow_mol[0].value == pytest.approx(
        478.82, 1e-3)
    assert model.fs.F101.vap_outlet.flow_mol[0].expr.value == pytest.approx(
        336.24, 1e-3)
    assert model.fs.F101.liq_outlet.flow_mol[0].expr.value == pytest.approx(
        142.59, 1e-3)

    assert not model.fs.M102.outlet.is_fixed()


@pytest.mark.unit
def test_unit_consistency(model):
    assert_units_consistent(model)


@pytest.mark.integration
def test_solve_flowsheet(model):
    solver = get_solver()
    optarg = {'tol': 1e-6,
              'max_iter': 5000}
    solver.options = optarg
    results = solver.solve(model, tee=False)
    assert results.solver.termination_condition == TerminationCondition.optimal

    assert value(model.fs.R101.rate_reaction_extent[0, "R1"]) == pytest.approx(
        237.6064, abs=1e-2)
    assert value(model.fs.R101.heat_duty[0])/1E6 == pytest.approx(
        -45.2203, abs=1e-2)
    assert value(model.fs.C101.work_mechanical[0])/1E6 == pytest.approx(
        -7.4506e-7, abs=1e-2)
    assert value(model.fs.T101.work_isentropic[0])/1E6 == pytest.approx(
        -0.9594, abs=1e-2)
    assert value(model.fs.F101.recovery*100) == pytest.approx(
        60.0063, abs=1e-2)
    assert value(model.fs.S101.split_fraction[0, "purge"]*100) == \
        pytest.approx(99.99, abs=1e-2)

    # check mass balance (inlets to outlets)
    feed_mass = value(model.fs.M101.H2_WGS.flow_mol[0] *
                      model.fs.thermo_params_vapor.H2.mw +
                      model.fs.M101.CO_WGS.flow_mol[0] *
                      model.fs.thermo_params_vapor.CO.mw)
    product_mass = value(sum(model.fs.F101.liq_outlet.flow_mol[0] *
                             model.fs.F101.liq_outlet.mole_frac_comp[0, comp]
                             * getattr(model.fs.thermo_params_VLE, comp).mw
                             for comp in ['CH4', 'CO', 'H2', 'CH3OH']))
    exhaust_mass = value(sum(model.fs.S101.purge.flow_mol[0] *
                             model.fs.S101.purge.mole_frac_comp[0, comp]
                             * getattr(model.fs.thermo_params_VLE, comp).mw
                             for comp in ['CH4', 'CO', 'H2', 'CH3OH']))

    assert value(feed_mass - product_mass - exhaust_mass) == \
        pytest.approx(0.0000, abs=1e-2)

    # check mass balances (on each unit)
    # mixer in is sum of named inlets, all others are total inlet
    # flash/splitter out are sum of named outlets, all others are total outlet
    for unit in ['M101', 'M102', 'C101', 'H101', 'R101',
                 'T101', 'H102', 'F101', 'S101']:
        block = getattr(model.fs, unit)

        # inlets
        if unit == 'M101':  # inlets are not named 'inlet'
            mass_in = feed_mass  # already have this
        elif unit == 'M102':
            mass_in = (feed_mass +  # already have this term
                       value(sum(model.fs.S101.recycle.flow_mol[0] *
                                 model.fs.S101.recycle.mole_frac_comp[0, comp]
                                 * getattr(model.fs.thermo_params_VLE, comp).mw
                                 for comp in ['CH4', 'CO', 'H2', 'CH3OH']))
                       )
        else:
            mass_in = sum(block.inlet.flow_mol[0] *
                          block.inlet.mole_frac_comp[0, comp]
                          * getattr(model.fs.thermo_params_VLE, comp).mw
                          for comp in ['CH4', 'CO', 'H2', 'CH3OH'])
        # outlets
        if unit == 'F101':  # outlets are not named 'outlet'
            mass_out = (
                product_mass +  # already have this term
                value(sum(model.fs.F101.vap_outlet.flow_mol[0] *
                          model.fs.F101.vap_outlet.mole_frac_comp[0, comp]
                          * getattr(model.fs.thermo_params_VLE, comp).mw
                          for comp in ['CH4', 'CO', 'H2', 'CH3OH']))
                      )
        elif unit == 'S101':  # outlets are not named 'outlet'
            mass_out = (
                exhaust_mass +  # already have this term
                value(sum(model.fs.S101.recycle.flow_mol[0] *
                          model.fs.S101.recycle.mole_frac_comp[0, comp]
                          * getattr(model.fs.thermo_params_VLE, comp).mw
                          for comp in ['CH4', 'CO', 'H2', 'CH3OH']))
                )
        else:
            mass_out = sum(block.outlet.flow_mol[0] *
                           block.outlet.mole_frac_comp[0, comp]
                           * getattr(model.fs.thermo_params_VLE, comp).mw
                           for comp in ['CH4', 'CO', 'H2', 'CH3OH'])

        assert value(mass_in - mass_out) == pytest.approx(0.0000, abs=1e-2)

    # check unit outlet temperatures and pressures

    assert value(model.fs.M101.mixed_state[0].temperature) == \
        pytest.approx(293.2017, abs=1e-2)
    assert value(model.fs.M102.mixed_state[0].temperature) == \
        pytest.approx(293.2017, abs=1e-2)
    assert value(model.fs.C101.control_volume.properties_out[0].temperature) == \
        pytest.approx(293.2017, abs=1e-2)
    assert value(model.fs.H101.control_volume.properties_out[0].temperature) == \
        pytest.approx(488.1500, abs=1e-2)
    assert value(model.fs.R101.control_volume.properties_out[0].temperature) == \
        pytest.approx(507.1500, abs=1e-2)
    assert value(model.fs.T101.control_volume.properties_out[0].temperature) == \
        pytest.approx(466.0282, abs=1e-2)
    assert value(model.fs.H102.control_volume.properties_out[0].temperature) == \
        pytest.approx(407.1500, abs=1e-2)
    assert value(model.fs.F101.control_volume.properties_out[0].temperature) == \
        pytest.approx(407.15, abs=1e-2)
    assert value(model.fs.S101.purge_state[0].temperature) == \
        pytest.approx(407.15, abs=1e-2)

    assert value(model.fs.M101.outlet.pressure[0]) == \
        pytest.approx(30e5, abs=1e-2)
    assert value(model.fs.M102.outlet.pressure[0]) == \
        pytest.approx(30e5, abs=1e-2)
    assert value(model.fs.C101.outlet.pressure[0]) == \
        pytest.approx(51e5, abs=1e-2)
    assert value(model.fs.H101.outlet.pressure[0]) == \
        pytest.approx(51e5, abs=1e-2)
    assert value(model.fs.R101.outlet.pressure[0]) == \
        pytest.approx(51e5, abs=1e-2)
    assert value(model.fs.T101.outlet.pressure[0]) == \
        pytest.approx(31e5, abs=1e-2)
    assert value(model.fs.H102.outlet.pressure[0]) == \
        pytest.approx(31e5, abs=1e-2)
    assert value(model.fs.F101.vap_outlet.pressure[0]) == \
        pytest.approx(31e5, abs=1e-2)
    assert value(model.fs.S101.purge.pressure[0]) == \
        pytest.approx(31e5, abs=1e-2)


@pytest.mark.unit
def test_add_costing(model):
    add_costing(model)

    assert hasattr(model.fs, "cooling_cost")
    assert hasattr(model.fs, "heating_cost")
    assert hasattr(model.fs, "electricity_cost")
    assert hasattr(model.fs, "operating_cost")
    assert hasattr(model.fs.R101, "L_eq")
    assert hasattr(model.fs.H101, "cost_heater")
    assert hasattr(model.fs.H102, "cost_heater")
    assert hasattr(model.fs, "annualized_capital_cost")
    assert hasattr(model.fs, "sales")
    assert hasattr(model.fs, "raw_mat_cost")
    assert hasattr(model.fs, "objective")
    assert degrees_of_freedom(model) == 0


@pytest.mark.integration
def test_optimize_with_costing(model):

    # Set up Optimization Problem (Maximize Revenue)
    # keep process pre-reaction fixed and unfix some post-process specs
    model.fs.R101.conversion.unfix()
    model.fs.R101.conversion_lb = Constraint(
        expr=model.fs.R101.conversion >= 0.75)
    model.fs.R101.conversion_ub = Constraint(
        expr=model.fs.R101.conversion <= 0.85)
    model.fs.R101.outlet_temp.deactivate()
    model.fs.R101.outlet_t_lb = Constraint(
        expr=model.fs.R101.control_volume.properties_out[0.0].temperature
        >= 405)
    model.fs.R101.outlet_t_ub = Constraint(
        expr=model.fs.R101.control_volume.properties_out[0.0].temperature
        <= 505)

    # Optimize turbine work (or delta P)
    model.fs.T101.deltaP.unfix()  # optimize turbine work/pressure drop
    model.fs.T101.outlet_p_lb = Constraint(
        expr=model.fs.T101.outlet.pressure[0] >= 10E5)
    model.fs.T101.outlet_p_ub = Constraint(
        expr=model.fs.T101.outlet.pressure[0] <= 51E5*0.8)

    # Optimize Cooler outlet temperature - unfix cooler outlet temperature
    model.fs.H102.outlet_temp.deactivate()
    model.fs.H102.outlet_t_lb = Constraint(
        expr=model.fs.H102.control_volume.properties_out[0.0].temperature
        >= 407.15*0.8)
    model.fs.H102.outlet_t_ub = Constraint(
        expr=model.fs.H102.control_volume.properties_out[0.0].temperature
        <= 480)

    model.fs.F101.deltaP.unfix()  # allow pressure change in streams

    model.fs.F101.isothermal = Constraint(
        expr=model.fs.F101.control_volume.properties_out[0].temperature ==
        model.fs.F101.control_volume.properties_in[0].temperature)

    model.fs.S101.split_fraction[0, "purge"].unfix()  # allow some gas recycle
    model.fs.S101.split_fraction_lb = Constraint(
        expr=model.fs.S101.split_fraction[0, "purge"] >= 0.10)  # min 10% purge
    model.fs.S101.split_fraction_ub = Constraint(
        expr=model.fs.S101.split_fraction[0, "purge"] <= 0.50)  # max 50% purge

    solver = get_solver()
    optarg = {'tol': 1e-6,
              'max_iter': 5000}
    solver.options = optarg
    results2 = solver.solve(model, tee=True)
    assert results2.solver.termination_condition == \
        TerminationCondition.optimal

    assert value(model.fs.R101.rate_reaction_extent[0, "R1"]) == pytest.approx(
        311.3070, abs=1e-2)
    assert value(model.fs.R101.heat_duty[0])/1E6 == pytest.approx(
        -59.3499, abs=1e-2)
    assert value(model.fs.C101.work_mechanical[0])/1E6 == pytest.approx(
        -1.7898e-24, abs=1e-2)
    assert value(model.fs.T101.work_isentropic[0])/1E6 == pytest.approx(
        -2.4215, abs=1e-2)
    assert value(model.fs.F101.recovery*100) == pytest.approx(
        95.1841, abs=1e-2)
    assert value(model.fs.S101.split_fraction[0, "purge"]*100) == \
        pytest.approx(10.0000, abs=1e-2)
    assert value(model.fs.objective)/1E6 == pytest.approx(
        104.7778, abs=1e-2)


@pytest.mark.unit
def test_report(model):
    # output is different for unit/integration test runs in order to test
    # report method during unit tests, and actual solution during integration

    if value(model.fs.R101.rate_reaction_extent[0, "R1"]) != 0.0:

        stream = StringIO()
        sys.stdout = stream
        report(model)
        sys.stdout = sys.__stdout__

        # this case contains solved solution values
        output = """
Extent of reaction:  311.3069854949999
Stoichiometry of each component normalized by the extent:
CH4 :  -0.0
H2 :  -2.0
CH3OH :  1.0
CO :  -1.0
These coefficients should follow 1*CO + 2*H2 => 1*CH3OH

Reaction conversion:  0.8500000099995839
Reactor duty (MW):  -59.350063183071555
Duty from Reaction (MW)): 28.21686516526679
Compressor work (MW):  -5.785653099883788e-25
Turbine work (MW):  -2.4215308912405678
Feed Mixer outlet temperature (C)):  20.051714213753655
Recycle Mixer outlet temperature (C)):  40.21858118534789
Feed Compressor outlet temperature (C)):  40.21858118534789
Feed Compressor outlet pressure (Pa)):  5100000.0
Heater outlet temperature (C)):  215.0
Reactor outlet temperature (C)):  231.8500047840716
Turbine outlet temperature (C)):  142.04194693695212
Turbine outlet pressure (Pa)):  1506510.241858315
Cooler outlet temperature (C)):  52.56999698828196
Flash outlet temperature (C)):  134.0
Purge percentage (amount of vapor vented to exhaust): 9.9999990000297  %
Methanol recovery(%):  95.18412474842862
annualized capital cost ($/year) = 284310.4928255764
operating cost ($/year) =  523119766.2459712
sales ($/year) =  140530623320.3268
raw materials cost ($/year) = 35229454878.16397
revenue (1000$/year)=  104777764.36542404


====================================================================================
Unit : fs.M101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
                                  Units        H2_WGS      CO_WGS     Outlet
    Total Molar Flowrate       mole / second     637.20      316.80     954.00
    Total Mole Fraction CH4    dimensionless 1.0000e-06  1.0000e-06 1.0000e-06
    Total Mole Fraction CO     dimensionless 1.0000e-06      1.0000    0.33208
    Total Mole Fraction H2     dimensionless     1.0000  1.0000e-06    0.66792
    Total Mole Fraction CH3OH  dimensionless 1.0000e-06  1.0000e-06 1.0000e-06
    Molar Enthalpy              joule / mole    -142.40 -1.1068e+05    -36848.
    Pressure                          pascal 3.0000e+06  3.0000e+06 3.0000e+06
====================================================================================

====================================================================================
Unit : fs.M102                                                             Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
                                  Units         feed      recycle    Outlet
    Total Molar Flowrate       mole / second     954.00     194.83     1148.8
    Total Mole Fraction CH4    dimensionless 1.0000e-06 4.4069e-05 8.3041e-06
    Total Mole Fraction CO     dimensionless    0.33208    0.25377    0.31880
    Total Mole Fraction H2     dimensionless    0.66792    0.67379    0.66892
    Total Mole Fraction CH3OH  dimensionless 1.0000e-06   0.072392   0.012278
    Molar Enthalpy              joule / mole    -36848.    -39572.    -37310.
    Pressure                          pascal 3.0000e+06 1.2103e+07 3.0000e+06
====================================================================================

====================================================================================
Unit : fs.F101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables:

    Key             : Value       : Units  : Fixed : Bounds
          Heat Duty : -9.9723e+06 :   watt : False : (None, None)
    Pressure Change :  1.0597e+07 : pascal : False : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                             Inlet    Vapor Outlet  Liquid Outlet
    flow_mol                  526.22       216.48        309.74
    mole_frac_comp CH4    1.8129e-05   4.4069e-05    1.0000e-08
    mole_frac_comp CO        0.10440      0.25377    1.0000e-08
    mole_frac_comp H2        0.27719      0.67379    1.0000e-08
    mole_frac_comp CH3OH     0.61839     0.072392        1.0000
    enth_mol             -1.3724e+05      -39572.   -2.3770e+05
    pressure              1.5065e+06   1.2103e+07    1.2103e+07
====================================================================================

====================================================================================
Unit : fs.S101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables:

    Key                           : Value   : Units         : Fixed : Bounds
    Split Fraction [('recycle',)] : 0.90000 : dimensionless : False : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                                  Units         Inlet
    Total Molar Flowrate       mole / second     216.48
    Total Mole Fraction CH4    dimensionless 4.4069e-05
    Total Mole Fraction CO     dimensionless    0.25377
    Total Mole Fraction H2     dimensionless    0.67379
    Total Mole Fraction CH3OH  dimensionless   0.072392
    Molar Enthalpy              joule / mole    -39572.
    Pressure                          pascal 1.2103e+07
====================================================================================
"""
        output_string = output.replace(" ", "").replace("\n", "")
        stream_string = str(stream.getvalue()).replace(" ", "").replace("\n", "")
        assert output_string == stream_string
    else:
        # model is not solved when integration tests are skipped, so need to
        # set a temporary nonzero extent value to allow unit testing of report
        # method and increase total code coverage from test module
        model.fs.R101.rate_reaction_extent[0, "R1"].fix(0.01)

        stream = StringIO()
        sys.stdout = stream
        report(model)
        sys.stdout = sys.__stdout__

        # this case is pre-initialization, pre-solve and contain defaults
        output = """
Extent of reaction:  0.01
Stoichiometry of each component normalized by the extent:
CH4 :  0.0
H2 :  0.0
CH3OH :  0.0
CO :  0.0
These coefficients should follow 1*CO + 2*H2 => 1*CH3OH

Reaction conversion:  0.75
Reactor duty (MW):  0.0
Duty from Reaction (MW)): 0.0009064
Compressor work (MW):  0.0
Turbine work (MW):  0.0
Feed Mixer outlet temperature (C)):  25.0
Recycle Mixer outlet temperature (C)):  25.0
Feed Compressor outlet temperature (C)):  25.0
Feed Compressor outlet pressure (Pa)):  5100000.0
Heater outlet temperature (C)):  25.0
Reactor outlet temperature (C)):  25.0
Turbine outlet temperature (C)):  25.0
Turbine outlet pressure (Pa)):  100000.0
Cooler outlet temperature (C)):  25.0
Flash outlet temperature (C)):  25.0
Purge percentage (amount of vapor vented to exhaust): 99.99  %
Methanol recovery(%):  1.0
annualized capital cost ($/year) = 121776.27770012307
operating cost ($/year) =  0.0
sales ($/year) =  5671299423.599999
raw materials cost ($/year) = 35229454878.16397
revenue (1000$/year)=  -29558277.23084167


====================================================================================
Unit : fs.M101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
                                  Units        H2_WGS      CO_WGS     Outlet
    Total Molar Flowrate       mole / second     637.20      316.80     100.00
    Total Mole Fraction CH4    dimensionless 1.0000e-06  1.0000e-06    0.25000
    Total Mole Fraction CO     dimensionless 1.0000e-06      1.0000    0.25000
    Total Mole Fraction H2     dimensionless     1.0000  1.0000e-06    0.25000
    Total Mole Fraction CH3OH  dimensionless 1.0000e-06  1.0000e-06    0.25000
    Molar Enthalpy              joule / mole    -142.40 -1.1068e+05     100.00
    Pressure                          pascal 3.0000e+06  3.0000e+06 1.0000e+05
====================================================================================

====================================================================================
Unit : fs.M102                                                             Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
                                  Units         feed      recycle    Outlet
    Total Molar Flowrate       mole / second     100.00     100.00     100.00
    Total Mole Fraction CH4    dimensionless    0.25000    0.25000    0.25000
    Total Mole Fraction CO     dimensionless    0.25000    0.25000    0.25000
    Total Mole Fraction H2     dimensionless    0.25000    0.25000    0.25000
    Total Mole Fraction CH3OH  dimensionless    0.25000    0.25000    0.25000
    Molar Enthalpy              joule / mole     100.00     100.00     100.00
    Pressure                          pascal 1.0000e+05 1.0000e+05 1.0000e+05
====================================================================================

====================================================================================
Unit : fs.F101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables:

    Key             : Value  : Units  : Fixed : Bounds
          Heat Duty : 0.0000 :   watt : False : (None, None)
    Pressure Change : 0.0000 : pascal :  True : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                            Inlet    Vapor Outlet  Liquid Outlet
    flow_mol                 100.00       50.000        50.000
    mole_frac_comp CH4      0.25000      0.25000    1.0000e-08
    mole_frac_comp CO       0.25000      0.25000    1.0000e-08
    mole_frac_comp H2       0.25000      0.25000    1.0000e-08
    mole_frac_comp CH3OH    0.25000      0.25000       0.25000
    enth_mol                 100.00      -97532.       -59600.
    pressure             1.0000e+05   1.0000e+05    1.0000e+05
====================================================================================

====================================================================================
Unit : fs.S101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables:

    Key                           : Value   : Units         : Fixed : Bounds
    Split Fraction [('recycle',)] : 0.50000 : dimensionless : False : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                                  Units         Inlet
    Total Molar Flowrate       mole / second     100.00
    Total Mole Fraction CH4    dimensionless    0.25000
    Total Mole Fraction CO     dimensionless    0.25000
    Total Mole Fraction H2     dimensionless    0.25000
    Total Mole Fraction CH3OH  dimensionless    0.25000
    Molar Enthalpy              joule / mole     100.00
    Pressure                          pascal 1.0000e+05
====================================================================================
"""
        output_string = output.replace(" ", "").replace("\n", "")
        stream_string = str(stream.getvalue()).replace(" ", "").replace("\n", "")
        assert output_string == stream_string
