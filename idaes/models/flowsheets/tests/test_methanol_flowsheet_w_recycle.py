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
from idaes.core.util import get_solver
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


@pytest.mark.integration
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


@pytest.mark.integration
def test_optimize_with_costing(model):
    add_costing(model)

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
    report(model)
