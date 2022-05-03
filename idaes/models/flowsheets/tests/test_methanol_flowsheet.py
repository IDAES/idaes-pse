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

from idaes.models.flowsheets.methanol_flowsheet import (
    build_model, set_inputs, initialize_flowsheet, add_costing, report)

from pyomo.environ import (Constraint,
                           ConcreteModel,
                           value)
from pyomo.network import Arc
from pyomo.environ import TerminationCondition
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

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
    Flash)


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
    assert isinstance(model.fs.C101, Compressor)
    assert isinstance(model.fs.H101, Heater)
    assert isinstance(model.fs.R101, StoichiometricReactor)
    assert isinstance(model.fs.T101, Turbine)
    assert isinstance(model.fs.H102, Heater)
    assert isinstance(model.fs.F101, Flash)

    assert isinstance(model.fs.s02, Arc)
    assert isinstance(model.fs.s03, Arc)
    assert isinstance(model.fs.s04, Arc)
    assert isinstance(model.fs.s05, Arc)
    assert isinstance(model.fs.s06, Arc)
    assert isinstance(model.fs.s07, Arc)

    assert degrees_of_freedom(model) == 23


@pytest.mark.unit
def test_set_inputs(model):
    set_inputs(model)

    assert degrees_of_freedom(model) == 0


@pytest.mark.integration
def test_initialize_flowsheet(model):
    initialize_flowsheet(model)

    assert degrees_of_freedom(model) == 0

    assert model.fs.M101.outlet.flow_mol[0].value == pytest.approx(954, 1e-3)
    assert model.fs.C101.outlet.flow_mol[0].value == pytest.approx(954, 1e-3)
    assert model.fs.H101.outlet.flow_mol[0].value == pytest.approx(954, 1e-3)
    assert model.fs.R101.outlet.flow_mol[0].value == pytest.approx(478.8, 1e-3)
    assert model.fs.T101.outlet.flow_mol[0].value == pytest.approx(478.8, 1e-3)
    assert model.fs.H102.outlet.flow_mol[0].value == pytest.approx(478.8, 1e-3)
    assert model.fs.F101.vap_outlet.flow_mol[0].expr.value == pytest.approx(
        336.23, 1e-3)
    assert model.fs.F101.liq_outlet.flow_mol[0].expr.value == pytest.approx(
        142.57, 1e-3)


@pytest.mark.integration
def test_unit_consistency(model):
    assert_units_consistent(model)


@pytest.mark.integration
def test_solve_flowsheet(model):
    solver = get_solver()
    optarg = {'tol': 1e-6,
              'max_iter': 5000}
    solver.options = optarg
    results = solver.solve(model, tee=True)
    assert results.solver.termination_condition == TerminationCondition.optimal

    assert value(model.fs.R101.rate_reaction_extent[0, "R1"]) == pytest.approx(
        237.6005, abs=1e-2)
    assert value(model.fs.R101.heat_duty[0])/1E6 == pytest.approx(
        -45.2192, abs=1e-2)
    assert value(model.fs.T101.work_isentropic[0])/1E6 == pytest.approx(
        -0.9593, abs=1e-2)
    assert value(model.fs.F101.recovery*100) == pytest.approx(
        60.0047, abs=1e-2)

    # check mass balance (inlets to outlets)
    feed_mass = value(model.fs.M101.H2_WGS.flow_mol[0] *
                      model.fs.thermo_params_vapor.H2.mw +
                      model.fs.M101.CO_WGS.flow_mol[0] *
                      model.fs.thermo_params_vapor.CO.mw)
    product_mass = value(sum(model.fs.F101.liq_outlet.flow_mol[0] *
                             model.fs.F101.liq_outlet.mole_frac_comp[0, comp]
                             * getattr(model.fs.thermo_params_VLE, comp).mw
                             for comp in ['CH4', 'CO', 'H2', 'CH3OH']))
    exhaust_mass = value(sum(model.fs.F101.vap_outlet.flow_mol[0] *
                             model.fs.F101.vap_outlet.mole_frac_comp[0, comp]
                             * getattr(model.fs.thermo_params_VLE, comp).mw
                             for comp in ['CH4', 'CO', 'H2', 'CH3OH']))

    assert value(feed_mass - product_mass - exhaust_mass) == \
        pytest.approx(0.0000, abs=1e-2)

    # check mass balances (on each unit)
    # mixer in is sum of named inlets, all others are total inlet
    # flash out is sum of named outlets, all others are total outlet
    for unit in ['M101', 'C101', 'H101', 'R101', 'T101', 'H102', 'F101']:
        block = getattr(model.fs, unit)

        # inlets
        if unit == 'M101':  # inlets are not named 'inlet'
            mass_in = feed_mass  # already have this
        else:
            mass_in = sum(block.inlet.flow_mol[0] *
                          block.inlet.mole_frac_comp[0, comp]
                          * getattr(model.fs.thermo_params_VLE, comp).mw
                          for comp in ['CH4', 'CO', 'H2', 'CH3OH'])

        # outlets
        if unit == 'F101':  # outlets are not named 'outlet'
            mass_out = product_mass + exhaust_mass  # already have these
        else:
            mass_out = sum(block.outlet.flow_mol[0] *
                           block.outlet.mole_frac_comp[0, comp]
                           * getattr(model.fs.thermo_params_VLE, comp).mw
                           for comp in ['CH4', 'CO', 'H2', 'CH3OH'])

        assert value(mass_in - mass_out) == pytest.approx(0.0000, abs=1e-2)

    # check unit outlet temperatures and pressures

    assert value(model.fs.M101.mixed_state[0].temperature) == \
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

    assert value(model.fs.M101.outlet.pressure[0]) == \
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

    solver = get_solver()
    optarg = {'tol': 1e-6,
              'max_iter': 5000}
    solver.options = optarg
    results2 = solver.solve(model, tee=True)
    assert results2.solver.termination_condition == \
        TerminationCondition.optimal

    assert value(model.fs.R101.rate_reaction_extent[0, "R1"]) == pytest.approx(
        269.2805, abs=1e-2)
    assert value(model.fs.R101.heat_duty[0])/1E6 == pytest.approx(
        -51.3636, abs=1e-2)
    assert value(model.fs.T101.work_isentropic[0])/1E6 == pytest.approx(
        -1.9905, abs=1e-2)
    assert value(model.fs.F101.recovery*100) == pytest.approx(
        95.5433, abs=1e-2)
    assert value(model.fs.objective)/1E6 == pytest.approx(
        81.0476, abs=1e-2)


@pytest.mark.unit
def test_report(model):
    report(model)
