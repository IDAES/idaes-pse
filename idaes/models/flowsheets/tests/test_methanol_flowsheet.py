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

from idaes.models.flowsheets.methanol_flowsheet import (
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
    solver = get_solver()
    optarg = {'tol': 1e-6,
              'max_iter': 5000}
    solver.options = optarg
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


@pytest.mark.integration
def test_report(model):
    stream = StringIO()
    sys.stdout = stream
    report(model)
    sys.stdout = sys.__stdout__

    output = """

Extent of reaction:  269.28054478798873
Stoichiometry of each component normalized by the extent:
CH4 :  0.0
H2 :  -2.0
CH3OH :  1.0
CO :  -1.0
These coefficients should follow 1*CO + 2*H2 => 1*CH3OH

Reaction conversion:  0.8500000099999442
Reactor duty (MW):  -51.36357357754515
Duty from Reaction (MW)): 24.407588579583297
Turbine work (MW):  -1.9904899177794635
Mixer outlet temperature (C)):  20.051714213753485
Compressor outlet temperature (C)):  20.051714213753485
Compressor outlet pressure (Pa)):  5100000.0
Heater outlet temperature (C)):  215.0
Reactor outlet temperature (C)):  231.85000468716584
Turbine outlet temperature (C)):  139.8588817267576
Turbine outlet pressure (Pa)):  1427653.3547821408
Cooler outlet temperature (C)):  52.56999709299214
Flash outlet temperature (C)):  134.0
Methanol recovery(%):  95.5432880783248
annualized capital cost ($/year) = 262011.81501550632
operating cost ($/year) =  451845877.93162763
sales ($/year) =  116729133888.84218
raw materials cost ($/year) = 35229454878.16397
revenue (1000$/year)=  81047571.12093157


====================================================================================
Unit : fs.M101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
                                H2_WGS      CO_WGS     Outlet  
    Total Molar Flowrate          637.20      316.80     954.00
    Total Mole Fraction CH4   1.0000e-06  1.0000e-06 1.0000e-06
    Total Mole Fraction CO    1.0000e-06      1.0000    0.33208
    Total Mole Fraction H2        1.0000  1.0000e-06    0.66792
    Total Mole Fraction CH3OH 1.0000e-06  1.0000e-06 1.0000e-06
    Molar Enthalpy               -142.40 -1.1068e+05    -36848.
    Pressure                  3.0000e+06  3.0000e+06 3.0000e+06
====================================================================================

====================================================================================
Unit : fs.F101                                                             Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key             : Value       : Fixed : Bounds
          Heat Duty : -8.3431e+06 : False : (None, None)
    Pressure Change :  1.0119e+07 : False : (None, None)

------------------------------------------------------------------------------------
    Stream Table
                             Inlet    Vapor Outlet  Liquid Outlet
    flow_mol                  415.44       158.16        257.28  
    mole_frac_comp CH4    2.2964e-06   6.0318e-06    1.0000e-08  
    mole_frac_comp CO        0.11438      0.30045    1.0000e-08  
    mole_frac_comp H2        0.23743      0.62366    1.0000e-08  
    mole_frac_comp CH3OH     0.64818     0.075879        1.0000  
    enth_mol             -1.4444e+05      -45435.   -2.3772e+05  
    pressure              1.4277e+06   1.1547e+07    1.1547e+07  
====================================================================================
"""
    assert output in stream.getvalue()
