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
Task: IDAES Support for ARPE-E Differentiate
Scenario: Methanol Synthesis From Syngas
Author: B. Paul and M. Zamarripa
"""

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Objective,
                           Var,
                           Expression,
                           ConcreteModel,
                           TransformationFactory,
                           value,
                           maximize,
                           units as pyunits)
from pyomo.environ import TerminationCondition
from pyomo.network import Arc, SequentialDecomposition

# Import IDAES core libraries
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state

# Import required models

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.generic.generic_reaction import \
    GenericReactionParameterBlock

from idaes.generic_models.properties.examples import \
    methanol_water_ideal as thermo_props
from idaes.generic_models.properties.examples import \
    methanol_reactions as reaction_props

from idaes.generic_models.unit_models import (
    Mixer,
    Heater,
    Compressor,
    Turbine,
    StoichiometricReactor,
    Flash,
    Separator as Splitter)
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
import idaes.core.util.unit_costing as costing
import idaes.logger as idaeslog


def build_model(m):
    # Define model components and blocks
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # lower upper state bound on temperature to allow convergence with recycle
    thermo_props.config_dict["state_bounds"]["temperature"] = \
        (198.15, 298.15, 512.15, pyunits.K)

    m.fs.thermo_params = GenericParameterBlock(
        default=thermo_props.config_dict)
    m.fs.reaction_params = GenericReactionParameterBlock(
        default={"property_package": m.fs.thermo_params,
                 **reaction_props.config_dict})

    # mixing feed streams
    m.fs.M101 = Mixer(
        default={"property_package": m.fs.thermo_params,
                 "momentum_mixing_type": MomentumMixingType.minimize,
                 "has_phase_equilibrium": True,
                 "inlet_list": ['H2_WGS', 'CO_WGS']})

    # mixing recycle back into feed streams
    m.fs.M102 = Mixer(
        default={"property_package": m.fs.thermo_params,
                 "momentum_mixing_type": MomentumMixingType.minimize,
                 "has_phase_equilibrium": True,
                 "inlet_list": ["feed", "recycle"]})

    # pre-compression
    m.fs.C101 = Compressor(
        default={"dynamic": False,
                 "property_package": m.fs.thermo_params,
                 "compressor": True,
                 "thermodynamic_assumption": ThermodynamicAssumption.isothermal
                 })

    # pre-heating
    m.fs.H101 = Heater(
        default={"property_package": m.fs.thermo_params,
                 "has_pressure_change": False,
                 "has_phase_equilibrium": False})

    # reactor
    m.fs.R101 = StoichiometricReactor(
        default={"has_heat_transfer": True,
                 "has_heat_of_reaction": True,
                 "has_pressure_change": False,
                 "property_package": m.fs.thermo_params,
                 "reaction_package": m.fs.reaction_params})

    # post-expansion
    m.fs.T101 = Turbine(
        default={"dynamic": False,
                 "property_package": m.fs.thermo_params})

    # post-cooling
    m.fs.H102 = Heater(
        default={"property_package": m.fs.thermo_params,
                 "has_pressure_change": False,
                 "has_phase_equilibrium": False})

    # product recovery
    m.fs.F101 = Flash(
        default={"property_package": m.fs.thermo_params,
                 "has_heat_transfer": True,
                 "has_pressure_change": True})

    # split waste gases for recycle
    m.fs.S101 = Splitter(
        default={"property_package": m.fs.thermo_params,
                 "ideal_separation": False,
                 "outlet_list": ["purge", "recycle"]})

    # Build the flowsheet
    print('Unit degrees of freedom')
    for unit in ('M101', 'C101', 'H101', 'R101', 'T101', 'H102', 'F101',
                 'M102', 'S101'):
        if unit == 'M101' or unit == 'M102':
            spec = 14
        else:
            spec = 7
        print(str(unit)+' '+str(degrees_of_freedom(getattr(m.fs, unit))-spec))

    # mixed feed to mix with recycle
    m.fs.s01 = Arc(source=m.fs.M101.outlet, destination=m.fs.M102.feed)

    # pre-compression
    m.fs.s02 = Arc(source=m.fs.M102.outlet, destination=m.fs.C101.inlet)

    # pre-heating
    m.fs.s03 = Arc(source=m.fs.C101.outlet, destination=m.fs.H101.inlet)

    # reactor feed
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)

    # post-expansion
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.T101.inlet)

    # post-cooling
    m.fs.s06 = Arc(source=m.fs.T101.outlet, destination=m.fs.H102.inlet)

    # product recovery
    m.fs.s07 = Arc(source=m.fs.H102.outlet, destination=m.fs.F101.inlet)

    # waste gases
    m.fs.s08 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)

    # recycle
    m.fs.s09 = Arc(source=m.fs.S101.recycle, destination=m.fs.M102.recycle)

    # connecting unit models
    TransformationFactory("network.expand_arcs").apply_to(m)

    # Add unit and stream specifications
    print('Total DOF: ', degrees_of_freedom(m))
    return m


def set_inputs(m):

    #  feed streams, post WGS
    m.fs.M101.H2_WGS.flow_mol[0].fix(637.2)  # mol/s, relative to 177 kmol/h
    m.fs.M101.H2_WGS.mole_frac_comp[0, "H2"].fix(1)
    m.fs.M101.H2_WGS.mole_frac_comp[0, "CO"].fix(1e-6)
    m.fs.M101.H2_WGS.mole_frac_comp[0, "CH3OH"].fix(1e-6)
    m.fs.M101.H2_WGS.mole_frac_comp[0, "CH4"].fix(1e-6)
    m.fs.M101.H2_WGS.mole_frac_comp[0, "H2O"].fix(1e-6)
    m.fs.M101.H2_WGS.enth_mol[0].fix(-142.4)  # J/mol
    m.fs.M101.H2_WGS.pressure.fix(30e5)  # Pa

    m.fs.M101.CO_WGS.flow_mol[0].fix(316.8)  # mol/s, relative to 88 kmol/h
    m.fs.M101.CO_WGS.mole_frac_comp[0, "H2"].fix(1e-6)
    m.fs.M101.CO_WGS.mole_frac_comp[0, "CO"].fix(1)
    m.fs.M101.CO_WGS.mole_frac_comp[0, "CH3OH"].fix(1e-6)
    m.fs.M101.CO_WGS.mole_frac_comp[0, "CH4"].fix(1e-6)
    m.fs.M101.CO_WGS.mole_frac_comp[0, "H2O"].fix(1e-6)
    m.fs.M101.CO_WGS.enth_mol[0].fix(-110676.4)  # J/mol
    m.fs.M101.CO_WGS.pressure.fix(30e5)  # Pa
    print('DOF after streams specified: ', degrees_of_freedom(m))

    # units specifications
    m.fs.C101.outlet.pressure.fix(51e5)  # Pa

    m.fs.H101.outlet_temp = Constraint(
        expr=m.fs.H101.control_volume.properties_out[0].temperature == 488.15)

    m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))
    m.fs.R101.conv_constraint = Constraint(
        expr=(m.fs.R101.conversion * m.fs.R101.inlet.flow_mol[0] *
              m.fs.R101.inlet.mole_frac_comp[0, "CO"] ==
              m.fs.R101.inlet.flow_mol[0] *
              m.fs.R101.inlet.mole_frac_comp[0, "CO"]
              - m.fs.R101.outlet.flow_mol[0] *
              m.fs.R101.outlet.mole_frac_comp[0, "CO"]))
    m.fs.R101.conversion.fix(0.75)
    m.fs.R101.outlet_temp = Constraint(
        expr=m.fs.R101.control_volume.properties_out[0].temperature == 507.15)
    m.fs.R101.heat_duty.setub(0)  # rxn is exothermic, so duty is cooling only

    m.fs.T101.deltaP.fix(-2e6)
    m.fs.T101.efficiency_isentropic.fix(0.9)

    m.fs.H102.outlet_temp = Constraint(
        expr=m.fs.H102.control_volume.properties_out[0].temperature == 407.15)

    m.fs.F101.recovery = Var(initialize=0.01, bounds=(0, 1))
    m.fs.F101.rec_constraint = Constraint(
        expr=(m.fs.F101.recovery == m.fs.F101.liq_outlet.flow_mol[0] *
              m.fs.F101.liq_outlet.mole_frac_comp[0, "CH3OH"] /
              (m.fs.F101.inlet.flow_mol[0] *
               m.fs.F101.inlet.mole_frac_comp[0, "CH3OH"])))
    m.fs.F101.deltaP.fix(0)  # Pa
    m.fs.F101.outlet_temp = Constraint(
        expr=m.fs.F101.control_volume.properties_out[0].temperature == 407.15)

    m.fs.S101.split_fraction[0, "purge"].fix(0.9999)  # initially no recycle

    print('DOF after units specified: ', degrees_of_freedom(m))


def initialize_flowsheet(m, solved_m=None):

    # Initialize and solve flowsheet

    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5

    # Using the SD tool
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)
    print()
    print('Tear Stream:')
    for o in heuristic_tear_set:
        print(o.name, ': ', o.source.name, ' to ', o.destination.name)
    print()
    print('Calculation order:')
    for o in order:
        print(o[0].name)
    print()

    # automatically build stream set for flowsheet and find the tear stream
    stream_set = [arc for arc in m.fs.component_data_objects(Arc)]

    # if available, pull tear guesses from solved model
    if solved_m is not None:
        # use the prior solution values for the tear stream
        tear_name = heuristic_tear_set[0].name.split(sep='.')
        tear = getattr(solved_m.fs, tear_name[1]).destination
        tear_guesses = {
            "flow_mol": {0: round(value(tear.flow_mol[0]), 2)},
            "mole_frac_comp": {
                (0, "CH4"): round(value(tear.mole_frac_comp[0, "CH4"]), 6),
                (0, "CO"): round(value(tear.mole_frac_comp[0, "CO"]), 5),
                (0, "H2"): round(value(tear.mole_frac_comp[0, "H2"]), 5),
                (0, "CH3OH"): round(value(tear.mole_frac_comp[0, "CH3OH"]), 6),
                (0, "H2O"): round(value(tear.mole_frac_comp[0, "H2O"]), 6)},
            "enth_mol": {0: round(value(tear.enth_mol[0]), 0)},
            "pressure": {0: round(value(tear.pressure[0]), 0)}}
        print('Tear guesses: ')
        print(tear_guesses)
    else:  # manually enter tear guesses - assumes set_inputs() above and
        # that solver selects s02 (to fs.C101.inlet) as the tear stream
        tear_guesses = {
            "flow_mol": {0: 954.00},
            "mole_frac_comp": {
                    (0, "CH4"): 1e-6,
                    (0, "CO"): 0.33207,
                    (0, "H2"): 0.66792,
                    (0, "CH3OH"): 1e-6,
                    (0, "H2O"): 1e-6},
            "enth_mol": {0: -36848},
            "pressure": {0: 3e6}}

    # set tear guesses
    for stream in stream_set:
        if stream in heuristic_tear_set:
            seq.set_guesses_for(stream.destination, tear_guesses)

    def function(unit):
        print('Solving ', str(unit))
        unit.initialize(outlvl=idaeslog.ERROR)  # no output unless it breaks
        for stream in stream_set:
            if stream.source.parent_block() == unit:
                propagate_state(arc=stream)  # this is an outlet of the unit
            stream.destination.unfix()
        print('DOF = ', degrees_of_freedom(m))
    print('Initial DOF = ', degrees_of_freedom(m))
    seq.run(m, function)
    for stream in stream_set:
        if stream.destination.is_fixed() is True:
            print('Unfixing ', stream.destination.name, '...')
            stream.destination.unfix()
    print('Final DOF = ', degrees_of_freedom(m))


def add_costing(m):

    from idaes.core.util.unit_costing import initialize as init_costing
    assert degrees_of_freedom(m) == 0

    # Expression to compute the total cooling cost (F/R cooling not assumed)
    m.fs.cooling_cost = Expression(
        expr=0.25e-7 * (-m.fs.F101.heat_duty[0])
        + 0.212e-7 * (-m.fs.H102.heat_duty[0])
        + 2.2e-7 * (-m.fs.R101.heat_duty[0]))

    # Expression to compute the total heating cost (F/R heating not assumed)
    m.fs.heating_cost = Expression(
        expr=2.2e-7 * m.fs.H101.heat_duty[0])

    # Expression to compute the total electricity cost (utilities - credit)
    m.fs.electricity_cost = Expression(
        expr=0.12e-5 * (m.fs.C101.work_mechanical[0])
        - 0.08e-5 * (m.fs.T101.work_isentropic[0]))

    # Expression to compute the total operating cost
    m.fs.operating_cost = Expression(
        expr=(3600 * 24 * 365 * (m.fs.heating_cost + m.fs.cooling_cost
                                 + m.fs.electricity_cost)))

    # Computing reactor capital cost
    m.fs.R101.get_costing()
    m.fs.R101.diameter.fix(2)
    m.fs.R101.length.fix(4)  # for initial problem at 75% conversion
    init_costing(m.fs.R101.costing)
    # Reactor length (size, and capital cost) is adjusted based on conversion
    # surrogate model which scales length linearly with conversion
    m.fs.R101.length.unfix()
    m.fs.R101.L_eq = Constraint(expr=m.fs.R101.length ==
                                13.2000*m.fs.R101.conversion - 5.9200)

    # Computing flash capital cost
    m.fs.F101.get_costing()
    m.fs.F101.diameter.fix(2)
    m.fs.F101.length.fix(4)
    init_costing(m.fs.F101.costing)

    # Computing heater/cooler capital costs
    # Surrogates prepared with IDAES shell and tube hx considering IP steam and
    # assuming steam outlet is condensed
    m.fs.H101.cost_heater = Expression(
        expr=0.036158*m.fs.H101.heat_duty[0] + 63931.475,
        doc='capital cost of heater in $')

    # Surrogates prepared with IDAES shell and tube hx considering cooling
    # water assuming that water inlet T is 25 deg C and outlet T is 40 deg C
    m.fs.H102.cost_heater = Expression(
        expr=0.10230*(-m.fs.H102.heat_duty[0]) + 100421.572,
        doc='capital cost of cooler in $')

    # Annualizing capital cost to same scale as operating costs (per year)
    m.fs.annualized_capital_cost = Expression(
        expr=(m.fs.R101.costing.purchase_cost
              + m.fs.F101.costing.purchase_cost
              + m.fs.H101.cost_heater
              + m.fs.H102.cost_heater)*5.4/15)

    # methanol price $449 us dollars per metric ton  - 32.042 g/mol
    # - 1 gr = 1e-6 MT  -- consider 1000
    # H2 $16.51 per kilogram - 2.016 g/mol
    # CO $62.00 per kilogram - 28.01 g/mol
    m.fs.sales = Expression(
        expr=(3600 * 24 * 365
              * m.fs.F101.liq_outlet.flow_mol[0]
              * m.fs.F101.liq_outlet.mole_frac_comp[0, "CH3OH"]
              * 32.042 * 1e-6 * 449 * 1000
              )
        )
    m.fs.raw_mat_cost = Expression(
        expr=(3600 * 24 * 365 * m.fs.M101.CO_WGS.flow_mol[0] *
              16.51 * 2.016 / 1000
              + 3600 * 24 * 365 * m.fs.M101.H2_WGS.flow_mol[0] *
              62.00 * 28.01 / 1000
              )
        )

    m.fs.objective = Objective(expr=(m.fs.sales
                                     - m.fs.operating_cost
                                     - m.fs.annualized_capital_cost
                                     - m.fs.raw_mat_cost)/1e3,
                               sense=maximize)
    assert degrees_of_freedom(m) == 0
    costing.calculate_scaling_factors(m.fs.F101.costing)
    costing.calculate_scaling_factors(m.fs.R101.costing)


def report(m):

    # Display some results

    print()
    print()
    extent = m.fs.R101.rate_reaction_extent[0, "R1"]  # shorter parameter alias
    print('Extent of reaction: ', value(extent))
    print('Stoichiometry of each component normalized by the extent:')
    complist = ('CH4', 'H2', 'H2O', 'CH3OH', 'CO')
    changelist = [value(m.fs.R101.outlet.mole_frac_comp[0, comp] *
                        m.fs.R101.outlet.flow_mol[0] -
                        m.fs.R101.inlet.mole_frac_comp[0, comp] *
                        m.fs.R101.inlet.flow_mol[0]) for comp in complist]
    normlist = [changelist[i]/value(extent) for i in range(len(complist))]
    change = dict(zip(complist, normlist))
    for entry in change:
        print(entry, ': ', round(change[entry], 2))
    print('These coefficients should follow 1*CO + 2*H2 => 1*CH3OH')

    print()
    print('Reaction conversion: ', m.fs.R101.conversion.value)
    print('Reactor duty (MW): ', m.fs.R101.heat_duty[0].value/1e6)
    print('Duty from Reaction (MW)):', value(extent) *
          -m.fs.reaction_params.reaction_R1.dh_rxn_ref.value/1e6)
    print('Compressor work (MW): ', m.fs.C101.work_mechanical[0].value/1e6)
    print('Turbine work (MW): ', m.fs.T101.work_isentropic[0].value/1e6)
    print('Feed Mixer outlet temperature (C)): ',
          m.fs.M101.mixed_state[0].temperature.value - 273.15)
    print('Recycle Mixer outlet temperature (C)): ',
          m.fs.M102.mixed_state[0].temperature.value - 273.15)
    print('Feed Compressor outlet temperature (C)): ',
          m.fs.C101.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Feed Compressor outlet pressure (Pa)): ',
          m.fs.C101.control_volume.properties_out[0].pressure.value)
    print('Heater outlet temperature (C)): ',
          m.fs.H101.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Reactor outlet temperature (C)): ',
          m.fs.R101.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Turbine outlet temperature (C)): ',
          m.fs.T101.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Turbine outlet pressure (Pa)): ',
          m.fs.T101.control_volume.properties_out[0].pressure.value)
    print('Cooler outlet temperature (C)): ',
          m.fs.H102.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Flash outlet temperature (C)): ',
          m.fs.F101.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Purge percentage (amount of vapor vented to exhaust):',
          100*value(m.fs.S101.split_fraction[0, "purge"]), ' %')
    print('Methanol recovery(%): ', value(100*m.fs.F101.recovery))
    print('annualized capital cost ($/year) =',
          value(m.fs.annualized_capital_cost))
    print('operating cost ($/year) = ', value(m.fs.operating_cost))
    print('sales ($/year) = ', value(m.fs.sales))
    print('raw materials cost ($/year) =', value(m.fs.raw_mat_cost))
    print('revenue (1000$/year)= ', value(m.fs.objective))

    print()
    m.fs.M101.report()
    m.fs.M102.report()
    m.fs.F101.report()
    m.fs.S101.report()
