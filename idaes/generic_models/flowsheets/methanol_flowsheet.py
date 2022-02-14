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
                           maximize)
from pyomo.environ import TerminationCondition
from pyomo.network import Arc

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
    Flash)
from idaes.generic_models.unit_models.mixer import MomentumMixingType
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
import idaes.core.util.unit_costing as costing


def build_model(m):
    # Define model components and blocks
    m.fs = FlowsheetBlock(default={"dynamic": False})

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

    # Build the flowsheet
    print('Unit degrees of freedom')
    for unit in ('M101', 'C101', 'H101', 'R101', 'T101', 'H102', 'F101'):
        if unit == 'M101':
            spec = 14
        else:
            spec = 7
        print(str(unit)+' '+str(degrees_of_freedom(getattr(m.fs, unit))-spec))

    # pre-compression
    m.fs.s02 = Arc(source=m.fs.M101.outlet, destination=m.fs.C101.inlet)

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

    print('DOF after units specified: ', degrees_of_freedom(m))


def initialize_flowsheet(m):

    # Initialize and solve flowsheet

    print('')
    m.fs.M101.initialize()
    propagate_state(arc=m.fs.s02)  # mixer to compressor

    m.fs.C101.initialize()
    propagate_state(arc=m.fs.s03)  # compressor to heater

    m.fs.H101.initialize()
    propagate_state(arc=m.fs.s04)  # heater to reactor

    m.fs.R101.initialize()
    propagate_state(arc=m.fs.s05)  # reactor to turbine

    m.fs.T101.initialize()
    propagate_state(arc=m.fs.s06)  # turbine to cooler

    m.fs.H102.initialize()
    propagate_state(arc=m.fs.s07)  # cooler to flash

    m.fs.F101.initialize()


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
    print('Turbine work (MW): ', m.fs.T101.work_isentropic[0].value/1e6)
    print('Mixer outlet temperature (C)): ',
          m.fs.M101.mixed_state[0].temperature.value - 273.15)
    print('Compressor outlet temperature (C)): ',
          m.fs.C101.control_volume.properties_out[0].temperature.value
          - 273.15)
    print('Compressor outlet pressure (Pa)): ',
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
    print('Methanol recovery(%): ', value(100*m.fs.F101.recovery))
    print('annualized capital cost ($/year) =',
          value(m.fs.annualized_capital_cost))
    print('operating cost ($/year) = ', value(m.fs.operating_cost))
    print('sales ($/year) = ', value(m.fs.sales))
    print('raw materials cost ($/year) =', value(m.fs.raw_mat_cost))
    print('revenue (1000$/year)= ', value(m.fs.objective))

    print()
    m.fs.M101.report()
    m.fs.F101.report()
