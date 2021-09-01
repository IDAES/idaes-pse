##############################################################################
# The development of this flowsheet/code is funded by the ARPA-E DIFFERENTIATE
# project: “Machine Learning for Natural Gas to Electric Power System Design”
# Project number: DE-FOA-0002107-1625.
# This project is a collaborative effort between the Pacific Northwest National
# Laboratory, the National Energy Technology Laboratory, and the University of
# Washington to design NGFC systems with high efficiencies and low CO2
# emissions.
##############################################################################
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
IDAES block for co2 purifucation process based on surrogate models

This is a black box (surrogate based) model for CO2 purification process
available in Aspen.

The unit has 1 inlet stream for rich CO2 inlet: inlet
the unti has 3 outlet streams: pureCO2, water, and vent

The state variables are flow_mol, mol_fraction, temperature, and pressure.
The state variables can be accessed thotough ports named:
    inlet
    pure_CO2
    water
    vent_gas

The surrogate models were prepared varying rich inlet flowrate. Expected value
11350 lbmol/hr, lower bound: 10850 lbmol/hr, and upper bound: 11850 lbmol/hr.
300 samples have been developed using Latin Hypercube Sampling Method, the data
set is trained by running the CPU rigorous model. ALAMO has been used to fit
the surrogate models.
Outputs are function of (inlet flowrate, inlet composition).
This model 1 is simple and temperature and pressure are fixed

pyomo block for co2 purifucation process based on surrogates

This is a black box (surrogate based) model for CO2 purification process
available in Aspen.

The unit has 1 inlet stream for rich CO2 inlet: inlet
Also, the unti has 3 outlet streams: pureCO2, water, and vent

The state variables are:
    (1) flow_mol
    (2) temperature
    (3) pressure
    (4) mole_frac_comp
The state variables can be accessed thotough ports named:
    inlet
    pureco2
    water
    vent

This model is only based on the inlet flow rate of the CO2 rich stream
The degrees of freedom are the Inlet states:
    (1) Inlet flow in mol/s
    (2) Inlet temperature in K
    (3) Inlet pressure in Pa
    (4) Inlet component mol fractions: CO2, H2O, N2, Ar, O2
"""
# Import IDAES cores
import numpy as np
from idaes.core import declare_process_block_class, UnitModelBlockData
import idaes.core.util.scaling as iscale
# Additional import for the unit operation
import pyomo.environ as pyo
from pyomo.environ import Var, value, units as pyunits, Constraint
from pyomo.network import Port
import idaes.logger as idaeslog

import idaes.power_generation.flowsheets.sofc.surrogates.cpu_surrogate_methods as sm

__author__ = "Differentiate Team (N. Susarla, A. Noring, M. Zamarripa)"
__version__ = "1.0.0"


# ----------------------------------------------------------------------------
@declare_process_block_class("CPU")
class CPUData(UnitModelBlockData):
    '''
    CO2Pure surrogate model based on total flow of CO2 rich stream flow only
    Assumptions:
    Fixed composition, temperature, and pressure of feed stream
    Fixed CO2 purity in the CO2 product stream

    '''

    def build(self):

        self.component_list = ['Ar', 'CO2', 'O2', 'H2O', 'N2']

        self.make_vars()
        self.add_material_balances()
        self.add_surrogates()

        # Add ports: 4 (1 for inlet and 3 for outlets)
        self.inlet = Port(noruleinit=True,
                          doc="A port for co2 rich inlet stream")
        self.pureco2 = Port(noruleinit=True,
                             doc="A port for pure CO2 outlet stream")
        self.water = Port(noruleinit=True,
                          doc="A port for water outlet stream")
        self.vent = Port(noruleinit=True,
                             doc="A port for vent gas outlet stream")

        # Add state vars to the ports
        self.inlet.add(self.inlet_flow_mol, "flow_mol")
        self.inlet.add(self.inlet_temperature, "temperature")
        self.inlet.add(self.inlet_pressure, "pressure")
        self.inlet.add(self.inlet_mole_frac_comp, "mole_frac_comp")

        self.pureco2.add(self.pureco2_flow_mol, "flow_mol")
        self.pureco2.add(self.pureco2_temperature, "temperature")
        self.pureco2.add(self.pureco2_pressure, "pressure")
        self.pureco2.add(self.pureco2_mole_frac_comp, "mole_frac_comp")

        self.water.add(self.water_flow_mol, "flow_mol")
        self.water.add(self.water_temperature, "temperature")
        self.water.add(self.water_pressure, "pressure")
        self.water.add(self.water_mole_frac_comp, "mole_frac_comp")

        self.vent.add(self.vent_flow_mol, "flow_mol")
        self.vent.add(self.vent_temperature, "temperature")
        self.vent.add(self.vent_pressure, "pressure")
        self.vent.add(self.vent_mole_frac_comp, "mole_frac_comp")

    def make_vars(self):
        ''' This section is for creating all the vars for this model.
            There are 1 inlet and 3 outlet streams.
            These streams are names as inlet, pureco2, water, vent
            For each of these streams, the following variables are defined:
            (1) Total mole flow [mol/s]: [stream_name]_flow_mol
            (2) Component molar fraction: [stream_name]_mole_frac_comp
            (3) Component mole flows [mol/s]: [stream_name]_flow_mol_comp
            (4) Temperature [K]: [stream_name]_temperature
            (5) Pressure [Pa]: [stream_name]_pressure

        '''
        # units declaration for vars
        flow_units = pyunits.mol/pyunits.s
        pressure_units = pyunits.Pa
        temperature_units = pyunits.K
        # heat_duty_units = pyunits.J/pyunits.s

        # Total mole flow [mol/s]
        self.inlet_flow_mol = Var(
            self.flowsheet().config.time,
            initialize=3600,
            # bounds=(1367.08, 5000),
            units=flow_units,
            doc='Total inlet mole flow [mol/s]')
        self.pureco2_flow_mol = Var(
            self.flowsheet().config.time,
            initialize=900,
            units=flow_units,
            doc='PureCO2 outlet stream mole flow [mol/s]')
        self.water_flow_mol = Var(
            self.flowsheet().config.time,
            initialize=300,
            units=flow_units,
            doc='Water outlet stream mole flow [mol/s]')
        self.vent_flow_mol = Var(
            self.flowsheet().config.time,
            initialize=2400,
            units=flow_units,
            doc='Vent gas stream mole flow [mol/s]')

        # Component molar fractions
        self.inlet_mole_frac_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=1/len(self.component_list),
            doc='Inlet stream: Component mole fraction')
        self.pureco2_mole_frac_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=1/len(self.component_list),
            doc='PureCO2 stream: Component mole fraction')
        self.water_mole_frac_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=1/len(self.component_list),
            doc='Water stream: Component mole fraction')
        self.vent_mole_frac_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=1/len(self.component_list),
            doc='Vent stream: Component mole fraction')

        # Component mole flows [mol/s]
        self.inlet_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=3600/len(self.component_list),
            units=flow_units,
            doc='Inlet stream: Component mole flow [mol/s]')
        self.pureco2_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=900/len(self.component_list),
            units=flow_units,
            doc='PureCO2 stream: Component mole flow [mol/s]')
        self.water_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=300/len(self.component_list),
            units=flow_units,
            doc='Water stream: Component mole flow [mol/s]')
        self.vent_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=2400/len(self.component_list),
            units=flow_units,
            doc='Vent stream: Component mole flow [mol/s]')

        # Temperature [K]
        self.inlet_temperature = Var(self.flowsheet().config.time,
                                     initialize=110,
                                     units=temperature_units,
                                     doc='Inlet temperature [K]')
        self.pureco2_temperature = Var(self.flowsheet().config.time,
                                       initialize=110,
                                       units=temperature_units,
                                       doc='PureCO2 temperature [K]')
        self.water_temperature = Var(self.flowsheet().config.time,
                                     initialize=110,
                                     units=temperature_units,
                                     doc='Water temperature [K]')
        self.vent_temperature = Var(self.flowsheet().config.time,
                                    initialize=110,
                                    units=temperature_units,
                                    doc='Vent temperature [K]')

        # Pressue [Pa]
        self.inlet_pressure = Var(self.flowsheet().config.time,
                                  initialize=17,
                                  units=pressure_units,
                                  doc='Inlet pressure [Pa]')
        self.pureco2_pressure = Var(self.flowsheet().config.time,
                                    initialize=17,
                                    units=pressure_units,
                                    doc='PureCO2 pressure [Pa]')
        self.water_pressure = Var(self.flowsheet().config.time,
                                  initialize=17,
                                  units=pressure_units,
                                  doc='Water pressure [Pa]')
        self.vent_pressure = Var(self.flowsheet().config.time,
                                 initialize=17,
                                 units=pressure_units,
                                 doc='Vent pressure [Pa]')

    def add_material_balances(self):
        ''' This section is for material balance constraints'''

        # Sum of all componenet mole fractions in a stream equals 1
        @self.Constraint(self.flowsheet().config.time,
                         doc="PureCO2 stream: component mole flow equation")
        def mole_frac_comp_pureco2_eqn(b, t):
            return (
                0 == 1 - sum(b.pureco2_mole_frac_comp[t, c]
                             for c in self.component_list))
        @self.Constraint(self.flowsheet().config.time,
                         doc="Water stream: component mole flow equation")
        def mole_frac_comp_water_eqn(b, t):
            return (
                0 == 1 - sum(b.water_mole_frac_comp[t, c]
                             for c in self.component_list))
        @self.Constraint(self.flowsheet().config.time,
                         doc="Vent stream: component mole flow equation")
        def mole_frac_comp_vent_eqn(b, t):
            return (
                0 == 1 - sum(b.vent_mole_frac_comp[t, c]
                             for c in self.component_list))

        # Component mole flow = total flow_mol * mole_frac_comp
        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="Inlet component mole flow eqn [mol/s]")
        def flow_mol_comp_inlet_eqn(b, t, c):
            return (
                b.inlet_flow_mol_comp[t, c] ==
                b.inlet_flow_mol[t] * b.inlet_mole_frac_comp[t, c])

        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="PureCO2 component mole flow eqn [mol/s]")
        def flow_mol_comp_pureco2_eqn(b, t, c):
            return (
                b.pureco2_flow_mol_comp[t, c] ==
                b.pureco2_flow_mol[t] * b.pureco2_mole_frac_comp[t, c])

        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="Water component mole flow eqn [mol/s]")
        def flow_mol_comp_water_eqn(b, t, c):
            return (
                b.water_flow_mol_comp[t, c] ==
                b.water_flow_mol[t] * b.water_mole_frac_comp[t, c])

        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="Vent component mole flow eqn [mol/s]")
        def flow_mol_comp_vent_eqn(b, t, c):
            return (
                b.vent_flow_mol_comp[t, c] ==
                b.vent_flow_mol[t] * b.vent_mole_frac_comp[t, c])

        # Total component balance across the unit (i. e. inlet and outlets)
        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="component material balance equation")
        def component_material_balance_eqn(b, t, c):
            return (
                0 == self.inlet_flow_mol_comp[t, c] -
                self.pureco2_flow_mol_comp[t, c] -
                self.water_flow_mol_comp[t, c] -
                self.vent_flow_mol_comp[t, c])

    def add_surrogates(self):
        ''' This section is to add the surrogate models'''


        # Compressure heat duty
        @self.Expression(self.flowsheet().config.time,
                          doc="Compressor train heat duty [J/s]")
        def heat_duty(b, t):
            return (
                sm.heat_duty_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Compressor train work
        @self.Expression(self.flowsheet().config.time,
                          doc="Compressor train work [W]")
        def work(b, t):
            return (
                sm.compressor_power_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Refrigeration duty
        @self.Expression(self.flowsheet().config.time,
                          doc="Refrigeration duty [W]")
        def refrigeration_duty(b, t):
            return (
                sm.refrigeration_duty_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Pure CO2 stream: Temperature
        @self.Constraint(self.flowsheet().config.time,
                          doc="Temperature: pureCO2 [K]")
        def pureco2_temperature_eq(b, t):
            return (
                b.pureco2_temperature[t] ==
                sm.pureco2_temperature_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Pure CO2 stream: Pressure
        @self.Constraint(self.flowsheet().config.time,
                          doc="Pressure: pureCO2 [Pa]")
        def pureco2_pressure_eq(b, t):
            return (
                b.pureco2_pressure[t] ==
                sm.pureco2_pressure_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))


        # Pure CO2 stream: Total flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="Flow_mol Total: pureCO2 [mol/s]")
        def pureco2_total_flow_eq(b, t):
            return (
                b.pureco2_flow_mol[t] ==
                sm.pureco2_flow_mol_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Pure CO2 stream: CO2 mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="co2 component flow_mol in pureco2")
        def pureco2_co2_flow_mol_comp_eq(b, t):
            return (
                b.pureco2_flow_mol_comp[t, 'CO2'] ==
                sm.pureco2_co2_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))


        # Pure CO2 stream: O2 mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="o2 component flow_mol in pureco2")
        def pureco2_o2_flow_mol_comp_eq(b, t):
            return (
                b.pureco2_flow_mol_comp[t, 'O2'] ==
                sm.pureco2_o2_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))


        # Pure CO2 stream: Ar mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="ar component flow_mol in pureco2")
        def pureco2_ar_flow_mol_comp_eq(b, t):
            return (
                b.pureco2_flow_mol_comp[t, 'Ar'] ==
                sm.pureco2_ar_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))


        # Pure CO2 stream: H2O mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="h2o component flow_mol in pureco2")
        def pureco2_h2o_flow_mol_comp_eq(b, t):
            return (
                b.pureco2_flow_mol_comp[t, 'H2O'] ==
                sm.pureco2_h2o_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))


        # Water stream: Total flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="Flow_mol Total: water [mol/s]")
        def water_total_flow_eq(b, t):
            return (
                b.water_flow_mol[t] ==
                sm.water_flow_mol_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Water stream: Temperature
        @self.Constraint(self.flowsheet().config.time,
                          doc="Temperature: water [K]")
        def water_temperature_eq(b, t):
            return (
                b.water_temperature[t] ==
                sm.water_temperature_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Water stream: Pressure
        @self.Constraint(self.flowsheet().config.time,
                          doc="Pressure: water [Pa]")
        def water_pressure_eq(b, t):
            return (
                b.water_pressure[t] ==
                sm.water_pressure_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Water stream: CO2 mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="co2 component flow_mol in water")
        def water_co2_flow_mol_comp_eq(b, t):
            return (
                b.water_flow_mol_comp[t, 'CO2'] ==
                sm.water_co2_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Water stream: O2 mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="o2 component flow_mol in water")
        def water_o2_flow_mol_comp_eq(b, t):
            return (
                b.water_flow_mol_comp[t, 'O2'] ==
                sm.water_o2_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Water stream: AR mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="ar component flow_mol in water")
        def water_ar_flow_mol_comp_eq(b, t):
            return (
                b.water_flow_mol_comp[t, 'Ar'] ==
                sm.water_ar_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Water stream: H2O mole flow
        @self.Constraint(self.flowsheet().config.time,
                          doc="h2o component flow_mol in water")
        def water_h2o_flow_mol_comp_eq(b, t):
            return (
                b.water_flow_mol_comp[t, 'H2O'] ==
                sm.water_h2o_flow_mol_comp_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))


        # Vent stream: Temperature
        @self.Constraint(self.flowsheet().config.time,
                          doc="Temperature: vent [K]")
        def vent_temperature_eq(b, t):
            return (
                b.vent_temperature[t] ==
                sm.vent_temperature_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Vent stream: Pressure
        @self.Constraint(self.flowsheet().config.time,
                          doc="Pressure: vent [Pa]")
        def vent_pressure_eq(b, t):
            return (
                b.vent_pressure[t] ==
                sm.vent_pressure_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

    def initialize(blk,
                   outlvl=idaeslog.NOTSET,
                   solver='ipopt',
                   optarg={'tol': 1e-6}):
        '''
        CO2 pure pyomo block initialization routine

        Keyword Arguments:
            outlvl : sets output level of initialisation routine

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''
        iscale.calculate_scaling_factors(blk)
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        opt = pyo.SolverFactory(solver)
        opt.options = optarg

        init_log.info_low("Starting initialization...")

        blk.inlet.flow_mol[0].fix()
        blk.inlet.mole_frac_comp[0, 'Ar'].fix()
        blk.inlet.mole_frac_comp[0, 'CO2'].fix()
        blk.inlet.mole_frac_comp[0, 'O2'].fix()
        blk.inlet.mole_frac_comp[0, 'H2O'].fix()
        blk.inlet.mole_frac_comp[0, 'N2'].fix()

        # solve model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 1 {}.".format(idaeslog.condition(res))
            )
        init_log.info_high("Initialization Step 1 Complete.")

        # check component material balances
        # deactivate component balance equation for inlet
        # fix all inlet flows
        # blk.inlet.mole_frac_comp[0, 'O2'].fix(0)
        # blk.mole_frac_comp_inlet_eqn.deactivate()
        # with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
        #     res = opt.solve(blk, tee=slc.tee)
        # init_log.info_high(
        #         "Initialization Step 2 {}.".format(idaeslog.condition(res))
        #     )

        # release all states and activate constraint
        blk.inlet.flow_mol[0].unfix()
        blk.inlet.mole_frac_comp[0, 'Ar'].unfix()
        blk.inlet.mole_frac_comp[0, 'CO2'].unfix()
        blk.inlet.mole_frac_comp[0, 'O2'].unfix()
        blk.inlet.mole_frac_comp[0, 'H2O'].unfix()
        blk.inlet.mole_frac_comp[0, 'N2'].unfix()
        # blk.mole_frac_comp_inlet_eqn.activate()

        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # for v in self.component_data_objects(Var):
        #     iscale.set_scaling_factor(v, 1)
        for c in self.component_data_objects(Constraint):
            iscale.set_scaling_factor(c, 1)

        iscale.set_scaling_factor(self.inlet_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0,'CO2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0,'H2O'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0,'O2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0,'Ar'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0,'N2'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0,'CO2'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0,'H2O'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0,'O2'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0,'Ar'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0,'N2'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0,'CO2'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0,'H2O'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0,'O2'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0,'Ar'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0,'N2'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0,'CO2'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0,'H2O'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0,'O2'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0,'Ar'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0,'N2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.pureco2_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.water_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.vent_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.inlet_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.pureco2_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.water_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.vent_pressure[0.0], 1e-5)
