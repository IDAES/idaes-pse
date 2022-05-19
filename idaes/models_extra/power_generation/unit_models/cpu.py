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
##############################################################################
# The development of this flowsheet/code is funded by the ARPA-E DIFFERENTIATE
# project: “Machine Learning for Natural Gas to Electric Power System Design”
# Project number: DE-FOA-0002107-1625.
# This project is a collaborative effort between the Pacific Northwest National
# Laboratory, the National Energy Technology Laboratory, and the University of
# Washington to design NGFC systems with high efficiencies and low CO2
# emissions.
##############################################################################
"""
IDAES block for surrogate models of a CO2 purifucation process

This is a black box (surrogate based) model for CO2 purification process
designed in Aspen.

The unit has 1 inlet stream for CO2-rich flue gas: inlet
The unit has 3 outlet streams: pure_CO2, water, and vent

The state variables are flow_mol, mol_fraction, temperature, and pressure.
The state variables can be accessed via ports named: inlet, pure_CO2, water,
and vent_gas

The surrogate models were prepared by varying the inlet flowrate and
composition. The expected value is 11350 lbmol/hr, lower bound: 10850 lbmol/hr,
and upper bound: 11850 lbmol/hr.
300 samples have been developed using Latin Hypercube Sampling
Method, the data set is trained by running the CPU rigorous model. ALAMO was
used to fit the surrogate models.
Outputs are a function of the inlet flowrate and species composition.
Inlet temperature and pressure are assumed to be fixed at 100 F and 14.7 psia.
The degrees of freedom are the inlet states:
(1) Inlet flow in mol/s
(2) Inlet component mol fractions: CO2, H2O, N2, Ar, O2
"""

# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
import idaes.core.util.scaling as iscale
# Additional import for the unit operation
import pyomo.environ as pyo
from pyomo.environ import Var, Constraint, units as pyunits
from pyomo.network import Port
import idaes.logger as idaeslog

__author__ = "Differentiate Team (N. Susarla, A. Noring, M. Zamarripa)"
__version__ = "1.1.0"


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
        work_units = pyunits.J/pyunits.s

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

        # Heat duty and work [W]
        self.heat_duty = Var(self.flowsheet().config.time,
                             initialize=0,
                             units=work_units,
                             doc='CO2 compressor cooling duty [W]')

        self.work = Var(self.flowsheet().config.time,
                        initialize=0,
                        units=work_units,
                        doc='CO2 compressor work [W]')

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

        # Compressor heat duty
        @self.Constraint(self.flowsheet().config.time,
                         doc="Compressor train heat duty [J/s]")
        def heat_duty_eq(b, t):
            return (
                b.heat_duty[t] ==
                heat_duty_fun.f(
                    b.inlet_flow_mol[t],
                    b.inlet_mole_frac_comp[t, 'Ar'],
                    b.inlet_mole_frac_comp[t, 'CO2'],
                    b.inlet_mole_frac_comp[t, 'O2'],
                    b.inlet_mole_frac_comp[t, 'H2O'],
                    b.inlet_mole_frac_comp[t, 'N2']
                    ))

        # Compressor train work
        @self.Constraint(self.flowsheet().config.time,
                         doc="Compressor train work [W]")
        def work_eq(b, t):
            return (
                b.work[t] ==
                compressor_power_fun.f(
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
                pureco2_temperature_fun.f(
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
                pureco2_pressure_fun.f(
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
                pureco2_flow_mol_fun.f(
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
                pureco2_co2_flow_mol_comp_fun.f(
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
                pureco2_o2_flow_mol_comp_fun.f(
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
                pureco2_ar_flow_mol_comp_fun.f(
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
                pureco2_h2o_flow_mol_comp_fun.f(
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
                water_flow_mol_fun.f(
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
                water_temperature_fun.f(
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
                water_pressure_fun.f(
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
                water_co2_flow_mol_comp_fun.f(
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
                water_o2_flow_mol_comp_fun.f(
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
                water_ar_flow_mol_comp_fun.f(
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
                water_h2o_flow_mol_comp_fun.f(
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
                vent_temperature_fun.f(
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
                vent_pressure_fun.f(
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
                   optarg={'tol': 1e-6},
                   release_state=True):
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

        # release all states and activate constraint
        if release_state:
            blk.inlet.flow_mol[0].unfix()
            blk.inlet.mole_frac_comp[0, 'Ar'].unfix()
            blk.inlet.mole_frac_comp[0, 'CO2'].unfix()
            blk.inlet.mole_frac_comp[0, 'O2'].unfix()
            blk.inlet.mole_frac_comp[0, 'H2O'].unfix()
            blk.inlet.mole_frac_comp[0, 'N2'].unfix()

        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        print("yo")
        super().calculate_scaling_factors()

        for c in self.component_data_objects(Constraint):
            iscale.set_scaling_factor(c, 1)

        iscale.set_scaling_factor(self.inlet_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol[0.0], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'CO2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'H2O'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'O2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'Ar'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'N2'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0, 'CO2'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0, 'H2O'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0, 'O2'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0, 'Ar'], 1e-3)
        iscale.set_scaling_factor(self.pureco2_flow_mol_comp[0.0, 'N2'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0, 'CO2'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0, 'H2O'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0, 'O2'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0, 'Ar'], 1e-3)
        iscale.set_scaling_factor(self.water_flow_mol_comp[0.0, 'N2'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0, 'CO2'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0, 'H2O'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0, 'O2'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0, 'Ar'], 1e-3)
        iscale.set_scaling_factor(self.vent_flow_mol_comp[0.0, 'N2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.pureco2_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.water_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.vent_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.inlet_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.pureco2_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.water_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.vent_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.heat_duty, 1e-6)
        iscale.set_scaling_factor(self.work, 1e-6)


class compressor_power_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 16592.309504827182536246 * x1 + 17806584.568199951201677 * x3 - 15529374.443028066307306


class heat_duty_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 27431.344430845900205895 * x1 + 25040167.356691692024469 * x3 - 21837134.557537198066711


class pureco2_flow_mol_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.84491492199656215156267 * x1 + 1381.8486099781443954271 * x3 - 1204.9518833417739642755


class pureco2_ar_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.87539452852376016910194E-007 * x1 + 0.31019041198388018487631E-001 * x2 - 0.16384578403083416373726E-002 * x5 - 0.11401484796994311183421E-002 * x6


class pureco2_co2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.84491226339466729555738 * x1 + 1381.8776912884593457420 * x3 - 1204.9772340007211823831


class pureco2_o2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.25557869821781261075201E-005 * x1 - 0.30160307289038564698691E-001 * x3 + 0.26298816057506283622169E-001


class pureco2_h2o_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.00


class pureco2_n2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.69246160597454559832190E-006 * x3 - 0.56005883473002498644045E-006


class pureco2_temperature_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 10.267585306350582641244 * x2 + 32.023906991174357017371 * x3 + 30.302544593098598824099 * x4 + 278.86088998572887476257


class pureco2_pressure_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.28666413419197245189076E-011 * x1 + 15271893.400000000372529


class water_flow_mol_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.64898275330398047522351E-001 * x1 + 995.23584715669267097837 * x4 - 64.595577093201725915605


class water_ar_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.23736607897806708417605E-008 * x1 + 0.83689507888515507049582E-003 * x2 - 0.14749984695782550662112E-003 * x5


class water_co2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.36504473074529521385892E-005 * x1 + 0.39821667908625070844697E-002 * x3 + 0.30283921974618616818065E-001 * x4 - 0.54375872162499300915828E-002


class water_o2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.46405004009352215417855E-007 * x1 - 0.99119387109367312720110E-003 * x3 + 0.86473506719933338534462E-003


class water_h2o_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.64894576227660510925332E-001 * x1 + 995.21157061827375400753 * x4 - 64.594007608678808196601


class water_n2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.71432912211586127618749E-009 * x1


class water_temperature_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.47345180660722969989695E-002 * x3 - 1.1291949419088465056404 * x4 - 0.13703395454789617582958E-001 * x5 + 311.08988060641769379799


class water_pressure_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 117210.91999999999825377


class vent_flow_mol_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.90250925128294370858306E-001 * x1 - 1373.8223230605847220431 * x3 + 1197.7852419317096064333


class vent_ar_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.27178896806122057258626E-002 * x1 + 953.43145974385458885081 * x2 - 168.06049873202732669597 * x5


class vent_co2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.27083917235717330090905E-001 * x1 - 407.85524106058562665567 * x3 + 355.59132055263916072363


class vent_o2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.44846813052226472406936E-001 * x1 - 956.47536271014735120843 * x3 + 834.42603654916342748038


class vent_h2o_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.00


class vent_n2_flow_mol_comp_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return 0.15024528735008041077648E-001 * x1


class vent_temperature_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.43741475554195015242121E-016 * x1 + 305.37222222222220580079


class vent_pressure_fun:
    def f(x1, x2, x3, x4, x5, x6):
        return - 0.22395635483747847803966E-013 * x1 + 115831.96800000000803266
