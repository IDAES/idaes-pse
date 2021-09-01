################################################################################
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
################################################################################
"""
IDAES carbon capture system block, hybrid mass balances and surrogate models

The unit has 1 inlet stream for rich CO2 inlet: inlet
the unit has 2 outlet streams: pureCO2 and exhaust_gas

The state variables are component flow mol, temperature, and pressure.
The state variables can be accessed thotough ports named:
    inlet - flue gas from NGCC plant
    pure_CO2 - stream for compression train
    exhaust_gas - exhaust gas to DAC system or stack

Surrogate models are used to compute operating costs as a function of solvent
flow rate, and specific reboiler duty. The inlet variables to estimate such
output variables are: lean loading and PZ molality.
Inlet Vars:

* solvent lean loading (0.1 to 0.6)
* solvent molality (PZ = 3, 5, 7 mol)

Output Vars:

* L/G ratio for absorber columns (0 - 12)
* lean solvent flowrate in kmol/s (0-100)
* specific reboiler duty (GJ/t CO2)

Reference:
Gaspar J. von Solms N., Thomsen K., Fosbol F. (2016) Multivariable Optimization
 of the Piperazine CO2 Post-Combustion Process. Energy Procedia 86(2016)229-238
"""
# Import IDAES cores
from idaes.core import declare_process_block_class, UnitModelBlockData
import idaes.core.util.scaling as iscale
# Additional import for the unit operation
import pyomo.environ as pyo
from pyomo.environ import Var, units as pyunits, Constraint, exp, log
from pyomo.network import Port
import idaes.logger as idaeslog

__author__ = "M. Zamarripa"
__version__ = "1.0.0"


# ----------------------------------------------------------------------------
@declare_process_block_class("CO2Capture")
class CO2CaptureData(UnitModelBlockData):
    '''
    CO2Capture surrogate model based on total flow of CO2 rich stream flow only
    Assumptions: (toDo: update Assumptions)
    Fixed composition, temperature, and pressure of feed stream
    Fixed CO2 purity in the CO2 product stream

    '''

    def build(self):

        self.component_list = ['CO2', 'H2O', 'O2', 'N2', 'Ar']

        self.make_vars()
        self.add_material_balances()
        # self.add_surrogates()
        # self.add_costing()

        # Add ports: 3 (1 for inlet and 2 for outlets)
        self.inlet = Port(noruleinit=True,
                          doc="A port for co2 rich inlet stream")
        self.pureCO2 = Port(noruleinit=True,
                            doc="A port for pure CO2 outlet stream")
        self.exhaust_gas = Port(noruleinit=True,
                                doc="A port for vent gas outlet stream")

        # Add state vars to the ports
        # self.inlet.add(self.inlet_flow_mol, "flow_mol")
        self.inlet.add(self.inlet_temperature, "temperature")
        self.inlet.add(self.inlet_pressure, "pressure")
        self.inlet.add(self.inlet_flow_mol_comp, "flow_mol_comp")

        # self.pureco2.add(self.pureco2_flow_mol, "flow_mol")
        self.pureCO2.add(self.pureCO2_temperature, "temperature")
        self.pureCO2.add(self.pureCO2_pressure, "pressure")
        self.pureCO2.add(self.pureCO2_flow_mol_comp, "flow_mol_comp")

        # self.exhaust_gas.add(self.exhaust_gas_flow_mol, "flow_mol")
        self.exhaust_gas.add(self.exhaust_gas_temperature, "temperature")
        self.exhaust_gas.add(self.exhaust_gas_pressure, "pressure")
        self.exhaust_gas.add(self.exhaust_gas_flow_mol_comp, "flow_mol_comp")

    def make_vars(self):
        '''
        This section builds port vars (Fc, T, P), CO2 capture rate

        '''
        # units declaration for vars
        flow_units = pyunits.mol/pyunits.s
        pressure_units = pyunits.Pa
        temperature_units = pyunits.K
        # heat_duty_units = pyunits.J/pyunits.s

        # Component mole flows [mol/s]
        self.inlet_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=1400/len(self.component_list),
            units=flow_units,
            doc='Inlet stream: Component mole flow [mol/s]')

        self.pureCO2_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=1200/len(self.component_list),
            units=flow_units,
            doc='PureCO2 stream: Component mole flow [mol/s]')

        self.exhaust_gas_flow_mol_comp = Var(
            self.flowsheet().config.time,
            self.component_list,
            initialize=100/len(self.component_list),
            units=flow_units,
            doc='exhaust_gas stream: Component mole flow [mol/s]')

        # Temperature [K]
        self.inlet_temperature = Var(self.flowsheet().config.time,
                                     initialize=110,
                                     units=temperature_units,
                                     doc='Inlet temperature [K]')

        self.pureCO2_temperature = Var(self.flowsheet().config.time,
                                       initialize=110,
                                       units=temperature_units,
                                       doc='PureCO2 temperature [K]')

        self.exhaust_gas_temperature = Var(self.flowsheet().config.time,
                                           initialize=110,
                                           units=temperature_units,
                                           doc='exhaust_gas temperature [K]')

        # Pressue [Pa]
        self.inlet_pressure = Var(self.flowsheet().config.time,
                                  initialize=17,
                                  units=pressure_units,
                                  doc='Inlet pressure [Pa]')

        self.pureCO2_pressure = Var(self.flowsheet().config.time,
                                    initialize=17,
                                    units=pressure_units,
                                    doc='PureCO2 pressure [Pa]')

        self.exhaust_gas_pressure = Var(self.flowsheet().config.time,
                                        initialize=17,
                                        units=pressure_units,
                                        doc='exhaust_gas pressure [Pa]')
        # CO2 Capture rate
        self.CO2_capture_rate = Var(self.flowsheet().config.time,
                                    initialize=0.9,
                                    doc='CO2 capture rate')
        # lean loading
        self.lean_loading = Var(self.flowsheet().config.time,
                                initialize=0.1,
                                bounds=(0.1, 6.0),
                                doc='lean loading')
        # Pz molality
        self.Pz_mol = Var(self.flowsheet().config.time,
                        # domain=pyo.Integers,
                          initialize=3,
                          bounds=(3, 7),
                          doc='Pz molality')

        self.SRD = Var(self.flowsheet().config.time,
                       initialize=8,
                       bounds=(0, 100),
                       doc='Specific reformer duty GJ/ton CO2')

    def add_material_balances(self):
        ''' This section is for material balance constraints'''

        # pureCO2 mass balance
        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="pureCO2 mass balances")
        def pureCO2_eqn(b, t, c):
            if c == "CO2":
                return b.pureCO2_flow_mol_comp[t, "CO2"] == \
                    b.inlet_flow_mol_comp[t, "CO2"] * b.CO2_capture_rate[t]
            else:
                return b.pureCO2_flow_mol_comp[t, c] == 0.0

        # water drop in flue gas
        @self.Expression(self.flowsheet().config.time, doc="water drop")
        def water_drop(b, t):
            return b.inlet_flow_mol_comp[0, 'H2O']*0.5

        # Overall mass balances
        @self.Constraint(self.flowsheet().config.time,
                         self.component_list,
                         doc="Inlet component mole flow eqn")
        def flow_mol_comp_inlet_eqn(b, t, c):
            if c == "H2O":
                return b.inlet_flow_mol_comp[t, c] == \
                    b.exhaust_gas_flow_mol_comp[t, c] + b.water_drop[t]
            elif c == "CO2":
                return b.inlet_flow_mol_comp[t, c] == \
                    b.exhaust_gas_flow_mol_comp[t, c] \
                    + b.pureCO2_flow_mol_comp[t, "CO2"]
            else:
                return b.inlet_flow_mol_comp[t, c] == \
                    b.exhaust_gas_flow_mol_comp[t, c]

        # Pressure equations
        @self.Constraint(self.flowsheet().config.time,
                         doc="Pressure drop")
        def exh_pressure_eqn(b, t):
            return b.inlet_pressure[t] == b.exhaust_gas_pressure[t]

        @self.Constraint(self.flowsheet().config.time,
                         doc="Pressure drop")
        def pureCO2_pressure_eqn(b, t):
            return b.inlet_pressure[t] == b.pureCO2_pressure[t]

        # Temperature equations
        @self.Constraint(self.flowsheet().config.time,
                         doc="Temperature")
        def pureCO2_temp_eqn(b, t):
            return b.inlet_temperature[t] == b.pureCO2_temperature[t]

        @self.Constraint(self.flowsheet().config.time,
                         doc="Temperature")
        def exh_temp_eqn(b, t):
            return b.inlet_temperature[t] == b.exhaust_gas_temperature[t]

        # surrogates and additional constraints/expressions
        ''' This section is to add the surrogate models'''

        @self.Constraint(self.flowsheet().config.time,
                         doc="Specific reboiler duty in GJ/t CO2 or MJ/ kg CO2")
        def SRD_eqn(b, t):
            x1 = b.lean_loading[t]
            x2 = b.Pz_mol[t]
            return b.SRD[t] == (- 15565.227507172588957474 * x1
                                + 8.9237367058697429911263 * log(x1)
                                + 15285.362653727905126289 * exp(x1)
                                - 6741.0217492965430210461 * x1**2
                                - 3865.5900513667870654899 * x1**3
                                + 0.17312045830047528404555E-002 * x2**3
                                - 1.2960483459315448317994 * x1*x2
                                - 15235.956076189686427824)

        @self.Expression(self.flowsheet().config.time,
                         doc="Reboiler duty in MW")
        def reboiler_duty(b, t):   # mol/s /1000 = kgmol/s * 44 kg/kgmol = kg / s
            return b.SRD[t] * (b.pureCO2_flow_mol_comp[t, "CO2"]
                               * 44.01 / 1000)

        @self.Expression(self.flowsheet().config.time,
                         doc="Lean loading flowrate - kmol.hr")
        def LL_flowrate(b, t):
            x1 = b.lean_loading[t]
            x2 = b.Pz_mol[t]
            return (- 4402167.2503122268244624 * x1
                    - 4.4210697568168892956919 * x2
                    + 42394.191672889020992443 * log(x1)
                    + 4021949.9049337119795382 * exp(x1)
                    - 1394566.6646474981680512 * x1**2
                    - 1225122.9083037786185741 * x1**3
                    - 3898754.8073137737810612)


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
        iscale.calculate_scaling_factors(blk)  # remove to solve using baron
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")
        opt = pyo.SolverFactory(solver)
        opt.options = optarg

        init_log.info_low("Starting initialization...")

        blk.inlet.flow_mol_comp[0, 'CO2'].fix()
        blk.inlet.flow_mol_comp[0, 'O2'].fix()
        blk.inlet.flow_mol_comp[0, 'Ar'].fix()
        blk.inlet.flow_mol_comp[0, 'H2O'].fix()
        blk.inlet.flow_mol_comp[0, 'N2'].fix()
        blk.CO2_capture_rate.fix()

        # solve model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
                "Initialization Step 1 {}.".format(idaeslog.condition(res))
            )
        init_log.info_high("Initialization Step 1 Complete.")

        # ToDo: release state
        init_log.info("Initialization Complete.")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'CO2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'H2O'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'O2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'N2'], 1e-3)
        iscale.set_scaling_factor(self.inlet_flow_mol_comp[0.0, 'Ar'], 1e-3)
        iscale.set_scaling_factor(self.pureCO2_flow_mol_comp[0.0, 'CO2'], 1e-3)
        iscale.set_scaling_factor(self.pureCO2_flow_mol_comp[0.0, 'H2O'], 1e-3)
        iscale.set_scaling_factor(self.pureCO2_flow_mol_comp[0.0, 'O2'], 1e-3)
        iscale.set_scaling_factor(self.pureCO2_flow_mol_comp[0.0, 'N2'], 1e-3)
        iscale.set_scaling_factor(self.pureCO2_flow_mol_comp[0.0, 'Ar'], 1e-3)
        iscale.set_scaling_factor(self.exhaust_gas_flow_mol_comp[0.0, 'CO2'],
                                  1e-3)
        iscale.set_scaling_factor(self.exhaust_gas_flow_mol_comp[0.0, 'H2O'],
                                  1e-3)
        iscale.set_scaling_factor(self.exhaust_gas_flow_mol_comp[0.0, 'O2'],
                                  1e-3)
        iscale.set_scaling_factor(self.exhaust_gas_flow_mol_comp[0.0, 'N2'],
                                  1e-3)
        iscale.set_scaling_factor(self.exhaust_gas_flow_mol_comp[0.0, 'Ar'],
                                  1e-3)
        iscale.set_scaling_factor(self.inlet_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.pureCO2_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.exhaust_gas_temperature[0.0], 1e-2)
        iscale.set_scaling_factor(self.inlet_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.pureCO2_pressure[0.0], 1e-5)
        iscale.set_scaling_factor(self.exhaust_gas_pressure[0.0], 1e-5)

        for t, c in self.exh_pressure_eqn.items():
            sf = iscale.get_scaling_factor(
                self.inlet_pressure[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.pureCO2_pressure_eqn.items():
            sf = iscale.get_scaling_factor(
                self.inlet_pressure[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.pureCO2_temp_eqn.items():
            sf = iscale.get_scaling_factor(
                self.inlet_temperature[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.exh_temp_eqn.items():
            sf = iscale.get_scaling_factor(
                self.inlet_temperature[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.flow_mol_comp_inlet_eqn.items():
            sf = iscale.get_scaling_factor(
                self.inlet_flow_mol_comp[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf)

        for t, c in self.pureCO2_eqn.items():
            sf = iscale.get_scaling_factor(
                self.inlet_flow_mol_comp[t], default=1, warning=True)
            iscale.constraint_scaling_transform(c, sf)
