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
Basic property package for flue gas.

Main assumptions:
    - ideal gas
    - components in flue gas: O2, N2, NO, CO2, H2O, SO2
"""
# Import Pyomo libraries
from pyomo.environ import (Constraint, Param, PositiveReals, Reals,
                           value, log, exp, sqrt, Var, Expression, Reference)
from pyomo.opt import SolverFactory, TerminationCondition

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    Component,
    VaporPhase
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom, number_activated_constraints)
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis
from idaes.core.util.initialization import fix_state_vars, revert_state_vars
from idaes.core.util import constants
import idaes.core.util.scaling as iscale

# Import Python libraries
import idaes.logger as idaeslog


# Some more inforation about this module
__author__ = "Boiler Subsystem Team  J. Ma, M. Zamarripa, T. Burgard"
__version__ = "3"

# Set up logger
_log = idaeslog.getLogger('idaes.unit_model.properties')


@declare_process_block_class("FlueGasParameterBlock")
class FlueGasParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for flue
    gas. The ideal gas assumption is applied.

    """

    def build(self):
        """ Add contents to the block."""
        super(FlueGasParameterData, self).build()
        self._state_block_class = FlueGasStateBlock

        # Create Component objects
        self.N2 = Component()
        self.O2 = Component()
        self.NO = Component()
        self.CO2 = Component()
        self.H2O = Component()
        self.SO2 = Component()

        # Create Phase object
        self.Vap = VaporPhase()

        # Molecular weight
        self.mw_comp = Param(self.component_list,
                        initialize={'O2': 0.031998,
                                    'N2': 0.0280134,
                                    'NO': 0.0300057,
                                    'CO2': 0.04401,
                                    'H2O': 0.01801528,
                                    'SO2': 0.0640588},
                        doc='Molecular Weight [kg/mol]')

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=1.01325e5,
                                  doc='Reference pressure [Pa]')

        self.temperature_ref = Param(within=PositiveReals,
                                     mutable=True,
                                     default=298.16,
                                     doc='Reference temperature [K]')

        # Critical Properties
        self.pressure_crit = Param(self.component_list,
                                   within=PositiveReals,
                                   initialize={'O2': 50.45985e5,
                                               'N2': 33.943875e5,
                                               'NO': 64.85e5,
                                               'CO2': 73.8e5,
                                               'H2O': 220.64e5,
                                               'SO2': 7.883e6},
                                   doc='Critical pressure [Pa]')

        self.temperature_crit = Param(self.component_list,
                                      within=PositiveReals,
                                      initialize={'O2': 154.58,
                                                  'N2': 126.19,
                                                  'NO': 180.0,
                                                  'CO2': 304.18,
                                                  'H2O': 647,
                                                  'SO2': 430.8},
                                      doc='Critical temperature [K]')

        # Constants for specific heat capacity, enthalpy, and entropy
        # calculations for ideal gas (from NIST 01/08/2020
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=1#Thermo-Gas)
        cp_mol_ig_comp_coeff_parameter_table = {
            ('A', 'N2'): 19.50583,
            ('B', 'N2'): 19.88705,
            ('C', 'N2'): -8.598535,
            ('D', 'N2'): 1.369784,
            ('E', 'N2'): 0.527601,
            ('F', 'N2'): -4.935202,
            ('G', 'N2'): 212.39,
            ('H', 'N2'): 0,
            ('A', 'O2'): 30.03235,
            ('B', 'O2'): 8.772972,
            ('C', 'O2'): -3.98813,
            ('D', 'O2'): 0.788313,
            ('E', 'O2'): -0.7416,
            ('F', 'O2'): -11.3247,
            ('G', 'O2'): 236.1663,
            ('H', 'O2'): 0,
            ('A', 'CO2'): 24.99735,
            ('B', 'CO2'): 55.18696,
            ('C', 'CO2'): -33.69137,
            ('D', 'CO2'): 7.948387,
            ('E', 'CO2'): -0.136638,
            ('F', 'CO2'): -403.6075,
            ('G', 'CO2'): 228.2431,
            ('H', 'CO2'): -393.5224,
            ('A', 'H2O'): 30.092,
            ('B', 'H2O'): 6.832514,
            ('C', 'H2O'): 6.793435,
            ('D', 'H2O'): -2.53448,
            ('E', 'H2O'): 0.082139,
            ('F', 'H2O'): -250.881,
            ('G', 'H2O'): 223.3967,
            ('H', 'H2O'): -241.8264,
            ('A', 'NO'): 23.83491,
            ('B', 'NO'): 12.58878,
            ('C', 'NO'): -1.139011,
            ('D', 'NO'): -1.497459,
            ('E', 'NO'): 0.214194,
            ('F', 'NO'): 83.35783,
            ('G', 'NO'): 237.1219,
            ('H', 'NO'): 90.29114,
            ('A', 'SO2'): 21.43049,
            ('B', 'SO2'): 74.35094,
            ('C', 'SO2'): -57.75217,
            ('D', 'SO2'): 16.35534,
            ('E', 'SO2'): 0.086731,
            ('F', 'SO2'): -305.7688,
            ('G', 'SO2'): 254.8872,
            ('H', 'SO2'): -296.8422}

        # cp_mol_ig_comp_coeff units: J/(gmol-K)
        self.cp_mol_ig_comp_coeff = Param(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                           self.component_list,
                           initialize=cp_mol_ig_comp_coeff_parameter_table,
                           doc='Constants for spec. heat capacity for ideal '
                           'gas J/(gmol-K)')

        # Vapor pressure coefficients, currently uses N2 data for
        # NO since the data for NO are not available NIST Webbook 01/08/2020
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Units=SI&Mask=4#Thermo-Phase
        iv_pvap = {
            ('N2', 'A'): 3.7362,
            ('N2', 'B'): 264.651,
            ('N2', 'C'): -6.788,
            ('O2', 'A'): 3.9523,
            ('O2', 'B'): 340.024,
            ('O2', 'C'): -4.144,
            ('H2O', 'A'): 3.55959,
            ('H2O', 'B'): 643.748,
            ('H2O', 'C'): -198.043,
            ('CO2', 'A'): 6.81228,
            ('CO2', 'B'): 1301.679,
            ('CO2', 'C'): -3.494,
            ('NO', 'A'): 3.7362,
            ('NO', 'B'): 264.651,
            ('NO', 'C'): -6.788,
            ('SO2', 'A'): 4.37798,
            ('SO2', 'B'): 966.575,
            ('SO2', 'C'): -42.071}

        self.vapor_pressure_coeff = Param(
            self.component_list,
            ['A', 'B', 'C'],
            initialize=iv_pvap,
            mutable=True,
            doc="Antoine coefficients for vapor pressure P in bar, T in K")

        self.set_default_scaling("flow_mol", 1e-4)
        self.set_default_scaling("flow_mass", 1e-3)
        self.set_default_scaling("flow_vol", 1e-3)
        # anything not explicitly listed
        self.set_default_scaling("mole_frac_comp", 1)
        self.set_default_scaling("mole_frac_comp", 1e3, index="NO")
        self.set_default_scaling("mole_frac_comp", 1e3, index="SO2")
        self.set_default_scaling("mole_frac_comp", 1e2, index="H2O")
        self.set_default_scaling("mole_frac_comp", 1e2, index="CO2")
        self.set_default_scaling("flow_vol", 1)

        # For flow_mol_comp, will calculate from flow_mol and mole_frac_comp
        # user should set a scale for both, and for each compoent of mole_frac_comp
        self.set_default_scaling("pressure", 1e-5)
        self.set_default_scaling("temperature", 1e-1)
        self.set_default_scaling("pressure_red", 1e-3)
        self.set_default_scaling("temperature_red", 1)
        self.set_default_scaling("enth_mol_phase", 1e-3)
        self.set_default_scaling("enth_mol", 1e-3)
        self.set_default_scaling("entr_mol", 1e-2)
        self.set_default_scaling("entr_mol_phase", 1e-2)
        self.set_default_scaling("cp_mol", 0.1)
        self.set_default_scaling("cp_mol_phase", 0.1)
        self.set_default_scaling("compress_fact", 1)
        self.set_default_scaling("dens_mol_phase", 1)
        self.set_default_scaling("pressure_sat", 1e-4)
        self.set_default_scaling("visc_d_comp", 1e4)
        self.set_default_scaling("therm_cond_comp", 1e2)
        self.set_default_scaling("visc_d", 1e4)
        self.set_default_scaling("therm_cond", 1e2)
        self.set_default_scaling("mw", 1)
        self.set_default_scaling("mw_comp", 1)

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
            'flow_mol_comp': {'method': None, 'units': 'mol/s'},
            'pressure': {'method': None, 'units': 'Pa'},
            'temperature': {'method': None, 'units': 'K'},
            'pressure_crit': {'method': None, 'units': 'Pa'},
            'temperature_crit': {'method': None, 'units': 'K'},
            'pressure_red': {'method': None, 'units': None},
            'temperature_red': {'method': None, 'units': None},
            'enth_mol_phase': {'method': '_enthalpy_calc', 'units': 'J/mol'},
            'entr_mol_phase': {'method': '_entropy_calc', 'units': 'J/mol/K'},
            'enth_mol': {'method': '_enthalpy_calc', 'units': 'J/mol'},
            'entr_mol': {'method': '_entropy_calc', 'units': 'J/mol.K'},
            'cp_mol': {'method': '_heat_cap_calc', 'units': 'J/mol.K'},
            'cp_mol_phase': {'method': '_heat_cap_calc', 'units': 'J/mol/K'},
            'compress_fact': {'method': '_compress_fact', 'units': None},
            'dens_mol_phase': {'method': '_dens_mol_phase', 'units': 'mol/m^3'},
            'pressure_sat': {'method': '_vapor_pressure', 'units': 'Pa'},
            'flow_vol': {'method': '_flow_volume', 'units': 'm^3/s'},
            'visc_d': {'method': '_therm_cond', 'units': 'kg/m-s'},
            'therm_cond': {'method': '_therm_cond', 'units': 'W/m-K'},
            'mw_comp': {'method': None, 'units': 'kg/mol'},
            'mw': {'method': None, 'units': 'kg/mol'},
        })

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'kg',
                               'amount': 'mol',
                               'temperature': 'K'})


class _FlueGasStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(
        self,
        state_args={
            "flow_mol_comp": {
                "N2": 1.0,
                "CO2": 1.0,
                "NO": 1.0,
                "O2": 1.0,
                "H2O": 1.0,
                "SO2": 1.0
            },
            "pressure": 1e5,
            "temperature": 495.0},
        hold_state=False,
        state_vars_fixed=False,
        outlvl=0,
        solver='ipopt',
        optarg={'tol': 1e-8}
    ):
        """Initialisation routine for property package.

        Key values for the state_args dict:
            flow_mol_comp : value at which to initialize component flows
                (default=27.5e3 mol/s)
            pressure : value at which to initialize pressure (default=2.97e7 Pa)
            temperature : value at which to initialize temperature
                (default=866.5 K)

        Args:
            outlvl: sets logging level
            state_vars_fixed: Flag to denote state vars have been fixed.
                - True - states have been fixed by the control volume 1D.
                         Control volume 0D does not fix the state vars, so will
                         be False if this state block is used with 0D blocks.
                - False - states have not been fixed. The state block will deal
                          with fixing/unfixing.
            optarg: solver options dictionary object (default=None)
            solver: str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state: flag indicating whether the initialization routine
                should unfix any state variables fixed during initialization
                (default=False).
                - True - states varaibles are not unfixed, and a dict of
                         returned containing flags for which states were fixed
                         during initialization.
                - False - state variables are unfixed after initialization by
                          calling the relase_state method

            Returns:
                If hold_states is True, returns a dict containing flags for
                which states were fixed during initialization.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="properties")

        opt = SolverFactory(solver)
        opt.options = optarg

        if state_vars_fixed is False:
            flags = fix_state_vars(self, state_args)
        # Check when the state vars are fixed already result in dof 0
        for b in self.values():
            if degrees_of_freedom(b) != 0:
                raise Exception(f"{self.name} initializtion error: State vars "
                                "fixed but degrees of freedom not equal to 0")
        # ---------------------------------------------------------------------
        # Solve 1st stage
        for k, b in self.items():
            if b.is_property_constructed("pressure_sat"):
                b.pressure_sat.value = value(exp(
                    b.vapor_pressure_coeff[1]
                    + b.vapor_pressure_coeff[2] / b.temperature
                    + b.vapor_pressure_coeff[3] * blk.temperature
                    + b.vapor_pressure_coeff[4] * log(blk.temperature)
                    + b.vapor_pressure_coeff[5].value * blk.temperature**2))

            deactivate_list = []
            if hasattr(b, 'enthalpy_correlation'):
                deactivate_list.append(b.enthalpy_correlation)
            if hasattr(b, "volumetric_flow_calculation"):
                deactivate_list.append(b.volumetric_flow_calculation)
            if hasattr(b, "entropy_correlation"):
                deactivate_list.append(b.entropy_correlation)
            for c in deactivate_list:
                c.deactivate()

            if number_activated_constraints(b) > 0:
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(b, tee=slc.tee)
            else:
                res = "skipped"
            init_log.info_high(
                "Initialization Step 1 {}.".format(idaeslog.condition(res)))

            for c in deactivate_list:
                c.activate()

            if number_activated_constraints(b) > 0:
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(b, tee=slc.tee)
            else:
                res = "skipped"
            init_log.info_high(
                "Initialization Step 2 {}.".format(idaeslog.condition(res)))

        init_log.info(
            'Initialisation Complete, {}.'.format(
                idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        # If input block, return flags, else release state
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                self.release_state(flags)

    def release_state(self, flags, outlvl=idaeslog.NOTSET):
        '''
        Method to relase state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        # Unfix state variables
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="properties")
        revert_state_vars(self, flags)
        init_log.info('{} State Released.'.format(self.name))


@declare_process_block_class("FlueGasStateBlock",
                             block_class=_FlueGasStateBlock)
class FlueGasStateBlockData(StateBlockData):
    """
    This is an example of a property package for calculating the thermophysical
    properties of flue gas using the ideal gas assumption.
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(FlueGasStateBlockData, self).build()
        comps = self.params.component_list
        # Add state variables
        self.flow_mol_comp = Var(
            comps,
            domain=Reals,
            initialize=1.0,
            bounds=(0, 1e6),
            doc='Component molar flowrate [mol/s]'
        )
        self.pressure = Var(
            domain=Reals,
            initialize=1.01325e5,
            bounds=(1, 5e7),
            doc='State pressure [Pa]'
        )
        self.temperature = Var(
            domain=Reals,
            initialize=500,
            bounds=(200, 1500),
            doc='State temperature [K]'
        )

        # Add expressions for some basic oft-used quantiies
        self.flow_mol = Expression(
            expr=sum(self.flow_mol_comp[j] for j in comps))

        def rule_mole_frac(b, c):
            return b.flow_mol_comp[c] / b.flow_mol
        self.mole_frac_comp = Expression(
            comps,
            rule=rule_mole_frac,
            doc='mole fraction of component i'
        )

        self.flow_mass = Expression(
            expr=sum(self.flow_mol_comp[j] * self.params.mw_comp[j] for j in comps),
            doc='total mass flow')

        def rule_mw_comp(b, j):
            return b.params.mw_comp[j]
        self.mw_comp = Expression(comps, rule=rule_mw_comp)

        def rule_mw(b):
            return sum(b.mw_comp[j] * b.mole_frac_comp[j] for j in comps)
        self.mw = Expression(rule=rule_mw)

        self.pressure_crit = Expression(
            expr=sum(
                self.params.pressure_crit[j] *
                self.mole_frac_comp[j] for j in comps))
        self.temperature_crit = Expression(
            expr=sum(
                self.params.temperature_crit[j] *
                self.mole_frac_comp[j] for j in comps))
        self.pressure_red = Expression(
            expr=self.pressure / self.pressure_crit)
        self.temperature_red = Expression(
            expr=self.temperature / self.temperature_crit)

        self.compress_fact = Expression(
            expr=1.0, doc='Vapor Compressibility Factor')

        def rule_dens_mol_phase(b, p):
            return b.pressure / b.compress_fact / \
                constants.Constants.gas_constant / b.temperature
        self.dens_mol_phase = Expression(
            self.params.phase_list,
            rule=rule_dens_mol_phase,
            doc='Molar Density')

        self.flow_vol = Expression(
            doc='Volumetric Flowrate',
            expr=self.flow_mol / self.dens_mol_phase["Vap"])

    def _heat_cap_calc(self):
        # heat capacity J/mol-K
        self.cp_mol = Var(initialize=1000, doc='heat capacity [J/mol-K]')

        def rule_cp_phase(b, p):
            # This property module only has one phase
            return self.cp_mol
        self.cp_mol_phase = Expression(self.params.phase_list, rule=rule_cp_phase)

        try:
            coeff = self.params.cp_mol_ig_comp_coeff
            ft = sum(self.flow_mol_comp[j] for j in self.params.component_list)
            t = self.temperature / 1000
            self.heat_cap_correlation = Constraint(expr=(
                self.cp_mol * ft ==
                sum(self.flow_mol_comp[j] * (
                    coeff['A', j] +
                    coeff['B', j] * t +
                    coeff['C', j] * t**2 +
                    coeff['D', j] * t**3 +
                    coeff['E', j] / t**2) for j in self.params.component_list)))
        except AttributeError:
            self.del_component(self.cp_mol)
            self.del_component(self.heat_cap_correlation)

    def _enthalpy_calc(self):
        self.enth_mol = Var(doc='Specific Enthalpy [J/mol]')

        def rule_enth_phase(b, p):
            # This property module only has one phase
            return self.enth_mol
        self.enth_mol_phase = Expression(
            self.params.phase_list, rule=rule_enth_phase)

        def enthalpy_correlation(b):
            coeff = self.params.cp_mol_ig_comp_coeff
            ft = sum(self.flow_mol_comp[j] for j in self.params.component_list)
            t = self.temperature / 1000
            kJ_to_J = 1000
            return self.enth_mol * ft == (
                sum(kJ_to_J * self.flow_mol_comp[j] * (
                    coeff['A', j] * t +
                    coeff['B', j] * t**2 / 2 +
                    coeff['C', j] * t**3 / 3 +
                    coeff['D', j] * t**4 / 4 -
                    coeff['E', j] / t +
                    coeff['F', j] -
                    coeff['H', j]) for j in self.params.component_list))
        try:
            self.enthalpy_correlation = Constraint(rule=enthalpy_correlation)
        except AttributeError:
            self.del_component(self.enth_mol_phase)
            self.del_component(self.enth_mol)
            self.del_component(self.enthalpy_correlation)

    def _entropy_calc(self):
        self.entr_mol = Var(doc='Specific Entropy [J/mol/K]')
        # Specific Entropy

        def rule_entr_phase(b, p):
            # This property module only has one phase
            return self.entr_mol
        self.entr_mol_phase = Expression(
            self.params.phase_list, rule=rule_entr_phase)

        def entropy_correlation(b):
            coeff = self.params.cp_mol_ig_comp_coeff
            ft = sum(self.flow_mol_comp[j] for j in self.params.component_list)
            t = self.temperature / 1000
            n = self.flow_mol_comp
            x = self.mole_frac_comp
            r_gas = constants.Constants.gas_constant
            return self.entr_mol * ft  == \
                sum(n[j] * (
                    coeff['A', j] * log(t) +
                    coeff['B', j] * t +
                    coeff['C', j] * t**2 / 2 +
                    coeff['D', j] * t**3 / 3 -
                    coeff['E', j] / t**2 / 2 +
                    coeff['G', j] +
                    r_gas * log(x[j])) for j in self.params.component_list)
        try:
            self.entropy_correlation = Constraint(rule=entropy_correlation)
        except AttributeError:
            self.del_component(self.entr_mol_phase)
            self.del_component(self.entropy_correlation)

    def _vapor_pressure(self):
        # Vapour Pressure
        self.pressure_sat = Var(initialize=101325, doc="Vapour pressure [Pa]")

        def vapor_pressure_correlation(b):
            return (log(b.pressure_sat) *
                    sum(b.flow_mol_comp[j] for j in b.params.component_list) ==
                    sum((b.params.vapor_pressure_coeff[j, 'A'] * b.temperature -
                         (b.params.vapor_pressure_coeff[j, 'B'] /
                          (b.temperature + b._param.vapor_pressure_coeff[j, 'C']))) *
                        b.flow_mol_comp[j]
                        for j in b.params.component_list))

        try:
            self.vapor_pressure_correlation = \
                Constraint(rule=vapor_pressure_correlation)
        except AttributeError:
            self.del_component(self.pressure_sat)
            self.del_component(self.vapor_pressure_correlation)

    def _therm_cond(self):
        comps = self.params.component_list
        self.therm_cond_comp = Var(
            comps, initialize=0.05, doc='thermal conductivity J/m-K-s')
        self.therm_cond = Var(
            initialize=0.05, doc='thermal conductivity of gas mixture J/m-K-s')
        self.visc_d_comp = Var(
            comps, initialize=2e-5, doc='dynamic viscocity of pure gas species')
        self.visc_d = Var(
            initialize=2e-5, doc='viscosity of gas mixture kg/m-s')
        self.sigma = Param(
            comps,
            initialize={
                'O2': 3.458,
                'N2': 3.621,
                'NO': 3.47,
                'CO2': 3.763,
                'H2O': 2.605,
                'SO2': 4.29},
            doc='collision diameter in Angstrom (10e-10 mts)'
        )
        self.ep_Kappa = Param(
            comps,
            initialize={
                'O2': 107.4,
                'N2': 97.53,
                'NO': 119.0,
                'CO2': 244.0,
                'H2O': 572.4,
                'SO2': 252.0},
            doc="characteristic energy of interaction between pair of molecules "
                "K = Boltzmann constant in Kelvin")
        try:
            def rule_therm_cond(b, c):
                return b.therm_cond_comp[c] == (
                    ((b.params.cp_mol_ig_comp_coeff['A', c] +
                      b.params.cp_mol_ig_comp_coeff['B', c] *
                      (b.temperature /
                       1000) +
                      b.params.cp_mol_ig_comp_coeff['C', c] *
                      (b.temperature /
                       1000)**2 +
                      b.params.cp_mol_ig_comp_coeff['D', c] *
                      (b.temperature /
                       1000)**3 +
                      b.params.cp_mol_ig_comp_coeff['E', c] /
                      (b.temperature /
                       1000)**2) /
                     b.params.mw_comp[c]) +
                    1.25 *
                    (constants.Constants.gas_constant /
                     b.params.mw_comp[c])) * b.visc_d_comp[c]
            self.therm_cond_con = Constraint(comps, rule=rule_therm_cond)

            def rule_theta(b, c):
                return b.temperature / b.ep_Kappa[c]
            self.theta = Expression(comps, rule=rule_theta)

            def rule_omega(b, c):
                return (1.5794145
                        + 0.00635771 * b.theta[c]
                        - 0.7314 * log(b.theta[c])
                        + 0.2417357 * log(b.theta[c])**2
                        - 0.0347045 * log(b.theta[c])**3)
            self.omega = Expression(comps, rule=rule_omega)

            # Pure gas viscocity
            def rule_visc_d(b, c):
                return (b.visc_d_comp[c] * b.sigma[c]**2 * b.omega[c] ==
                        2.6693e-6 * sqrt(b.params.mw_comp[c] * 1000 * b.temperature))
            self.visc_d_con = Constraint(comps, rule=rule_visc_d)

            # section to calculate viscosity of gas mixture
            def rule_phi(b, i, j):
                return (1 / 2.8284
                        * (1 + (b.params.mw_comp[i] / b.params.mw_comp[j]))**(-0.5)
                        * (1 + sqrt(b.visc_d_comp[i] / b.visc_d_comp[j])
                            * (b.params.mw_comp[j] / b.params.mw_comp[i])**0.25)**2)
            self.phi_ij = Expression(
                comps,
                comps,
                rule=rule_phi
            )

            # viscosity of Gas mixture kg/m-s
            def rule_visc_d_mix(b):
                return b.visc_d == sum(
                    (b.mole_frac_comp[i] * b.visc_d_comp[i])
                    / sum(b.mole_frac_comp[j] * b.phi_ij[i, j] for j in comps)
                    for i in comps)
            self.vis_d_mix_con = Constraint(rule=rule_visc_d_mix)

            # thermal conductivity of gas mixture in kg/m-s
            def rule_therm_mix(b):
                return b.therm_cond == sum(
                    (b.mole_frac_comp[i] * b.therm_cond_comp[i])
                    / sum(b.mole_frac_comp[j] * b.phi_ij[i, j] for j in comps)
                    for i in comps)
            self.therm_mix_con = Constraint(rule=rule_therm_mix)

        except AttributeError:
            self.del_component(self.therm_cond_comp)
            self.del_component(self.therm_cond)
            self.del_component(self.visc_d_comp)
            self.del_component(self.visc_d)
            self.del_component(self.omega)
            self.del_component(self.theta)
            self.del_component(self.phi_ij)
            self.del_component(self.sigma)
            self.del_component(self.ep_Kappa)
            self.del_component(self.therm_cond_con)
            self.del_component(self.theta_con)
            self.del_component(self.omega_con)
            self.del_component(self.visc_d_con)
            self.del_component(self.phi_con)

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_terms(self, p, j):
        return self.flow_mol_comp[j]

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def get_enthalpy_flow_terms(self, p):
        if not self.is_property_constructed("enthalpy_flow_terms"):
            try:
                def rule_enthalpy_flow_terms(b, p):
                    return self.enth_mol_phase[p] * self.flow_mol
                self.enthalpy_flow_terms = Expression(
                    self.params.phase_list,
                    rule=rule_enthalpy_flow_terms
                )
            except AttributeError:
                self.del_component(enthalpy_flow_terms)
        return self.enthalpy_flow_terms[p]

    def get_material_density_terms(self, p, j):
        return self.dens_mol_phase[p]

    def get_energy_density_terms(self, p):
        if not self.is_property_constructed("energy_density_terms"):
            try:
                def rule_energy_density_terms(b, p):
                    return self.enth_mol_phase[p] * \
                        self.dens_mol_phase[p] - self.pressure
                self.energy_density_terms = Expression(
                    self.params.phase_list,
                    rule=rule_energy_density_terms
                )
            except AttributeError:
                self.del_component(energy_density_terms)
        return self.energy_density_terms[p]

    def define_state_vars(self):
        return {
            "flow_mol_comp": self.flow_mol_comp,
            "temperature": self.temperature,
            "pressure": self.pressure
        }

    def model_check(self):
        """
        Model checks for property block
        """
        # Check temperature bounds
        for v in self.compoent_object_data(Var, descend_into=True):
            if value(v) < v.lb:
                _log_error(f"{v} is below lower bound in {self.name}")
            if value(v) > v.ub:
                _log_error(f"{v} is above upper bound in {self.name}")

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # Get some scale factors that are frequently used to calculate others
        sf_flow = iscale.get_scaling_factor(self.flow_mol)
        sf_mol_fraction = {}
        comps = self.params.component_list
        for i in comps:
            sf_mol_fraction[i] = iscale.get_scaling_factor(self.mole_frac_comp[i])
        # calculate flow_mol_comp scale factors
        for i, c in self.flow_mol_comp.items():
            iscale.set_scaling_factor(c, sf_flow * sf_mol_fraction[i])

        if self.is_property_constructed("energy_density_terms"):
            for i, c in self.energy_density_terms.items():
                sf1 = iscale.get_scaling_factor(self.enth_mol_phase[i])
                sf2 = iscale.get_scaling_factor(self.dens_mol_phase[i])
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("enthalpy_flow_terms"):
            for i, c in self.enthalpy_flow_terms.items():
                sf1 = iscale.get_scaling_factor(self.enth_mol_phase[i])
                sf2 = iscale.get_scaling_factor(self.flow_mol)
                iscale.set_scaling_factor(c, sf1 * sf2)

        if self.is_property_constructed("heat_cap_correlation"):
            iscale.constraint_scaling_transform(
                self.heat_cap_correlation,
                iscale.get_scaling_factor(self.cp_mol) *
                iscale.get_scaling_factor(self.flow_mol)
            )
        if self.is_property_constructed("enthalpy_correlation"):
            for p, c in self.enthalpy_correlation.items():
                iscale.constraint_scaling_transform(
                    c,
                    iscale.get_scaling_factor(self.enth_mol) *
                    iscale.get_scaling_factor(self.flow_mol)
                )
        if self.is_property_constructed("entropy_correlation"):
            iscale.constraint_scaling_transform(
                self.entropy_correlation,
                iscale.get_scaling_factor(self.entr_mol)*
                iscale.get_scaling_factor(self.flow_mol)
            )
        if self.is_property_constructed("vapor_pressure_correlation"):
            iscale.constraint_scaling_transform(
                self.vapor_pressure_correlation,
                log(iscale.get_scaling_factor(self.pressure_sat)) *
                iscale.get_scaling_factor(self.flow_mol)
            )
        if self.is_property_constructed("therm_cond_con"):
            for i, c in self.therm_cond_con.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.therm_cond_comp[i]))
        if self.is_property_constructed("therm_mix_con"):
            iscale.constraint_scaling_transform(
                self.therm_mix_con,
                iscale.get_scaling_factor(self.therm_cond))
        if self.is_property_constructed("visc_d_con"):
            for i, c in self.visc_d_con.items():
                iscale.constraint_scaling_transform(
                    c, iscale.get_scaling_factor(self.visc_d_comp[i]))
        if self.is_property_constructed("visc_d_mix_con"):
            iscale.constraint_scaling_transform(
                self.visc_d_mix_con,
                iscale.get_scaling_factor(self.visc_d))
