##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Example property package for the combustion of methane in air using
Gibbs energy minimisation.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, log, Param, \
                          PositiveReals, Reals, Set, value, Var
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.opt import SolverFactory, TerminationCondition

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom

# Some more inforation about this module
__author__ = "Andrew Lee, Jinliang Ma"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("MethaneCombustionParameterBlock")
class PhysicalParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(PhysicalParameterData, self).build()

        self.state_block_class = MethaneCombustionStateBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2', 'N2', 'O2', 'CH4',
                                              'CO', 'CO2', 'H2O', 'NH3'])

        # List of components in each phase (optional)
        self.phase_comp = {"Vap": self.component_list}

        # List of all chemical elements that constitute the chemical species
        self.element_list = Set(initialize=['H', 'N', 'O', 'C'])

        # Elemental composition of all species
        self.element_comp = {'H2': {'H': 2, 'N': 0, 'O': 0, 'C': 0},
                             'N2': {'H': 0, 'N': 2, 'O': 0, 'C': 0},
                             'O2': {'H': 0, 'N': 0, 'O': 2, 'C': 0},
                             'CH4': {'H': 4, 'N': 0, 'O': 0, 'C': 1},
                             'CO': {'H': 0, 'N': 0, 'O': 1, 'C': 1},
                             'CO2': {'H': 0, 'N': 0, 'O': 2, 'C': 1},
                             'H2O': {'H': 2, 'N': 0, 'O': 1, 'C': 0},
                             'NH3': {'H': 3, 'N': 1, 'O': 0, 'C': 0}}

        # Heat capacity parameters - Shomate eqn (from NIST webbook)
        # Parameter 9 is the enthalpy at 1500K (Tref) in J/mol
        cp_param_dict = {('H2', 1): 18.563083,
                         ('H2', 2): 12.257357,
                         ('H2', 3): -2.859786,
                         ('H2', 4): 0.268238,
                         ('H2', 5): 1.977990,
                         ('H2', 6): -1.147438,
                         ('H2', 7): 156.288133,
                         ('H2', 8): 0.0,
                         ('H2', 9): 36290.28,
                         ('N2', 1): 19.50583,
                         ('N2', 2): 19.88705,
                         ('N2', 3): -8.598535,
                         ('N2', 4): 1.369784,
                         ('N2', 5): 0.527601,
                         ('N2', 6): -4.935202,
                         ('N2', 7): 212.3900,
                         ('N2', 8): 0,
                         ('N2', 9): 38405.02,
                         ('O2', 1): 30.03235,
                         ('O2', 2): 8.772972,
                         ('O2', 3): -3.988133,
                         ('O2', 4): 0.788313,
                         ('O2', 5): -0.741599,
                         ('O2', 6): -11.32468,
                         ('O2', 7): 236.1663,
                         ('O2', 8): 0,
                         ('O2', 9): 40598.9,
                         ('CH4', 1): 85.81217,
                         ('CH4', 2): 11.26467,
                         ('CH4', 3): -2.114146,
                         ('CH4', 4): 0.138190,
                         ('CH4', 5): -26.42221,
                         ('CH4', 6): -153.5327,
                         ('CH4', 7): 224.4143,
                         ('CH4', 8): -74.87310,
                         ('CH4', 9): 78142.7,
                         ('CO', 1): 35.15070,
                         ('CO', 2): 1.300095,
                         ('CO', 3): -0.205921,
                         ('CO', 4): 0.013550,
                         ('CO', 5): -3.282780,
                         ('CO', 6): -127.8375,
                         ('CO', 7): 231.7120,
                         ('CO', 8): -110.5271,
                         ('CO', 9): 38852.26,
                         ('CO2', 1): 58.16639,
                         ('CO2', 2): 2.720074,
                         ('CO2', 3): -0.492289,
                         ('CO2', 4): 0.038844,
                         ('CO2', 5): -6.447293,
                         ('CO2', 6): -425.9186,
                         ('CO2', 7): 263.6125,
                         ('CO2', 8): -393.5224,
                         ('CO2', 9): 61707.0,
                         ('H2O', 1): 41.96426,
                         ('H2O', 2): 8.622053,
                         ('H2O', 3): -1.499780,
                         ('H2O', 4): 0.098119,
                         ('H2O', 5): -11.15764,
                         ('H2O', 6): -272.1797,
                         ('H2O', 7): 219.7809,
                         ('H2O', 8): -241.8264,
                         ('H2O', 9): 48150.13,
                         ('NH3', 1): 52.02427,
                         ('NH3', 2): 18.48801,
                         ('NH3', 3): -3.765128,
                         ('NH3', 4): 0.248541,
                         ('NH3', 5): -12.45799,
                         ('NH3', 6): -85.53895,
                         ('NH3', 7): 223.8022,
                         ('NH3', 8): -45.89806,
                         ('NH3', 9): 63578.64}
        self.cp_params = Param(self.component_list,
                               range(1, 10),
                               mutable=False,
                               initialize=cp_param_dict,
                               doc="Shomate equation heat capacity parameters")

        # Heat of Formation
        hf_dict = {"CH4": -74600,
                   "CO": -110530,
                   "CO2": -393520,
                   "H2": 0,
                   "H2O": -241830,
                   "N2": 0,
                   "NH3": -45900,
                   "O2": 0}
        self.enth_mol_form = Param(
                self.component_list,
                mutable=False,
                initialize=hf_dict,
                doc="Component molar heats of formation [J/mol]")

        # Gas Constant
        self.gas_const = Param(within=PositiveReals,
                               mutable=False,
                               default=8.314,
                               doc='Gas Constant [J/mol.K]')

        # Thermodynamic reference state
        self.pressure_ref = Param(within=PositiveReals,
                                  mutable=True,
                                  default=101325.0,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(within=PositiveReals,
                                     mutable=True,
                                     default=1500.0,
                                     doc='Reference temperature [K]')

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'flow_mol_comp': {'method': None, 'units': 'mol/s'},
                'pressure': {'method': None, 'units': 'Pa'},
                'temperature': {'method': None, 'units': 'K'},
                'mole_frac_comp': {'method': None, 'units': None},
                'cp_mol': {'method': '_cp_mol', 'units': 'J/mol.K'},
                'cp_mol_comp': {'method': '_cp_mol_comp',
                                'units': 'J/mol.K'},
                'dens_mol_phase': {'method': '_dens_mol_phase',
                                   'units': 'mol/m^3'},
                'energy_internal_mol': {'method': '_energy_internal_mol',
                                        'units': 'J/mol'},
                'energy_internal_mol_phase_comp': {
                        'method': '_energy_internal_mol_phase_comp',
                        'units': 'J/mol'},
                'enth_mol': {'method': '_enth_mol', 'units': 'J/mol'},
                'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp',
                                        'units': 'J/mol'},
                'entr_mol': {'method': '_entr_mol', 'units': 'J/mol'},
                'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp',
                                        'units': 'J/mol'},
                'flow_mol': {'method': '_flow_mol', 'units': 'mol/s'},
                'gibbs_mol': {'method': '_gibbs_mol', 'units': 'J/mol'},
                'gibbs_mol_phase_comp': {'method': '_gibbs_mol_phase_comp',
                                         'units': 'J/mol'}})
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _StateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """
    def initialize(blk, state_args=None, hold_state=False, outlvl=0,
                   state_vars_fixed=False, solver='ipopt',
                   optarg={'tol': 1e-8}):
        '''
        Initialisation routine for property package.

        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provied at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:

                         flow_mol_comp : value at which to initialize component
                                         flows (default=None)
                         pressure : value at which to initialize pressure
                                    (default=None)
                         temperature : value at which to initialize temperature
                                      (default=None)
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states varaibles are not unfixed, and
                                 a dict of returned containing flags for
                                 which states were fixed during
                                 initialization.
                        - False - state variables are unfixed after
                                 initialization by calling the
                                 relase_state method

        Returns:
            If hold_states is True, returns a dict containing flags for
            which states were fixed during initialization.
        '''
        if state_vars_fixed is False:
            # Fix state variables if not already fixed
            Fcflag = {}
            Pflag = {}
            Tflag = {}

            for k in blk.keys():
                for j in blk[k]._params.component_list:
                    if blk[k].flow_mol_comp[j].fixed is True:
                        Fcflag[k, j] = True
                    else:
                        Fcflag[k, j] = False
                        if state_args is None:
                            blk[k].flow_mol_comp[j].fix(1.0)
                        else:
                            blk[k].flow_mol_comp[j].fix(state_args["flow_mol_comp"][j])

                if blk[k].pressure.fixed is True:
                    Pflag[k] = True
                else:
                    Pflag[k] = False
                    if state_args is None:
                        blk[k].pressure.fix(101325.0)
                    else:
                        blk[k].pressure.fix(state_args["pressure"])

                if blk[k].temperature.fixed is True:
                    Tflag[k] = True
                else:
                    Tflag[k] = False
                    if state_args is None:
                        blk[k].temperature.fix(1500.0)
                    else:
                        blk[k].temperature.fix(state_args["temperature"])

                for j in blk[k]._params.component_list:
                    blk[k].mole_frac_comp[j] = \
                        (value(blk[k].flow_mol_comp[j]) /
                         sum(value(blk[k].flow_mol_comp[i])
                             for i in blk[k]._params.component_list))

            # If input block, return flags, else release state
            flags = {"Fcflag": Fcflag, "Pflag": Pflag,
                     "Tflag": Tflag}

        else:
            # Check when the state vars are fixed already result in dof 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of freedom "
                                    "for state block is not zero during "
                                    "initialization.")
        # Set solver options
        if outlvl > 1:
            stee = True
        else:
            stee = False

        opt = SolverFactory(solver)
        opt.options = optarg

        # ---------------------------------------------------------------------
        # Initialise values
        for k in blk.keys():
            for j in blk[k]._params.component_list:

                if hasattr(blk[k], "cp_shomate_eqn"):
                    calculate_variable_from_constraint(blk[k].cp_mol_comp[j],
                                                       blk[k].
                                                       cp_shomate_eqn[j])

                if hasattr(blk[k], "enthalpy_shomate_eqn"):
                    calculate_variable_from_constraint(
                            blk[k].enth_mol_phase_comp["Vap", j],
                            blk[k].enthalpy_shomate_eqn[j])

                if hasattr(blk[k], "entropy_shomate_eqn"):
                    calculate_variable_from_constraint(
                            blk[k].entr_mol_phase_comp["Vap", j],
                            blk[k].entropy_shomate_eqn[j])

                if hasattr(blk[k], "partial_gibbs_energy_eqn"):
                    calculate_variable_from_constraint(
                            blk[k].gibbs_mol_phase_comp["Vap", j],
                            blk[k].partial_gibbs_energy_eqn[j])

            if hasattr(blk[k], "ideal_gas"):
                calculate_variable_from_constraint(
                            blk[k].dens_mol_phase["Vap"],
                            blk[k].ideal_gas)

            if hasattr(blk[k], "mixture_heat_capacity_eqn"):
                calculate_variable_from_constraint(
                            blk[k].cp_mol,
                            blk[k].mixture_heat_capacity_eqn)

            if hasattr(blk[k], "mixture_enthalpy_eqn"):
                calculate_variable_from_constraint(
                            blk[k].enth_mol,
                            blk[k].mixture_enthalpy_eqn)

            if hasattr(blk[k], "mixture_entropy_eqn"):
                calculate_variable_from_constraint(
                            blk[k].entr_mol,
                            blk[k].mixture_entropy_eqn)

            if hasattr(blk[k], "total_flow_eqn"):
                calculate_variable_from_constraint(
                            blk[k].flow_mol,
                            blk[k].total_flow_eqn)

            if hasattr(blk[k], "mixture_gibbs_eqn"):
                calculate_variable_from_constraint(
                            blk[k].gibbs_mol,
                            blk[k].mixture_gibbs_eqn)

        results = solve_indexed_blocks(opt, blk, tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info('{} Initialisation Step 1 Complete.'
                          .format(blk.name))
            else:
                _log.warning('{} Initialisation Step 1 Failed.'
                             .format(blk.name))

        # ---------------------------------------------------------------------
        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} Initialisation Complete.'.format(blk.name))

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of of logging
        '''
        if flags is None:
            return

        # Unfix state variables
        for k in blk.keys():
            for j in blk[k]._params.component_list:
                if flags['Fcflag'][k, j] is False:
                    blk[k].flow_mol_comp[j].unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} State Released.'.format(blk.name))


@declare_process_block_class("MethaneCombustionStateBlock",
                             block_class=_StateBlock)
class MethaneCombustionStateBlockData(StateBlockData):
    """
    An example property package for ideal gas properties with Gibbs energy
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(MethaneCombustionStateBlockData, self).build()

        # Create state variables
        self.flow_mol_comp = Var(self._params.component_list,
                                 initialize=1.0,
                                 doc='Component molar flowrate [mol/s]')
        self.pressure = Var(domain=Reals,
                            initialize=101325.0,
                            bounds=(1e3, 1e6),
                            doc='State pressure [Pa]')
        self.temperature = Var(domain=Reals,
                               initialize=1500,
                               bounds=(1450, 3000),
                               doc='State temperature [K]')
        self.mole_frac_comp = Var(self._params.component_list,
                             domain=Reals,
                             initialize=0.0,
                             doc='State component mole fractions [-]')

        # Create standard constraints
        # Mole fractions
        def mole_frac_constraint(b, j):
            return b.flow_mol_comp[j] == (
                       b.mole_frac_comp[j] *
                       sum(b.flow_mol_comp[k]
                           for k in b._params.component_list))
        self.mole_frac_constraint = Constraint(self._params.component_list,
                                               rule=mole_frac_constraint)

    def _dens_mol_phase(self):
        # Molar density
        self.dens_mol_phase = Var(self._params.phase_list,
                                  doc="Molar density")

        def ideal_gas(b, p):
            return (b.dens_mol_phase[p]*b._params.gas_const*b.temperature ==
                    b.pressure)
        try:
            # Try to build constraint
            self.ideal_gas = Constraint(self._params.phase_list,
                                        rule=ideal_gas)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.dens_mol_phase)
            self.del_component(self.ideal_gas)
            raise

    def _cp_mol_comp(self):
        # Pure component vapour heat capacities
        self.cp_mol_comp = Var(self._params.component_list,
                               domain=Reals,
                               initialize=1.0,
                               doc="Pure component vapour heat capacities "
                               "[J/mol.K]")

        def pure_component_cp_mol(b, j):
            return b.cp_mol_comp[j] == (
                        b._params.cp_params[j, 1] +
                        b._params.cp_params[j, 2]*(b.temperature*1e-3) +
                        b._params.cp_params[j, 3]*(b.temperature*1e-3)**2 +
                        b._params.cp_params[j, 4]*(b.temperature*1e-3)**3 +
                        b._params.cp_params[j, 5]/(b.temperature*1e-3)**2)
        try:
            # Try to build constraint
            self.cp_shomate_eqn = Constraint(self._params.component_list,
                                             rule=pure_component_cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol_comp)
            self.del_component(self.cp_shomate_eqn)
            raise

    def _cp_mol(self):
        # Mixture heat capacities
        self.cp_mol = Var(domain=Reals,
                          initialize=1.0,
                          doc="Mixture heat capacity [J/mol.K]")

        def cp_mol(b):
            return b.cp_mol == sum(b.cp_mol_comp[j]*b.mole_frac_comp[j]
                                   for j in b._params.component_list)
        try:
            # Try to build constraint
            self.mixture_heat_capacity_eqn = Constraint(
                                            rule=cp_mol)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.cp_mol)
            self.del_component(self.mixture_heat_capacity_eqn)
            raise

    def _energy_internal_mol_phase_comp(self):
        # Pure component vapour internal energies
        self.energy_internal_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Pure component internal energies [J/mol]")

        def pure_comp_int_energy(b, j):
            return b.energy_internal_mol_phase_comp["Vap", j] == (1e3*(
                    b._params.cp_params[j, 1]*(b.temperature*1e-3) +
                    b._params.cp_params[j, 2]*(b.temperature*1e-3)**2/2 +
                    b._params.cp_params[j, 3]*(b.temperature*1e-3)**3/3 +
                    b._params.cp_params[j, 4]*(b.temperature*1e-3)**4/4 -
                    b._params.cp_params[j, 5]/(b.temperature*1e-3) +
                    b._params.cp_params[j, 6] -
                    b._params.cp_params[j, 8]) -
                    b._params.cp_params[j, 9] +
                    b._params.enth_mol_form[j] -
                    b._params.gas_const*(b.temperature -
                                         b._params.temperature_ref))
        try:
            # Try to build constraint
            self.internal_energy_shomate_eqn = Constraint(
                    self._params.component_list,
                    rule=pure_comp_int_energy)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.energy_internal_phase_mol_comp)
            self.del_component(self.internal_energy_shomate_eqn)
            raise

    def _energy_internal_mol(self):
        # Mixture molar internal energy
        self.energy_internal_mol = Var(
                domain=Reals,
                initialize=0.0,
                doc='Mixture specific internal energy [J/mol]')
        try:
            # Try to build constraint
            self.mixture_energy_internal_eqn = Constraint(expr=(
                        self.energy_internal_mol == sum(
                                self.mole_frac_comp[j] *
                                self.energy_internal_mol_phase_comp["Vap", j]
                                for j in self._params.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.energy_internal_mol)
            self.del_component(self.mixture_energy_internal_eqn)
            raise

    def _enth_mol_phase_comp(self):
        # Pure component vapour enthalpies
        self.enth_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Pure component enthalpies [J/mol]")

        def pure_comp_enthalpy(b, j):
            return b.enth_mol_phase_comp["Vap", j] == (1e3*(
                    b._params.cp_params[j, 1]*(b.temperature*1e-3) +
                    b._params.cp_params[j, 2]*(b.temperature*1e-3)**2/2 +
                    b._params.cp_params[j, 3]*(b.temperature*1e-3)**3/3 +
                    b._params.cp_params[j, 4]*(b.temperature*1e-3)**4/4 -
                    b._params.cp_params[j, 5]/(b.temperature*1e-3) +
                    b._params.cp_params[j, 6] -
                    b._params.cp_params[j, 8]) -
                    b._params.cp_params[j, 9] +
                    b._params.enth_mol_form[j])
        try:
            # Try to build constraint
            self.enthalpy_shomate_eqn = Constraint(self._params.component_list,
                                                   rule=pure_comp_enthalpy)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol_phase_comp)
            self.del_component(self.enthalpy_shomate_eqn)
            raise

    def _enth_mol(self):
        # Mixture molar enthalpy
        self.enth_mol = Var(domain=Reals,
                            initialize=0.0,
                            doc='Mixture specific enthalpy [J/mol]')
        try:
            # Try to build constraint
            self.mixture_enthalpy_eqn = Constraint(expr=(
                        self.enth_mol == sum(
                                self.mole_frac_comp[j] *
                                self.enth_mol_phase_comp["Vap", j]
                                for j in self._params.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.enth_mol)
            self.del_component(self.mixture_enthalpy_eqn)
            raise

    def _entr_mol_phase_comp(self):
        # Partial component entropies
        self.entr_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                domain=Reals,
                initialize=1.0,
                doc="Partial component entropies [J/mol.K]")

        def partial_comp_entropy(b, j):
            return b.entr_mol_phase_comp["Vap", j] == (
                    b._params.cp_params[j, 1]*log((b.temperature*1e-3)) +
                    b._params.cp_params[j, 2]*(b.temperature*1e-3) +
                    b._params.cp_params[j, 3]*(b.temperature*1e-3)**2/2 +
                    b._params.cp_params[j, 4]*(b.temperature*1e-3)**3/3 -
                    b._params.cp_params[j, 5]/(2*(b.temperature*1e-3)**2) +
                    b._params.cp_params[j, 7] -
                    b._params.gas_const*log(b.mole_frac_comp[j] *
                                            b.pressure/b._params.pressure_ref))
        try:
            # Try to build constraint
            self.entropy_shomate_eqn = Constraint(self._params.component_list,
                                                  rule=partial_comp_entropy)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.entr_mol_phase_comp)
            self.del_component(self.entropy_shomate_eqn)
            raise

    def _entr_mol(self):
        # Mixture molar entropy
        self.entr_mol = Var(domain=Reals,
                            initialize=0.0,
                            doc='Mixture specific entropy [J/mol.K]')
        try:
            # Try to build constraint
            self.mixture_entropy_eqn = Constraint(expr=(
                        self.entr_mol == sum(
                                self.mole_frac_comp[j] *
                                self.entr_mol_phase_comp["Vap", j]
                                for j in self._params.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.entr_mol)
            self.del_component(self.mixture_entropy_eqn)
            raise

    def _flow_mol(self):
        # Total molar material flow
        self.flow_mol = Var(domain=Reals,
                            initialize=1.0,
                            doc="Total mixture molar flowrate [mol/s]")

        try:
            # Try to build constraint
            self.total_flow_eqn = Constraint(expr=(
                        self.flow_mol == sum(
                                self.flow_mol_comp[j]
                                for j in self._params.component_list)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.flow_mol)
            self.del_component(self.total_flow_eqn)
            raise

    def _gibbs_mol_phase_comp(self):
        # Partial component Gibbs free energy
        self.gibbs_mol_phase_comp = Var(
                            self._params.phase_list,
                            self._params.component_list,
                            domain=Reals,
                            initialize=0.0,
                            doc="Partial component Gibbs free energy [J/mol]")

        # Assume constant cp_mol for simplicity and vapour phase only
        def comp_gibbs_energy_equation(b, j):
            return b.gibbs_mol_phase_comp["Vap", j] == (
                        b.enth_mol_phase_comp["Vap", j] -
                        b.temperature*b.entr_mol_phase_comp["Vap", j])
        try:
            # Try to build constraint
            self.partial_gibbs_energy_eqn = Constraint(
                                            self._params.component_list,
                                            rule=comp_gibbs_energy_equation)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.gibbs_mol_comp)
            self.del_component(self.partial_gibbs_energy_eqn)
            raise

    def _gibbs_mol(self):
        # Mixture Gibbs energy
        self.gibbs_mol = Var(domain=Reals,
                             initialize=0.0,
                             doc='Mixture Gibbs energy [J/mol]')
        try:
            # Try to build constraint
            self.mixture_gibbs_eqn = Constraint(expr=(
                        self.gibbs_mol == (self.enth_mol -
                                           self.temperature*self.entr_mol)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.gibbs_mol)
            self.del_component(self.mixture_gibbs_eqn)
            raise

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return b.flow_mol*b.enth_mol

    def get_material_density_terms(b, p, j):
        return b.dens_mol_phase[p]*b.mole_frac_comp[j]

    def get_energy_density_terms(b, p):
        return b.dens_mol_phase[p]*b.energy_internal_mol

    def default_material_balance_type(self):
        return MaterialBalanceType.elementTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(b):
        return {"flow_mol_comp": b.flow_mol_comp,
                "temperature": b.temperature,
                "pressure": b.pressure}

    def model_check(blk):
        """
        Model checks for property block
        """
        # Check temperature bounds
        if value(blk.temperature) < blk.temperature.lb:
            _log.error('{} Temperature set below lower bound.'
                       .format(blk.name))
        if value(blk.temperature) > blk.temperature.ub:
            _log.error('{} Temperature set above upper bound.'
                       .format(blk.name))

        # Check pressure bounds
        if value(blk.pressure) < blk.pressure.lb:
            _log.error('{} Pressure set below lower bound.'.format(blk.name))
        if value(blk.pressure) > blk.pressure.ub:
            _log.error('{} Pressure set above upper bound.'.format(blk.name))
