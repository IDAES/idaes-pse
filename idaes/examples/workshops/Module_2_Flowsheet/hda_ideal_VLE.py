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
Example ideal parameter block for the VLE calucations for a
Benzene-Toluene-o-Xylene system.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, log, NonNegativeReals,\
    Var, Set, Param, sqrt, log10
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        MaterialBalanceType,
                        EnergyBalanceType)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import degrees_of_freedom, \
                                             number_unfixed_variables
from idaes.core.util.misc import extract_data
from idaes.core.util.constants import Constants as const

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("HDAParameterBlock")
class HDAParameterData(PhysicalParameterBlock):
    CONFIG = PhysicalParameterBlock.CONFIG()

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(HDAParameterData, self).build()

        self.state_block_class = IdealStateBlock

        self.component_list = Set(initialize=['benzene',
                                              'toluene',
                                              'hydrogen',
                                              'methane'])

        self.phase_list = Set(initialize=['Liq', 'Vap'],
                              ordered=True)

        # List of components in each phase (optional)
        self.phase_comp = {"Liq": self.component_list,
                           "Vap": self.component_list}

        # List of phase equilibrium index
        self.phase_equilibrium_idx = Set(initialize=[1, 2, 3, 4, 5])

        self.phase_equilibrium_list = \
            {1: ["benzene", ("Vap", "Liq")],
             2: ["toluene", ("Vap", "Liq")],
             3: ["hydrogen", ("Vap", "Liq")],
             4: ["methane", ("Vap", "Liq")],
             5: ["diphenyl", ("Vap", "Liq")]}

        # Thermodynamic reference state
        self.pressure_ref = Param(mutable=True,
                                  default=101325,
                                  doc='Reference pressure [Pa]')
        self.temperature_ref = Param(mutable=True,
                                     default=298.15,
                                     doc='Reference temperature [K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        pressure_crit_data = {'benzene': 48.9e5,
                              'toluene': 41e5,
                              'hydrogen': 12.9e5,
                              'methane': 46e5,
                              'diphenyl': 38.5e5
                              }

        self.pressure_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(pressure_crit_data),
            doc='Critical pressure [Pa]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        temperature_crit_data = {'benzene': 562.2,
                                 'toluene': 591.8,
                                 'hydrogen': 33.0,
                                 'methane': 190.4,
                                 'diphenyl': 789
                                 }

        self.temperature_crit = Param(
            self.component_list,
            within=NonNegativeReals,
            mutable=False,
            initialize=extract_data(temperature_crit_data),
            doc='Critical temperature [K]')

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        mw_comp_data = {'benzene': 78.1136E-3,
                        'toluene': 92.1405E-3,
                        'hydrogen': 2.016e-3,
                        'methane': 16.043e-3,
                        'diphenyl': 154.212e-4}

        self.mw_comp = Param(self.component_list,
                             mutable=False,
                             initialize=extract_data(mw_comp_data),
                             doc="molecular weight Kg/mol")

        # Constants for liquid densities
        # Source: Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        dens_liq_data = {('benzene', '1'): 1.0162,
                         ('benzene', '2'): 0.2655,
                         ('benzene', '3'): 562.16,
                         ('benzene', '4'): 0.28212,
                         ('toluene', '1'): 0.8488,
                         ('toluene', '2'): 0.26655,
                         ('toluene', '3'): 591.8,
                         ('toluene', '4'): 0.2878,
                         ('hydrogen', '1'): 5.414,
                         ('hydrogen', '2'): 0.34893,
                         ('hydrogen', '3'): 33.19,
                         ('hydrogen', '4'): 0.2706,
                         ('methane', '1'): 2.9214,
                         ('methane', '2'): 0.28976,
                         ('methane', '3'): 190.56,
                         ('methane', '4'): 0.28881,
                         ('diphenyl', '1'): 0.5039,
                         ('diphenyl', '2'): 0.25273,
                         ('diphenyl', '3'): 789.26,
                         ('diphenyl', '4'): 0.281}

        self.dens_liq_params = Param(
                self.component_list,
                ['1', '2', '3', '4'],
                mutable=False,
                initialize=extract_data(dens_liq_data),
                doc="Parameters to compute liquid densities")

        # Boiling point at standard pressure
        # Source: Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        bp_data = {('benzene'): 353.25,
                   ('toluene'): 383.95,
                   ('hydrogen'): 20.45,
                   ('methane'): 111.75,
                   ('diphenyl'): 528.05}

        self.temperature_boil = Param(
                self.component_list,
                mutable=False,
                initialize=extract_data(bp_data),
                doc="Pure component boiling points at standard pressure [K]")

        # Constants for specific heat capacity, enthalpy
        # Sources: The Properties of Gases and Liquids (1987)
        #         4th edition, Chemical Engineering Series - Robert C. Reid
        #         Perry's Chemical Engineers Handbook
        #         - Robert H. Perry (Cp_liq)
        cp_ig_data = {('Liq', 'benzene', '1'): 1.29E5,
                      ('Liq', 'benzene', '2'): -1.7E2,
                      ('Liq', 'benzene', '3'): 6.48E-1,
                      ('Liq', 'benzene', '4'): 0,
                      ('Liq', 'benzene', '5'): 0,
                      ('Vap', 'benzene', '1'): -3.392E1,
                      ('Vap', 'benzene', '2'): 4.739E-1,
                      ('Vap', 'benzene', '3'): -3.017E-4,
                      ('Vap', 'benzene', '4'): 7.130E-8,
                      ('Vap', 'benzene', '5'): 0,
                      ('Liq', 'toluene', '1'): 1.40E5,
                      ('Liq', 'toluene', '2'): -1.52E2,
                      ('Liq', 'toluene', '3'): 6.95E-1,
                      ('Liq', 'toluene', '4'): 0,
                      ('Liq', 'toluene', '5'): 0,
                      ('Vap', 'toluene', '1'): -2.435E1,
                      ('Vap', 'toluene', '2'): 5.125E-1,
                      ('Vap', 'toluene', '3'): -2.765E-4,
                      ('Vap', 'toluene', '4'): 4.911E-8,
                      ('Vap', 'toluene', '5'): 0,
                      ('Liq', 'hydrogen', '1'): 0,  # 6.6653e1,
                      ('Liq', 'hydrogen', '2'): 0,  # 6.7659e3,
                      ('Liq', 'hydrogen', '3'): 0,  # -1.2363e2,
                      ('Liq', 'hydrogen', '4'): 0,  # 4.7827e2, # Eqn 2
                      ('Liq', 'hydrogen', '5'): 0,
                      ('Vap', 'hydrogen', '1'): 2.714e1,
                      ('Vap', 'hydrogen', '2'): 9.274e-3,
                      ('Vap', 'hydrogen', '3'): -1.381e-5,
                      ('Vap', 'hydrogen', '4'): 7.645e-9,
                      ('Vap', 'hydrogen', '5'): 0,
                      ('Liq', 'methane', '1'): 0,  # 6.5708e1,
                      ('Liq', 'methane', '2'): 0,  # 3.8883e4,
                      ('Liq', 'methane', '3'): 0,  # -2.5795e2,
                      ('Liq', 'methane', '4'): 0,  # 6.1407e2, # Eqn 2
                      ('Liq', 'methane', '5'): 0,
                      ('Vap', 'methane', '1'): 1.925e1,
                      ('Vap', 'methane', '2'): 5.213e-2,
                      ('Vap', 'methane', '3'): 1.197e-5,
                      ('Vap', 'methane', '4'): -1.132e-8,
                      ('Vap', 'methane', '5'): 0,
                      ('Liq', 'diphenyl', '1'): 1.2177e5,
                      ('Liq', 'diphenyl', '2'): 4.2930e2,
                      ('Liq', 'diphenyl', '3'): 0,
                      ('Liq', 'diphenyl', '4'): 0,
                      ('Liq', 'diphenyl', '5'): 0,
                      ('Vap', 'diphenyl', '1'): -9.707e1,
                      ('Vap', 'diphenyl', '2'): 1.106e0,
                      ('Vap', 'diphenyl', '3'): -8.855e-4,
                      ('Vap', 'diphenyl', '4'): 2.790e-7,
                      ('Vap', 'diphenyl', '5'): 0}

        self.cp_ig = Param(self.phase_list, self.component_list,
                           ['1', '2', '3', '4', '5'],
                           mutable=False,
                           initialize=extract_data(cp_ig_data),
                           doc="parameters to compute Cp_comp")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        # fitted to Antoine form
        # H2, Methane from NIST webbook
        pressure_sat_coeff_data = {('benzene', 'A'): 4.202,
                                   ('benzene', 'B'): 1322,
                                   ('benzene', 'C'): -38.56,
                                   ('toluene', 'A'): 4.216,
                                   ('toluene', 'B'): 1435,
                                   ('toluene', 'C'): -43.33,
                                   ('hydrogen', 'A'): 3.543,
                                   ('hydrogen', 'B'): 99.40,
                                   ('hydrogen', 'C'): 7.726,
                                   ('methane', 'A'): 3.990,
                                   ('methane', 'B'): 443.0,
                                   ('methane', 'C'): -0.49,
                                   ('diphenyl', 'A'): 4.345,
                                   ('diphenyl', 'B'): 1988,
                                   ('diphenyl', 'C'): -70.82}

        self.pressure_sat_coeff = Param(
            self.component_list,
            ['A', 'B', 'C'],
            mutable=False,
            initialize=extract_data(pressure_sat_coeff_data),
            doc="parameters to compute Cp_comp")

        # Source: The Properties of Gases and Liquids (1987)
        # 4th edition, Chemical Engineering Series - Robert C. Reid
        dh_vap = {'benzene': 3.387e4, 'toluene': 3.8262e4,
                  'hydrogen': 0, 'methane': 0, "diphenyl": 6.271e4}

        self.dh_vap = Param(self.component_list,
                            mutable=False,
                            initialize=extract_data(dh_vap),
                            doc="heat of vaporization")

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mol': {'method': None, 'units': 'mol/s'},
             'flow_mol_phase_comp': {'method': None, 'units': 'mol/s'},
             'mole_frac_comp': {'method': None, 'units': 'none'},
             'temperature': {'method': None, 'units': 'K'},
             'pressure': {'method': None, 'units': 'Pa'},
             'flow_mol_phase': {'method': None, 'units': 'mol/s'},
             'dens_mol_phase': {'method': '_dens_mol_phase',
                                'units': 'mol/m^3'},
             'pressure_sat': {'method': '_pressure_sat', 'units': 'Pa'},
             'mole_frac_phase_comp': {'method': '_mole_frac_phase',
                                      'units': 'no unit'},
             'energy_internal_mol_phase_comp': {
                     'method': '_energy_internal_mol_phase_comp',
                     'units': 'J/mol'},
             'energy_internal_mol_phase': {
                     'method': '_enenrgy_internal_mol_phase',
                     'units': 'J/mol'},
             'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp',
                                     'units': 'J/mol'},
             'enth_mol_phase': {'method': '_enth_mol_phase',
                                'units': 'J/mol'},
             'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp',
                                     'units': 'J/mol'},
             'entr_mol_phase': {'method': '_entr_mol_phase',
                                'units': 'J/mol'},
             'temperature_bubble': {'method': '_temperature_bubble',
                                    'units': 'K'},
             'temperature_dew': {'method': '_temperature_dew',
                                 'units': 'K'},
             'pressure_bubble': {'method': '_pressure_bubble',
                                 'units': 'Pa'},
             'pressure_dew': {'method': '_pressure_dew',
                              'units': 'Pa'},
             'fug_vap': {'method': '_fug_vap', 'units': 'Pa'},
             'fug_liq': {'method': '_fug_liq', 'units': 'Pa'},
             'dh_vap': {'method': '_dh_vap', 'units': 'J/mol'},
             'ds_vap': {'method': '_ds_vap', 'units': 'J/mol.K'}})

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _IdealStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(blk, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=1,
                   solver='ipopt', optarg={'tol': 1e-8}):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : Dictionary with initial guesses for the state vars
                         chosen. Note that if this method is triggered
                         through the control volume, and if initial guesses
                         were not provied at the unit model level, the
                         control volume passes the inlet values as initial
                         guess.The keys for the state_args dictionary are:

                         flow_mol_phase_comp : value at which to initialize
                                               phase component flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
            outlvl : sets output level of initialization routine
                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)
            optarg : solver options dictionary object (default=None)
            state_vars_fixed: Flag to denote if state vars have already been
                              fixed.
                              - True - states have already been fixed by the
                                       control volume 1D. Control volume 0D
                                       does not fix the state vars, so will
                                       be False if this state block is used
                                       with 0D blocks.
                             - False - states have not been fixed. The state
                                       block will deal with fixing/unfixing.
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
        """

        _log.info('Starting {} initialization'.format(blk.name))

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flags = fix_state_vars(blk, state_args)

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

        if optarg is None:
            sopt = {'tol': 1e-8}
        else:
            sopt = optarg

        opt = SolverFactory('ipopt')
        opt.options = sopt

        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations
        for k in blk.keys():
            if hasattr(blk[k], "eq_temperature_dew"):
                calculate_variable_from_constraint(blk[k].temperature_dew,
                                                   blk[k].eq_temperature_dew)

            if hasattr(blk[k], "eq_pressure_dew"):
                calculate_variable_from_constraint(blk[k].pressure_dew,
                                                   blk[k].eq_pressure_dew)

        if outlvl > 0:
            _log.info("Dew and bubble points initialization for "
                      "{} completed".format(blk.name))

        # ---------------------------------------------------------------------
        # If flash, initialize T1 and Teq
        for k in blk.keys():
            if (blk[k].config.has_phase_equilibrium and
                    not blk[k].config.defined_state):
                blk[k]._t1.value = max(blk[k].temperature.value,
                                       blk[k].temperature_bubble.value)
                blk[k]._teq.value = min(blk[k]._t1.value,
                                        blk[k].temperature_dew.value)

        if outlvl > 0:
            _log.info("Equilibrium temperature initialization for "
                      "{} completed".format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        # TODO : This will need ot be generalised more when we move to a
        # modular implementation
        for k in blk.keys():
            # Deactivate equilibrium constraints, as state is fixed
            if hasattr(blk[k], 'equilibrium_constraint'):
                blk[k].equilibrium_constraint.deactivate()

        free_vars = 0
        for k in blk.keys():
            free_vars += number_unfixed_variables(blk[k])
        if free_vars > 0:
            try:
                results = solve_indexed_blocks(opt, [blk], tee=stee)
            except:
                results = None
        else:
            results = None

        for k in blk.keys():
            # Reactivate equilibrium constraints
            if hasattr(blk[k], 'equilibrium_constraint'):
                blk[k].equilibrium_constraint.activate()

        if outlvl > 0:
            if results is None or results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Property initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Property initialization for "
                             "{} failed".format(blk.name))

        # ---------------------------------------------------------------------
        # Return state to initial conditions
        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        if outlvl > 0:
            _log.info("Initialization completed for {}".format(blk.name))

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialization.
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
        revert_state_vars(blk, flags)

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} states released.'.format(blk.name))


@declare_process_block_class("IdealStateBlock",
                             block_class=_IdealStateBlock)
class IdealStateBlockData(StateBlockData):
    """An example property package for ideal VLE."""

    def build(self):
        """Callable method for Block construction."""
        super(IdealStateBlockData, self).build()

        # Add state variables
        self.flow_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                initialize=0.5,
                bounds=(1e-8, 100),
                doc='Phase-component molar flow rates [mol/s]')

        self.pressure = Var(initialize=101325,
                            bounds=(101325, 400000),
                            domain=NonNegativeReals,
                            doc='State pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               bounds=(298.15, 1000),
                               domain=NonNegativeReals,
                               doc='State temperature [K]')

        # Add supporting variables
        def flow_mol_phase(b, p):
            return sum(b.flow_mol_phase_comp[p, j]
                       for j in b._params.component_list)
        self.flow_mol_phase = Expression(self._params.phase_list,
                                         rule=flow_mol_phase,
                                         doc='Phase molar flow rates [mol/s]')

        def flow_mol(b):
            return sum(b.flow_mol_phase_comp[p, j]
                       for j in b._params.component_list
                       for p in b._params.phase_list)
        self.flow_mol = Expression(rule=flow_mol,
                                   doc='Total molar flowrate [mol/s]')

        def mole_frac_phase_comp(b, p, j):
            return b.flow_mol_phase_comp[p, j]/b.flow_mol_phase[p]
        self.mole_frac_phase_comp = Expression(
                self._params.phase_list,
                self._params.component_list,
                rule=mole_frac_phase_comp,
                doc='Phase mole fractions [-]')

        def mole_frac_comp(b, j):
            return (sum(b.flow_mol_phase_comp[p, j]
                        for p in b._params.phase_list) / b.flow_mol)
        self.mole_frac_comp = Expression(self._params.component_list,
                                         rule=mole_frac_comp,
                                         doc='Mixture mole fractions [-]')

        # Reaction Stoichiometry
        add_object_reference(self, "phase_equilibrium_list_ref",
                             self._params.phase_equilibrium_list)

        if (self.config.has_phase_equilibrium and
                self.config.defined_state is False):
            # Definition of equilibrium temperature for smooth VLE
            self._teq = Var(
                    initialize=self.temperature.value,
                    doc='Temperature for calculating phase equilibrium')
            self._t1 = Var(initialize=self.temperature.value,
                           doc='Intermediate temperature for calculating Teq')

            self.eps_1 = Param(default=0.01,
                               mutable=True,
                               doc='Smoothing parameter for Teq')
            self.eps_2 = Param(default=0.0005,
                               mutable=True,
                               doc='Smoothing parameter for Teq')

            # PSE paper Eqn 13
            def rule_t1(b):
                return b._t1 == 0.5*(
                        b.temperature + b.temperature_bubble +
                        sqrt((b.temperature-b.temperature_bubble)**2 +
                             b.eps_1**2))
            self._t1_constraint = Constraint(rule=rule_t1)

            # PSE paper Eqn 14
            # TODO : Add option for supercritical extension
            def rule_teq(b):
                return b._teq == 0.5*(b._t1 + b.temperature_dew -
                                      sqrt((b._t1-b.temperature_dew)**2 +
                                           b.eps_2**2))
            self._teq_constraint = Constraint(rule=rule_teq)

            def rule_tr_eq(b, i):
                return b._teq / b._params.temperature_crit[i]
            self._tr_eq = Expression(
                    self._params.component_list,
                    rule=rule_tr_eq,
                    doc='Component reduced temperatures [-]')

            def rule_equilibrium(b, i):
                return b.fug_vap[i] == b.fug_liq[i]
            self.equilibrium_constraint = Constraint(
                    self._params.component_list, rule=rule_equilibrium)

# -----------------------------------------------------------------------------
# Property Methods
    def _dens_mol_phase(self):
        self.dens_mol_phase = Var(self._params.phase_list,
                                  initialize=1.0,
                                  doc="Molar density [mol/m^3]")

        def rule_dens_mol_phase(b, p):
            if p == 'Vap':
                return b._dens_mol_vap()
            else:
                return b._dens_mol_liq()
        self.eq_dens_mol_phase = Constraint(self._params.phase_list,
                                            rule=rule_dens_mol_phase)

    def _energy_internal_mol_phase_comp(self):
        self.energy_internal_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                doc="Phase-component molar specific internal energies [J/mol]")

        def rule_energy_internal_mol_phase_comp(b, p, j):
            if p == 'Vap':
                return b.energy_internal_mol_phase_comp[p, j] == \
                        b.enth_mol_phase_comp[p, j] - \
                        const.gas_constant*(b.temperature -
                                            b._params.temeprature_ref)
            else:
                return b.energy_internal_mol_phase_comp[p, j] == \
                        b.enth_mol_phase_comp[p, j]
        self.eq_energy_internal_mol_phase_comp = Constraint(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_energy_internal_mol_phase_comp)

    def _energy_internal_mol_phase(self):
        self.energy_internal_mol_phase = Var(
            self._params.phase_list,
            doc='Phase molar specific internal energies [J/mol]')

        def rule_energy_internal_mol_phase(b, p):
            return b.energy_internal_mol_phase[p] == sum(
                b.energy_internal_mol_phase_comp[p, i] *
                b.mole_frac_phase_comp[p, i]
                for i in b._params.component_list)
        self.eq_energy_internal_mol_phase = Constraint(
                self._params.phase_list,
                rule=rule_energy_internal_mol_phase)

    def _enth_mol_phase_comp(self):
        self.enth_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                initialize=7e5,
                doc='Phase-component molar specific enthalpies [J/mol]')

        def rule_enth_mol_phase_comp(b, p, j):
            if p == 'Vap':
                return b._enth_mol_comp_vap(j)
            else:
                return b._enth_mol_comp_liq(j)
        self.eq_enth_mol_phase_comp = Constraint(
                self._params.phase_list,
                self._params.component_list,
                rule=rule_enth_mol_phase_comp)

    def _enth_mol_phase(self):
        self.enth_mol_phase = Var(
                self._params.phase_list,
                initialize=7e5,
                doc='Phase molar specific enthalpies [J/mol]')

        def rule_enth_mol_phase(b, p):
            return b.enth_mol_phase[p] == sum(
                    b.enth_mol_phase_comp[p, i] *
                    b.mole_frac_phase_comp[p, i]
                    for i in b._params.component_list)
        self.eq_enth_mol_phase = Constraint(self._params.phase_list,
                                            rule=rule_enth_mol_phase)

    def _entr_mol_phase_comp(self):
        self.entr_mol_phase_comp = Var(
                self._params.phase_list,
                self._params.component_list,
                doc='Phase-component molar specific entropies [J/mol.K]')

        def rule_entr_mol_phase_comp(b, p, j):
            if p == 'Vap':
                return b._entr_mol_comp_vap(j)
            else:
                return b._entr_mol_comp_liq(j)
        self.eq_entr_mol_phase_comp = Constraint(
                self._params.phase_list,
                self._params.component_list,
                rule=rule_entr_mol_phase_comp)

    def _entr_mol_phase(self):
        self.entr_mol_phase = Var(
                self._params.phase_list,
                doc='Phase molar specific enthropies [J/mol.K]')

        def rule_entr_mol_phase(b, p):
            return b.entr_mol_phase[p] == sum(
                    b.entr_mol_phase_comp[p, i] *
                    b.mole_frac_phase_comp[p, i]
                    for i in b._params.component_list)
        self.eq_entr_mol_phase = Constraint(self._params.phase_list,
                                            rule=rule_entr_mol_phase)

# -----------------------------------------------------------------------------
# General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if j in self._params.component_list:
            return self.flow_mol_phase_comp[p, j]
        else:
            return 0

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.flow_mol_phase[p] * self.enth_mol_phase[p]

    def get_material_density_terms(self, p, j):
        """Create material density terms."""
        if j in self._params.component_list:
            return self.dens_mol_phase[p] * self.mole_frac_phase_comp[p, j]
        else:
            return 0

    def get_enthalpy_density_terms(self, p):
        """Create enthalpy density terms."""
        return self.dens_mol_phase[p] * self.energy_internal_mol_phase[p]

    def default_material_balance_type(self):
        return MaterialBalanceType.componentPhase

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_mol_phase_comp": self.flow_mol_phase_comp,
                "temperature": self.temperature,
                "pressure": self.pressure}

    # Property package utility functions
    def calculate_bubble_point_temperature(self, clear_components=True):
        """"To compute the bubble point temperature of the mixture."""

        if hasattr(self, "eq_temperature_bubble"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.temperature_bubble,
                                           self.eq_temperature_bubble)

        return self.temperature_bubble.value

        if clear_components is True:
            self.del_component(self.eq_temperature_bubble)
            self.del_component(self._p_sat_bubbleT)
            self.del_component(self.temperature_bubble)

    def calculate_dew_point_temperature(self, clear_components=True):
        """"To compute the dew point temperature of the mixture."""

        if hasattr(self, "eq_temperature_dew"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.temperature_dew,
                                           self.eq_temperature_dew)

        return self.temperature_dew.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_temperature_dew)
            self.del_component(self._p_sat_dewT)
            self.del_component(self.temperature_dew)

    def calculate_bubble_point_pressure(self, clear_components=True):
        """"To compute the bubble point pressure of the mixture."""

        if hasattr(self, "eq_pressure_bubble"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.pressure_bubble,
                                           self.eq_pressure_bubble)

        return self.pressure_bubble.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_pressure_bubble)
            self.del_component(self._p_sat_bubbleP)
            self.del_component(self.pressure_bubble)

    def calculate_dew_point_pressure(self, clear_components=True):
        """"To compute the dew point pressure of the mixture."""

        if hasattr(self, "eq_pressure_dew"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.pressure_dew,
                                           self.eq_pressure_dew)

        return self.pressure_dew.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_pressure_dew)
            self.del_component(self._p_sat_dewP)
            self.del_component(self.pressure_dew)

# -----------------------------------------------------------------------------
# Bubble and Dew Points
# Ideal-Ideal properties allow for the simplifications below
# Other methods require more complex equations with shadow compositions

# For future work, propose the following:
# Core class writes a set of constraints Phi_L_i == Phi_V_i
# Phi_L_i and Phi_V_i make calls to submethods which add shadow compositions
# as needed
    def _temperature_bubble(self):
        self.temperature_bubble = Param(initialize=33.0,
                                        doc="Bubble point temperature (K)")

    def _temperature_dew(self):

        self.temperature_dew = Var(initialize=298.15,
                                   doc="Dew point temperature (K)")

        def rule_psat_dew(b, j):
            return 1e5*10**(b._params.pressure_sat_coeff[j, 'A'] -
                            b._params.pressure_sat_coeff[j, 'B'] /
                            (b.temperature_dew +
                             b._params.pressure_sat_coeff[j, 'C']))

        try:
            # Try to build expression
            self._p_sat_dewT = Expression(self._params.component_list,
                                          rule=rule_psat_dew)

            def rule_temp_dew(b):
                return b.pressure * sum(b.mole_frac_comp[i] /
                                        b._p_sat_dewT[i]
                                        for i in ['toluene', 'benzene']) \
                    - 1 == 0
            self.eq_temperature_dew = Constraint(rule=rule_temp_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_dew)
            self.del_component(self._p_sat_dewT)

    def _pressure_bubble(self):
        self.pressure_bubble = Param(initialize=1e8,
                                     doc="Bubble point pressure (Pa)")

    def _pressure_dew(self):
        self.pressure_dew = Var(initialize=298.15,
                                doc="Dew point pressure (Pa)")

        def rule_psat_dew(b, j):
            return 1e5*10**(b._params.pressure_sat_coeff[j, 'A'] -
                            b._params.pressure_sat_coeff[j, 'B'] /
                            (b.temperature +
                             b._params.pressure_sat_coeff[j, 'C']))

        try:
            # Try to build expression
            self._p_sat_dewP = Expression(self._params.component_list,
                                          rule=rule_psat_dew)

            def rule_pressure_dew(b):
                return b.pressure_dew * \
                    sum(b.mole_frac_comp[i] / b._p_sat_dewP[i]
                        for i in ['toluene', 'benzene']) \
                    - 1 == 0
            self.eq_pressure_dew = Constraint(rule=rule_pressure_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.pressure_dew)
            self.del_component(self._p_sat_dewP)

# -----------------------------------------------------------------------------
# Liquid phase properties
    def _dens_mol_liq(b):
        return b.dens_mol_phase['Liq'] == 1e3*sum(
                b.mole_frac_phase_comp['Liq', j] *
                b._params.dens_liq_params[j, '1'] /
                b._params.dens_liq_params[j, '2'] **
                (1 + (1-b.temperature /
                      b._params.dens_liq_params[j, '3']) **
                 b._params.dens_liq_params[j, '4'])
                for j in ['benzene', 'toluene'])                                # TODO: Need to include diphenyl here later

    def _fug_liq(self):
        def fug_liq_rule(b, i):
            if i in ['hydrogen', 'methane']:
                return b.mole_frac_phase_comp['Liq', i]
            else:
                return b.pressure_sat[i] * b.mole_frac_phase_comp['Liq', i]
        self.fug_liq = Expression(self._params.component_list,
                                  rule=fug_liq_rule)

    def _pressure_sat(self):
        self.pressure_sat = Var(self._params.component_list,
                                initialize=101325,
                                doc="Vapor pressure [Pa]")

        def rule_P_sat(b, j):
            return ((log10(b.pressure_sat[j]*1e-5) -
                     b._params.pressure_sat_coeff[j, 'A']) *
                    (b._teq + b._params.pressure_sat_coeff[j, 'C'])) == \
                   -b._params.pressure_sat_coeff[j, 'B']
        self.eq_pressure_sat = Constraint(self._params.component_list,
                                          rule=rule_P_sat)

    def _enth_mol_comp_liq(b, j):
        return b.enth_mol_phase_comp['Liq', j] * 1E3 == \
                ((b._params.cp_ig['Liq', j, '5'] / 5) *
                    (b.temperature**5 - b._params.temperature_ref**5)
                    + (b._params.cp_ig['Liq', j, '4'] / 4) *
                      (b.temperature**4 - b._params.temperature_ref**4)
                    + (b._params.cp_ig['Liq', j, '3'] / 3) *
                      (b.temperature**3 - b._params.temperature_ref**3)
                    + (b._params.cp_ig['Liq', j, '2'] / 2) *
                      (b.temperature**2 - b._params.temperature_ref**2)
                    + b._params.cp_ig['Liq', j, '1'] *
                      (b.temperature - b._params.temperature_ref))

    def _entr_mol_comp_liq(b, j):
        return b.entr_mol_phase_comp['Liq', j] * 1E3 == (
                ((b._params.cp_ig['Liq', j, '5'] / 4) *
                    (b.temperature**4 - b._params.temperature_ref**4)
                    + (b._params.cp_ig['Liq', j, '4'] / 3) *
                      (b.temperature**3 - b._params.temperature_ref**3)
                    + (b._params.cp_ig['Liq', j, '3'] / 2) *
                      (b.temperature**2 - b._params.temperature_ref**2)
                    + b._params.cp_ig['Liq', j, '2'] *
                      (b.temperature - b._params.temperature_ref)
                    + b._params.cp_ig['Liq', j, '1'] *
                      log(b.temperature / b._params.temperature_ref)) -
                const.gas_constant * log(
                        b.mole_frac_phase_comp['Liq', j]*b.pressure /
                        b._params.pressure_ref))

# -----------------------------------------------------------------------------
# Vapour phase properties
    def _dens_mol_vap(b):
        return b.pressure == (b.dens_mol_phase['Vap'] *
                              const.gas_constant *
                              b.temperature)

    def _fug_vap(self):
        def fug_vap_rule(b, i):
            if i in ['hydrogen', 'methane']:
                return 1e-6
            else:
                return b.mole_frac_phase_comp['Vap', i] * b.pressure
        self.fug_vap = Expression(self._params.component_list,
                                  rule=fug_vap_rule)

    def _dh_vap(self):
        # heat of vaporization
        add_object_reference(self, "dh_vap",
                             self._params.dh_vap)

    def _ds_vap(self):
        # entropy of vaporization = dh_Vap/T_boil
        # TODO : something more rigorous would be nice
        self.ds_vap = Var(self._params.component_list,
                          initialize=86,
                          doc="Entropy of vaporization [J/mol.K]")

        def rule_ds_vap(b, j):
            return b.dh_vap[j] == (b.ds_vap[j] *
                                   b._params.temperature_boil[j])
        self.eq_ds_vap = Constraint(self._params.component_list,
                                    rule=rule_ds_vap)

    def _enth_mol_comp_vap(b, j):
        return b.enth_mol_phase_comp['Vap', j] == b.dh_vap[j] + \
                ((b._params.cp_ig['Vap', j, '5'] / 5) *
                    (b.temperature**5 - b._params.temperature_ref**5)
                    + (b._params.cp_ig['Vap', j, '4'] / 4) *
                      (b.temperature**4 - b._params.temperature_ref**4)
                    + (b._params.cp_ig['Vap', j, '3'] / 3) *
                      (b.temperature**3 - b._params.temperature_ref**3)
                    + (b._params.cp_ig['Vap', j, '2'] / 2) *
                      (b.temperature**2 - b._params.temperature_ref**2)
                    + b._params.cp_ig['Vap', j, '1'] *
                      (b.temperature - b._params.temperature_ref))

    def _entr_mol_comp_vap(b, j):
        return b.entr_mol_phase_comp['Vap', j] == (
                b.ds_vap[j] +
                ((b._params.cp_ig['Vap', j, '5'] / 4) *
                    (b.temperature**4 - b._params.temperature_ref**4)
                    + (b._params.cp_ig['Vap', j, '4'] / 3) *
                      (b.temperature**3 - b._params.temperature_ref**3)
                    + (b._params.cp_ig['Vap', j, '3'] / 2) *
                      (b.temperature**2 - b._params.temperature_ref**2)
                    + b._params.cp_ig['Vap', j, '2'] *
                      (b.temperature - b._params.temperature_ref)
                    + b._params.cp_ig['Vap', j, '1'] *
                      log(b.temperature / b._params.temperature_ref)) -
                const.gas_constant * log(
                        b.mole_frac_phase_comp['Vap', j]*b.pressure /
                        b._params.pressure_ref))
