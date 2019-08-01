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
Ideal gas + Ideal/Non-ideal liquid property package.

VLE calucations assuming an ideal gas for the gas phase. For the liquid phase,
options include ideal liquid or non-ideal liquid using an activity
coefficient model; options include Non Random Two Liquid Model (NRTL) or the
Wilson model to compute the activity coefficient. This property package
supports the following combinations for gas-liquid mixtures:
1. Ideal (vapor) - Ideal (liquid)
2. Ideal (vapor) - NRTL (liquid)
3. Ideal (vapor) - Wilson (liquid)

This property package currently supports the flow_mol, temperature, pressure
and mole_frac as state variables (mole basis). Support for other combinations
will be available in the future.

Please note that the parameters required to compute the activity coefficient
for the component needs to be provided by the user in the parameter block or
can be estimated by the user if VLE data is available. Please see the
documentation for more details.

SI units.

References:

1. "The properties of gases and liquids by Robert C. Reid"
2. "Perry's Chemical Engineers Handbook by Robert H. Perry".
3. H. Renon and J.M. Prausnitz, "Local compositions in thermodynamic excess
   functions for liquid mixtures.", AIChE Journal Vol. 14, No.1, 1968.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, log, NonNegativeReals, value, Var, exp,\
    Set, Expression, Param, sqrt
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util.exceptions import ConfigurationError
from idaes.ui.report import degrees_of_freedom

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.2"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("ActivityCoeffParameterBlock")
class ActivityCoeffParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    Contains parameters and indexing sets associated with properties for
    BTX system.
    """
    # Config block for the _IdealStateBlock
    CONFIG = PhysicalParameterBlock.CONFIG()

    CONFIG.declare("activity_coeff_model", ConfigValue(
        default="Ideal",
        domain=In(['Ideal', 'NRTL', 'Wilson']),
        description="Flag indicating the activity coefficient model",
        doc="""Flag indicating the activity coefficient model to be used
for the non-ideal liquid, and thus corresponding constraints  should be
included,
**default** - Ideal liquid.
**Valid values:** {
**'NRTL'** - Non Random Two Liquid Model,
**'Wilson'** - Wilson Liquid Model,}"""))

    CONFIG.declare("valid_phase", ConfigValue(
        default=('Vap', 'Liq'),
        domain=In(['Liq', 'Vap', ('Vap', 'Liq'), ('Liq', 'Vap')]),
        description="Flag indicating the valid phase",
        doc="""Flag indicating the valid phase for a given set of
conditions, and thus corresponding constraints  should be included,
**default** - ('Vap', 'Liq').
**Valid values:** {
**'Liq'** - Liquid only,
**'Vap'** - Vapor only,
**('Vap', 'Liq')** - Vapor-liquid equilibrium,
**('Liq', 'Vap')** - Vapor-liquid equilibrium,}"""))

    def build(self):
        """Callable method for Block construction."""
        super(ActivityCoeffParameterData, self).build()

        self.state_block_class = ActivityCoeffStateBlock

        # List of valid phases in property package
        if self.config.valid_phase == ('Liq', 'Vap') or \
                self.config.valid_phase == ('Vap', 'Liq'):
            self.phase_list = Set(initialize=['Liq', 'Vap'],
                                  ordered=True)
        elif self.config.valid_phase == 'Liq':
            self.phase_list = Set(initialize=['Liq'])
        else:
            self.phase_list = Set(initialize=['Vap'])

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        obj.add_properties(
            {'flow_mol': {'method': None, 'units': 'mol/s'},
             'mole_frac': {'method': None, 'units': 'no unit'},
             'temperature': {'method': None, 'units': 'K'},
             'pressure': {'method': None, 'units': 'Pa'},
             'flow_mol_phase': {'method': None, 'units': 'mol/s'},
             'density_mol': {'method': '_density_mol',
                             'units': 'mol/m^3'},
             'pressure_sat': {'method': '_pressure_sat', 'units': 'Pa'},
             'mole_frac_phase': {'method': '_mole_frac_phase',
                                 'units': 'no unit'},
             'enthalpy_comp_liq': {'method': '_enthalpy_comp_liq',
                                   'units': 'J/mol'},
             'enthalpy_comp_vap': {'method': '_enthalpy_comp_vap',
                                   'units': 'J/mol'},
             'enthalpy_liq': {'method': '_enthalpy_liq',
                              'units': 'J/mol'},
             'enthalpy_vap': {'method': '_enthalpy_vap',
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
             'fug_liq': {'method': '_fug_liq', 'units': 'Pa'}})

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _ActivityCoeffStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def initialize(blk, flow_mol=None, mole_frac=None,
                   temperature=None, pressure=None,
                   hold_state=False, state_vars_fixed=False, outlvl=1,
                   solver='ipopt', optarg={'tol': 1e-8}):
        """
        Initialisation routine for property package.
        Keyword Arguments:
            flow_mol_comp : value at which to initialize component flows
                             (default=None)
            pressure : value at which to initialize pressure (default=None)
            temperature : value at which to initialize temperature
                          (default=None)
            outlvl : sets output level of initialisation routine
                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = include solver output infomation (tee=True)
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
        """
        # Deactivate the constraints specific for outlet block i.e.
        # when defined state is False
        for k in blk.keys():
            if (blk[k].config.defined_state is False):
                blk[k].eq_mol_frac_out.deactivate()

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            Fflag = {}
            Xflag = {}
            Pflag = {}
            Tflag = {}

            for k in blk.keys():
                if blk[k].flow_mol.fixed is True:
                    Fflag[k] = True
                else:
                    Fflag[k] = False
                    if flow_mol is None:
                        blk[k].flow_mol.fix(1.0)
                    else:
                        blk[k].flow_mol.fix(flow_mol)

                for j in blk[k]._params.component_list:
                    if blk[k].mole_frac[j].fixed is True:
                        Xflag[k, j] = True
                    else:
                        Xflag[k, j] = False
                        if mole_frac is None:
                            blk[k].mole_frac[j].fix(1 / len(blk[k].
                                                    _params.component_list))
                        else:
                            blk[k].mole_frac[j].fix(mole_frac[j])

                if blk[k].pressure.fixed is True:
                    Pflag[k] = True
                else:
                    Pflag[k] = False
                    if pressure is None:
                        blk[k].pressure.fix(101325.0)
                    else:
                        blk[k].pressure.fix(pressure)

                if blk[k].temperature.fixed is True:
                    Tflag[k] = True
                else:
                    Tflag[k] = False
                    if temperature is None:
                        blk[k].temperature.fix(300)
                    else:
                        blk[k].temperature.fix(temperature)

            # ---------------------------------------------------------------------
            # If input block, return flags, else release state
            flags = {"Fflag": Fflag,
                     "Xflag": Xflag,
                     "Pflag": Pflag,
                     "Tflag": Tflag}
        else:
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
        # Initialization sequence: Deactivating certain constraints
        # for 1st solve
        for k in blk.keys():

            blk[k].eq_total.deactivate()
            blk[k].eq_comp.deactivate()
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Vap', 'Liq')):
                blk[k].eq_sum_mol_frac.deactivate()
                blk[k].eq_phase_equilibrium.deactivate()
                try:
                    blk[k].eq_h_liq.deactivate()
                except AttributeError:
                    pass
                try:
                    blk[k].eq_h_vap.deactivate()
                except AttributeError:
                    pass

                # Deactivate activity coefficient constraints
                if blk[k].config.parameters.config.activity_coeff_model \
                        != "Ideal":
                    blk[k].eq_Gij_coeff.deactivate()
                    blk[k].eq_A.deactivate()
                    blk[k].eq_B.deactivate()
                    blk[k].eq_activity_coeff.deactivate()

            # Deactivate liquid phase specific constraints
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k].config.parameters.config.valid_phase == "Liq":
                try:
                    blk[k].eq_h_liq.deactivate()
                except AttributeError:
                    pass

            # Deactivate liquid phase specific constraints
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k].config.parameters.config.valid_phase == "Vap":
                try:
                    blk[k].eq_h_vap.deactivate()
                except AttributeError:
                    pass
        # First solve for the active constraints remaining
        if (blk[k].config.has_phase_equilibrium) or \
                (blk[k].config.parameters.config.valid_phase ==
                    ('Liq', 'Vap')) or \
                (blk[k].config.parameters.config.valid_phase ==
                    ('Vap', 'Liq')):
            results = solve_indexed_blocks(opt, [blk], tee=stee)

            if outlvl > 0:
                if results.solver.termination_condition \
                        == TerminationCondition.optimal:
                    _log.info("Initialisation step 1 for "
                              "{} completed".format(blk.name))
                else:
                    _log.warning("Initialisation step 1 for "
                                 "{} failed".format(blk.name))

        else:
            if outlvl > 0:
                _log.info("Initialisation step 1 for "
                          "{} skipped".format(blk.name))

        # Continue initialization sequence and activate select constraints
        for k in blk.keys():
            blk[k].eq_total.activate()
            blk[k].eq_comp.activate()
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Vap', 'Liq')):
                blk[k].eq_phase_equilibrium.activate()
                if blk[k].config.parameters.config.activity_coeff_model \
                        != "Ideal":
                    # assume ideal and solve
                    blk[k].activity_coeff_comp.fix(1)
                blk[k].eq_sum_mol_frac.activate()

        # Second solve for the active constraints
        results = solve_indexed_blocks(opt, [blk], tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Initialisation step 2 for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Initialisation step 2 for "
                             "{} failed".format(blk.name))

        # Activate activity coefficient specific constraints
        if blk[k].config.parameters.config.activity_coeff_model \
                != "Ideal":
            for k in blk.keys():
                blk[k].eq_Gij_coeff.activate()
                blk[k].eq_A.activate()
                blk[k].eq_B.activate()

            results = solve_indexed_blocks(opt, [blk], tee=stee)
            if outlvl > 0:
                if results.solver.termination_condition \
                        == TerminationCondition.optimal:
                    _log.info("Initialisation step 3 for "
                              "{} completed".format(blk.name))
                else:
                    _log.warning("Initialisation step 3 for "
                                 "{} failed".format(blk.name))

            for k in blk.keys():
                blk[k].eq_activity_coeff.activate()
                blk[k].activity_coeff_comp.unfix()

            results = solve_indexed_blocks(opt, [blk], tee=stee)
            if outlvl > 0:
                if results.solver.termination_condition \
                        == TerminationCondition.optimal:
                    _log.info("Initialisation step 4 for "
                              "{} completed".format(blk.name))
                else:
                    _log.warning("Initialisation step 4 for "
                                 "{} failed".format(blk.name))

        for k in blk.keys():
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k].config.parameters.config.valid_phase == "Liq":
                try:
                    blk[k].eq_h_liq.activate()
                except AttributeError:
                    pass
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k].config.parameters.config.valid_phase == "Vap":
                try:
                    blk[k].eq_h_vap.activate()
                except AttributeError:
                    pass
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Vap', 'Liq')):
                try:
                    blk[k].eq_h_liq.activate()
                except AttributeError:
                    pass
                try:
                    blk[k].eq_h_vap.activate()
                except AttributeError:
                    pass

        results = solve_indexed_blocks(opt, [blk], tee=stee)
        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Initialisation step 5 for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Initialisation step 5 for "
                             "{} failed".format(blk.name))

        for k in blk.keys():
            if (blk[k].config.defined_state is False):
                blk[k].eq_mol_frac_out.activate()

        if state_vars_fixed is False:
            if hold_state is True:
                return flags
            else:
                blk.release_state(flags)

        if outlvl > 0:
            if outlvl > 0:
                _log.info('Initialisation completed for {}.'.format(blk.name))

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
            if flags['Fflag'][k] is False:
                blk[k].flow_mol.unfix()
            for j in blk[k]._params.component_list:
                if flags['Xflag'][k, j] is False:
                    blk[k].mole_frac[j].unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} State Released.'.format(blk.name))


@declare_process_block_class("ActivityCoeffStateBlock",
                             block_class=_ActivityCoeffStateBlock)
class ActivityCoeffStateBlockData(StateBlockData):
    """An example property package for ideal VLE."""

    def build(self):
        """Callable method for Block construction."""
        super(ActivityCoeffStateBlockData, self).build()

        # Check for valid phase indicator and consistent flags
        if self.config.has_phase_equilibrium and \
                self.config.parameters.config.valid_phase in ['Vap', 'Liq']:
            raise ConfigurationError("Inconsistent inputs. Valid phase"
                                     " flag not set to VL for the state"
                                     " block but has_phase_equilibrium"
                                     " is set to True.")
        self._make_state_vars()
        self._make_vars()
        if not self.config.has_phase_equilibrium and \
                self.config.parameters.config.valid_phase == "Liq":
            self._make_liq_phase_eq()

        if (self.config.has_phase_equilibrium) or \
                (self.config.parameters.config.valid_phase ==
                    ('Liq', 'Vap')) or \
                (self.config.parameters.config.valid_phase ==
                    ('Vap', 'Liq')):
            if self.config.parameters.config.activity_coeff_model == 'NRTL':
                self._make_NRTL_eq()
            if self.config.parameters.config.activity_coeff_model == "Wilson":
                self._make_Wilson_eq()
            self._make_flash_eq()

        if not self.config.has_phase_equilibrium and \
                self.config.parameters.config.valid_phase == "Vap":
            self._make_vap_phase_eq()

    def _make_state_vars(self):
        """List the necessary state variable objects."""
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Component molar flowrate [mol/s]')
        self.mole_frac = Var(self._params.component_list,
                             bounds=(0, 1),
                             initialize=1 / len(self._params.component_list))
        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            doc='State pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               doc='State temperature [K]')

    def _make_vars(self):
        self.flow_mol_phase = Var(self._params.phase_list,
                                  initialize=0.5)

        self.mole_frac_phase = \
            Var(self._params.phase_list,
                self._params.component_list,
                initialize=1 / len(self._params.component_list),
                bounds=(0, 1))

    def _make_liq_phase_eq(self):
        def rule_total_mass_balance(self):
            return self.flow_mol_phase['Liq'] == self.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(self, i):
            return self.flow_mol * self.mole_frac[i] == \
                self.flow_mol_phase['Liq'] * self.mole_frac_phase['Liq', i]
        self.eq_comp = Constraint(self._params.component_list,
                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = \
                Constraint(expr=sum(self.mole_frac[i]
                           for i in self._params.component_list) == 1)

    def _make_vap_phase_eq(self):
        def rule_total_mass_balance(self):
            return self.flow_mol_phase['Vap'] == self.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(self, i):
            return self.flow_mol * self.mole_frac[i] == \
                self.flow_mol_phase['Vap'] * self.mole_frac_phase['Vap', i]
        self.eq_comp = Constraint(self._params.component_list,
                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = \
                Constraint(expr=sum(self.mole_frac[i]
                                    for i in self._params.component_list) == 1)

    def _make_flash_eq(self):

        # Total mole balance
        def rule_total_mass_balance(self):
            return self.flow_mol_phase['Liq'] + \
                self.flow_mol_phase['Vap'] == self.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        # Component mole balance
        def rule_comp_mass_balance(self, i):
            return self.flow_mol * self.mole_frac[i] == \
                self.flow_mol_phase['Liq'] * self.mole_frac_phase['Liq', i] + \
                self.flow_mol_phase['Vap'] * self.mole_frac_phase['Vap', i]
        self.eq_comp = Constraint(self._params.component_list,
                                  rule=rule_comp_mass_balance)

        # sum of mole fractions constraint (sum(x_i)-sum(y_i)=0)
        def rule_mole_frac(self):
            return sum(self.mole_frac_phase['Liq', i]
                       for i in self._params.component_list) -\
                sum(self.mole_frac_phase['Vap', i]
                    for i in self._params.component_list) == 0
        self.eq_sum_mol_frac = Constraint(rule=rule_mole_frac)

        if self.config.defined_state is False:
            # applied at outlet only as complete state information is unknown
            self.eq_mol_frac_out = \
                Constraint(expr=sum(self.mole_frac[i]
                           for i in self._params.component_list) == 1)

        # Smooth Flash Formulation

        # Please refer to Burgard et al., "A Smooth, Square Flash
        # Formulation for Equation Oriented Flowsheet Optimization",
        # Computer Aided Chemical Engineering 44, 871-876, 2018.

        self._temperature_equilibrium = \
            Var(initialize=self.temperature.value,
                doc='Temperature for calculating '
                    'phase equilibrium')

        self._t1 = Var(initialize=self.temperature.value,
                       doc='Intermediate temperature for calculating '
                           'the equilibrium temperature')

        self.eps_1 = Param(default=0.01,
                           mutable=True,
                           doc='Smoothing parameter for equilibrium '
                               'temperature')
        self.eps_2 = Param(default=0.0005,
                           mutable=True,
                           doc='Smoothing parameter for equilibrium '
                               'temperature')

        # Equation #13 in reference cited above
        # Approximation for max(temperature, temperature_bubble)
        def rule_t1(b):
            return b._t1 == 0.5 * \
                (b.temperature + b.temperature_bubble +
                 sqrt((b.temperature - b.temperature_bubble)**2 +
                      b.eps_1**2))
        self._t1_constraint = Constraint(rule=rule_t1)

        # Equation #14 in reference cited above
        # Approximation for min(_t1, temperature_dew)
        # TODO : Add option for supercritical extension
        def rule_teq(b):
            return b._temperature_equilibrium == 0.5 * \
                (b._t1 + b.temperature_dew -
                 sqrt((b._t1 - b.temperature_dew)**2 +
                      b.eps_2**2))
        self._teq_constraint = Constraint(rule=rule_teq)

        def rule_phase_eq(self, i):
            return self.fug_vap[i] == self.fug_liq[i]
        self.eq_phase_equilibrium = Constraint(self._params.component_list,
                                               rule=rule_phase_eq)

    def _make_NRTL_eq(self):

        # NRTL Model specific variables (values to be fixed by user or need to
        # be estimated based on VLE data)
        # See documentation for suggested or typical values.
        self.alpha = Var(self._params.component_list,
                         self._params.component_list,
                         initialize=0.3,
                         doc="Non-randomness parameter for NRTL model")

        self.tau = Var(self._params.component_list,
                       self._params.component_list,
                       initialize=1.0,
                       doc="Binary interaction parameter for NRTL model")

        # NRTL model variables
        self.Gij_coeff = Var(self._params.component_list,
                             self._params.component_list,
                             initialize=1.0,
                             doc="Gij coefficient for use in NRTL model ")

        self.activity_coeff_comp = Var(self._params.component_list,
                                       initialize=1.0,
                                       doc="Activity coefficient of component")

        self.A = Var(self._params.component_list,
                     initialize=1.0,
                     doc="Intermediate variable to compute activity"
                     " coefficient")

        self.B = Var(self._params.component_list,
                     initialize=1.0,
                     doc="Intermediate variable to compute activity"
                     " coefficient")

        def rule_Gij_coeff(self, i, j):
            # i,j component
            if i != j:
                return self.Gij_coeff[i, j] == exp(-self.alpha[i, j] *
                                                   self.tau[i, j])
            else:
                self.Gij_coeff[i, j].fix(1)
                return Constraint.Skip

        self.eq_Gij_coeff = Constraint(self._params.component_list,
                                       self._params.component_list,
                                       rule=rule_Gij_coeff)

        # First sum part in the NRTL equation
        def rule_A(self, i):
            value_1 = sum(self.mole_frac_phase['Liq', j] *
                          self.tau[j, i] * self.Gij_coeff[j, i]
                          for j in self._params.component_list)
            value_2 = sum(self.mole_frac_phase['Liq', k] *
                          self.Gij_coeff[k, i]
                          for k in self._params.component_list)
            return self.A[i] == value_1 / value_2
        self.eq_A = Constraint(self._params.component_list, rule=rule_A)

        # Second sum part in the NRTL equation
        def rule_B(self, i):
            value = sum((self.mole_frac_phase['Liq', j] *
                         self.Gij_coeff[i, j] /
                         sum(self.mole_frac_phase['Liq', k] *
                             self.Gij_coeff[k, j]
                             for k in self._params.component_list)) *
                        (self.tau[i, j] - sum(self.mole_frac_phase['Liq', m] *
                                              self.tau[m, j] *
                         self.Gij_coeff[m, j]
                         for m in self._params.component_list) /
                         sum(self.mole_frac_phase['Liq', k] *
                         self.Gij_coeff[k, j]
                         for k in self._params.component_list))
                        for j in self._params.component_list)
            return self.B[i] == value
        self.eq_B = Constraint(self._params.component_list, rule=rule_B)

        # Activity coefficient using NRTL
        def rule_activity_coeff(self, i):
            return log(self.activity_coeff_comp[i]) == self.A[i] + self.B[i]
        self.eq_activity_coeff = Constraint(self._params.component_list,
                                            rule=rule_activity_coeff)

    def _make_Wilson_eq(self):

        # Wilson Model specific variables (values to be fixed by user or need
        # to be estimated based on VLE data)
        self.vol_mol_comp = Var(self._params.component_list,
                                initialize=1.0,
                                doc="Molar volume of component")

        self.tau = Var(self._params.component_list, self._params.component_list,
                       initialize=1.0,
                       doc="Binary interaction parameter for Wilson model")

        # Wilson model variables
        self.Gij_coeff = Var(self._params.component_list,
                             self._params.component_list,
                             initialize=1.0,
                             doc="Gij coefficient for use in Wilson model ")

        self.activity_coeff_comp = Var(self._params.component_list,
                                       initialize=1.0,
                                       doc="Activity coefficient of component")

        self.A = Var(self._params.component_list,
                     initialize=1.0,
                     doc="Intermediate variable to compute activity"
                     " coefficient")

        self.B = Var(self._params.component_list,
                     initialize=1.0,
                     doc="Intermediate variable to compute activity"
                     " coefficient")

        def rule_Gij_coeff(self, i, j):
            # component i,j
            if i != j:
                return self.Gij_coeff[i, j] == \
                    (self.vol_mol_comp[i] /
                     self.vol_mol_comp[j]) * exp(-self.tau[i, j])
            else:
                self.Gij_coeff[i, j].fix(1)
                return Constraint.Skip

        self.eq_Gij_coeff = Constraint(self._params.component_list,
                                       self._params.component_list,
                                       rule=rule_Gij_coeff)

        # First sum part in Wilson equation
        def rule_A(self, i):
            value_1 = log(sum(self.mole_frac_phase['Liq', j] *
                              self.Gij_coeff[j, i]
                              for j in self._params.component_list))
            return self.A[i] == value_1
        self.eq_A = Constraint(self._params.component_list, rule=rule_A)

        # Second sum part in Wilson equation
        def rule_B(self, i):
            value = sum((self.mole_frac_phase['Liq', j] *
                         self.Gij_coeff[i, j] /
                         sum(self.mole_frac_phase['Liq', k] *
                             self.Gij_coeff[k, j]
                             for k in self._params.component_list))
                        for j in self._params.component_list)
            return self.B[i] == value
        self.eq_B = Constraint(self._params.component_list, rule=rule_B)

        # Activity coefficient using Wilson equation
        def rule_activity_coeff(self, i):
            return log(self.activity_coeff_comp[i]) == \
                1 - self.A[i] - self.B[i]
        self.eq_activity_coeff = Constraint(self._params.component_list,
                                            rule=rule_activity_coeff)

    def _pressure_sat(self):
        self.pressure_sat = Var(self._params.component_list,
                                initialize=101325,
                                doc="vapor pressure ")

        def rule_reduced_temp(self, i):
            # reduced temperature is variable "x" in the documentation
            return (self._params.temperature_critical[i] -
                    self._temperature_equilibrium) / \
                self._params.temperature_critical[i]
        self._reduced_temp = Expression(self._params.component_list,
                                        rule=rule_reduced_temp)

        def rule_P_vap(self, j):
            return (1 - self._reduced_temp[j]) * \
                log(self.pressure_sat[j] /
                    self._params.pressure_critical[j]) == \
                (self._params.pressure_sat_coeff[j, 'A'] *
                 self._reduced_temp[j] +
                 self._params.pressure_sat_coeff[j, 'B'] *
                 self._reduced_temp[j]**1.5 +
                 self._params.pressure_sat_coeff[j, 'C'] *
                 self._reduced_temp[j]**3 +
                 self._params.pressure_sat_coeff[j, 'D'] *
                 self._reduced_temp[j]**6)
        self.eq_P_vap = Constraint(self._params.component_list, rule=rule_P_vap)

    def _fug_vap(self):
        def rule_fug_vap(self, i):
            return self.mole_frac_phase['Vap', i] * self.pressure
        self.fug_vap = Expression(self._params.component_list,
                                  rule=rule_fug_vap)

    def _fug_liq(self):
        def rule_fug_liq(self, i):
            if self.config.parameters.config.\
                    activity_coeff_model == "Ideal":
                return self.mole_frac_phase['Liq', i] * \
                    self.pressure_sat[i]
            else:
                return self.mole_frac_phase['Liq', i] * \
                    self.activity_coeff_comp[i] * self.pressure_sat[i]
        self.fug_liq = Expression(self._params.component_list,
                                  rule=rule_fug_liq)

    def _density_mol(self):
        self.density_mol = Var(self._params.phase_list, doc="Molar density")

        def density_mol_calculation(self, p):
            if p == "Vap":
                return self.pressure == (self.density_mol[p] *
                                         self._params.gas_const *
                                         self.temperature)
            elif p == "Liq":  # TODO: Add a correlation to compute liq density
                _log.warning("Using a place holder for liquid density "
                             "{}. Please provide value or expression to "
                             "compute the liquid density".format(self.name))
                return self.density_mol[p] == 11.1E3  # mol/m3

        try:
            # Try to build constraint
            self.density_mol_calculation = Constraint(
                self._params.phase_list, rule=density_mol_calculation)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.density_mol)
            self.del_component(self.density_mol_calculation)
            raise

    def _enthalpy_comp_liq(self):
        # Liquid phase comp enthalpy (J/mol)
        self.enthalpy_comp_liq = Var(self._params.component_list,
                                     initialize=10000)

        def rule_hl_ig_pc(b, j):
            # 1E3 conversion factor to convert from J/kmol to J/mol
            return self.enthalpy_comp_liq[j] * 1E3 == \
                ((self._params.CpIG['Liq', j, '5'] / 5) *
                    (self.temperature**5 -
                     self._params.temperature_reference**5)
                    + (self._params.CpIG['Liq', j, '4'] / 4) *
                      (self.temperature**4 -
                       self._params.temperature_reference**4)
                    + (self._params.CpIG['Liq', j, '3'] / 3) *
                      (self.temperature**3 -
                       self._params.temperature_reference**3)
                    + (self._params.CpIG['Liq', j, '2'] / 2) *
                      (self.temperature**2 -
                       self._params.temperature_reference**2)
                    + self._params.CpIG['Liq', j, '1'] *
                      (self.temperature - self._params.temperature_reference))
        self.eq_hl_ig_pc = Constraint(self._params.component_list,
                                      rule=rule_hl_ig_pc)

    def _enthalpy_liq(self):
        # Liquid phase enthalpy (J/mol)
        self.enthalpy_liq = Var(initialize=10000)

        def rule_hliq(self):
            return self.enthalpy_liq == \
                sum(self.enthalpy_comp_liq[i] * self.mole_frac_phase['Liq', i]
                    for i in self._params.component_list)
        self.eq_h_liq = Constraint(rule=rule_hliq)

    def _enthalpy_comp_vap(self):
        # Vapor phase component enthalpy (J/mol)
        self.enthalpy_comp_vap = Var(self._params.component_list,
                                     initialize=40000)

        def rule_hv_ig_pc(b, j):
            return self.enthalpy_comp_vap[j] == self._params.dh_vap[j] + \
                ((self._params.CpIG['Vap', j, '5'] / 5) *
                    (self.temperature**5 -
                     self._params.temperature_reference**5)
                    + (self._params.CpIG['Vap', j, '4'] / 4) *
                      (self.temperature**4 -
                       self._params.temperature_reference**4)
                    + (self._params.CpIG['Vap', j, '3'] / 3) *
                      (self.temperature**3 -
                       self._params.temperature_reference**3)
                    + (self._params.CpIG['Vap', j, '2'] / 2) *
                      (self.temperature**2 -
                       self._params.temperature_reference**2)
                    + self._params.CpIG['Vap', j, '1'] *
                      (self.temperature -
                       self._params.temperature_reference))
        self.eq_hv_ig_pc = Constraint(self._params.component_list,
                                      rule=rule_hv_ig_pc)

    def _enthalpy_vap(self):
        # Vapor phase total enthalpy (J/mol)
        self.enthalpy_vap = Var(initialize=10000)

        def rule_hvap(self):
            return self.enthalpy_vap == \
                sum(self.enthalpy_comp_vap[i] * self.mole_frac_phase['Vap', i]
                    for i in self._params.component_list)
        self.eq_h_vap = Constraint(rule=rule_hvap)

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if (p == "Vap") and (j in self._params.component_list):
            return self.flow_mol_phase['Vap'] * self.mole_frac_phase['Vap', j]
        elif (p == "Liq") and (j in self._params.component_list):
            return self.flow_mol_phase['Liq'] * self.mole_frac_phase['Liq', j]
        else:
            return 0

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        if p == "Vap":
            return self.flow_mol_phase['Vap'] * self.enthalpy_vap
        elif p == "Liq":
            return self.flow_mol_phase['Liq'] * self.enthalpy_liq

    def get_material_density_terms(self, p, j):
        """Create material density terms."""
        if p == "Liq":
            if j in self._params.component_list:
                return self.density_mol[p] * self.mole_frac_phase['Liq', j]
            else:
                return 0
        elif p == "Vap":
            if j in self._params.component_list:
                return self.density_mol[p] * self.mole_frac_phase['Vap', j]
            else:
                return 0

    def get_enthalpy_density_terms(self, p):
        """Create enthalpy density terms."""
        if p == "Liq":
            return self.density_mol[p] * self.enthalpy_liq
        elif p == "Vap":
            return self.density_mol[p] * self.enthalpy_vap

    def get_material_flow_basis(self):
        """Declare material flow basis."""
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_mol": self.flow_mol,
                "mole_frac": self.mole_frac,
                "temperature": self.temperature,
                "pressure": self.pressure}

    def model_check(blk):
        """Model checks for property block."""
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
    def _temperature_bubble(self):
        self.temperature_bubble = Var(initialize=298.15,
                                      doc="Bubble point temperature (K)")

        def rule_psat_bubble(m, j):
            return self._params.pressure_critical[j] * \
                exp((self._params.pressure_sat_coeff[j, 'A'] *
                    (1 - self.temperature_bubble /
                    self._params.temperature_critical[j]) +
                    self._params.pressure_sat_coeff[j, 'B'] *
                    (1 - self.temperature_bubble /
                    self._params.temperature_critical[j])**1.5 +
                    self._params.pressure_sat_coeff[j, 'C'] *
                    (1 - self.temperature_bubble /
                    self._params.temperature_critical[j])**3 +
                    self._params.pressure_sat_coeff[j, 'D'] *
                    (1 - self.temperature_bubble /
                    self._params.temperature_critical[j])**6) /
                    (1 - (1 - self.temperature_bubble /
                          self._params.temperature_critical[j])))

        try:
            # Try to build expression
            self._p_sat_bubbleT = Expression(self._params.component_list,
                                             rule=rule_psat_bubble)

            def rule_temp_bubble(self):
                if self.config.parameters.config.activity_coeff_model == \
                        "Ideal":

                    return sum(self.mole_frac[i] *
                               self._p_sat_bubbleT[i]
                               for i in self._params.component_list) - \
                        self.pressure == 0
                elif self.config.parameters.config.\
                        activity_coeff_model == 'NRTL':
                    # NRTL model variables
                    def rule_Gij_coeff_bubble(self, i, j):
                        if i != j:
                            return exp(-self.alpha[i, j] * self.tau[i, j])
                        else:
                            return 1

                    self.Gij_coeff_bubble = Expression(
                        self._params.component_list,
                        self._params.component_list,
                        rule=rule_Gij_coeff_bubble)

                    def rule_A_bubble(self, i):
                        value_1 = sum(self.mole_frac[j] *
                                      self.tau[j, i] *
                                      self.Gij_coeff_bubble[j, i]
                                      for j in self._params.component_list)
                        value_2 = sum(self.mole_frac[k] *
                                      self.Gij_coeff_bubble[k, i]
                                      for k in self._params.component_list)
                        return value_1 / value_2
                    self.A_bubble = Expression(self._params.component_list,
                                               rule=rule_A_bubble)

                    def rule_B_bubble(self, i):
                        value = sum((self.mole_frac[j] *
                                     self.Gij_coeff_bubble[i, j] /
                                     sum(self.mole_frac[k] *
                                         self.Gij_coeff_bubble[k, j]
                                         for k in self._params.component_list)) *
                                    (self.tau[i, j] - sum(self.mole_frac[m] *
                                                          self.tau[m, j] *
                                     self.Gij_coeff_bubble[m, j]
                                     for m in self._params.component_list) /
                                     sum(self.mole_frac[k] *
                                     self.Gij_coeff_bubble[k, j]
                                     for k in self._params.component_list))
                                    for j in self._params.component_list)
                        return value
                    self.B_bubble = Expression(self._params.component_list,
                                               rule=rule_B_bubble)

                    def rule_activity_coeff_bubble(self, i):
                        return exp(self.A_bubble[i] + self.B_bubble[i])
                    self.activity_coeff_comp_bubble = \
                        Expression(self._params.component_list,
                                   rule=rule_activity_coeff_bubble)

                    return sum(self.mole_frac[i] *
                               self.activity_coeff_comp_bubble[i] *
                               self._p_sat_bubbleT[i]
                               for i in self._params.component_list) - \
                        self.pressure == 0
                else:
                    def rule_Gij_coeff_bubble(self, i, j):
                        if i != j:
                            return (self.vol_mol_comp[i] /
                                    self.vol_mol_comp[j]) * \
                                exp(-self.tau[i, j])
                        else:
                            return 1

                    self.Gij_coeff_bubble = \
                        Expression(self._params.component_list,
                                   self._params.component_list,
                                   rule=rule_Gij_coeff_bubble)

                    def rule_A_bubble(self, i):
                        value_1 = log(sum(self.mole_frac[j] *
                                          self.Gij_coeff_bubble[j, i]
                                          for j in self._params.component_list))
                        return value_1
                    self.A_bubble = Expression(self._params.component_list,
                                               rule=rule_A_bubble)

                    def rule_B_bubble(self, i):
                        value = sum((self.mole_frac[j] *
                                     self.Gij_coeff_bubble[i, j] /
                                     sum(self.mole_frac[k] *
                                         self.Gij_coeff_bubble[k, j]
                                         for k in self._params.component_list))
                                    for j in self._params.component_list)
                        return value
                    self.B_bubble = Expression(self._params.component_list,
                                               rule=rule_B_bubble)

                    def rule_activity_coeff_bubble(self, i):
                        return exp(1 - self.A_bubble[i] - self.B_bubble[i])
                    self.activity_coeff_comp_bubble = \
                        Expression(self._params.component_list,
                                   rule=rule_activity_coeff_bubble)

                    return sum(self.mole_frac[i] *
                               self.activity_coeff_comp_bubble[i] *
                               self._p_sat_bubbleT[i]
                               for i in self._params.component_list) - \
                        self.pressure == 0
            self.eq_temperature_bubble = Constraint(rule=rule_temp_bubble)

        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_bubble)
            self.del_component(self._p_sat_bubbleT)

    def _temperature_dew(self):

        self.temperature_dew = Var(initialize=298.15,
                                   doc="Dew point temperature (K)")

        def rule_psat_dew(m, j):
            return self._params.pressure_critical[j] * \
                exp((self._params.pressure_sat_coeff[j, 'A'] *
                    (1 - self.temperature_dew /
                    self._params.temperature_critical[j]) +
                    self._params.pressure_sat_coeff[j, 'B'] *
                    (1 - self.temperature_dew /
                    self._params.temperature_critical[j])**1.5 +
                    self._params.pressure_sat_coeff[j, 'C'] *
                    (1 - self.temperature_dew /
                    self._params.temperature_critical[j])**3 +
                    self._params.pressure_sat_coeff[j, 'D'] *
                    (1 - self.temperature_dew /
                    self._params.temperature_critical[j])**6) /
                    (1 - (1 - self.temperature_dew /
                          self._params.temperature_critical[j])))

        try:
            # Try to build expression
            self._p_sat_dewT = Expression(self._params.component_list,
                                          rule=rule_psat_dew)

            def rule_temp_dew(self):
                if self.config.parameters.config.activity_coeff_model == \
                        "Ideal":
                    return self.pressure * \
                        sum(self.mole_frac[i] /
                            self._p_sat_dewT[i]
                            for i in self._params.component_list) - 1 == 0
                elif self.config.parameters.config.\
                        activity_coeff_model == 'NRTL':
                    # NRTL model variables
                    def rule_Gij_coeff_dew(self, i, j):
                        if i != j:
                            return exp(-self.alpha[i, j] * self.tau[i, j])
                        else:
                            return 1

                    self.Gij_coeff_dew = Expression(
                        self._params.component_list,
                        self._params.component_list,
                        rule=rule_Gij_coeff_dew)

                    def rule_A_dew(self, i):
                        value_1 = sum(self.mole_frac[j] *
                                      self.tau[j, i] * self.Gij_coeff_dew[j, i]
                                      for j in self._params.component_list)
                        value_2 = sum(self.mole_frac[k] *
                                      self.Gij_coeff_dew[k, i]
                                      for k in self._params.component_list)
                        return value_1 / value_2
                    self.A_dew = Expression(self._params.component_list,
                                            rule=rule_A_dew)

                    def rule_B_dew(self, i):
                        value = sum((self.mole_frac[j] *
                                     self.Gij_coeff_dew[i, j] /
                                     sum(self.mole_frac[k] *
                                         self.Gij_coeff_dew[k, j]
                                         for k in self._params.component_list)) *
                                    (self.tau[i, j] - sum(self.mole_frac[m] *
                                                          self.tau[m, j] *
                                     self.Gij_coeff_dew[m, j]
                                     for m in self._params.component_list) /
                                     sum(self.mole_frac[k] *
                                     self.Gij_coeff_dew[k, j]
                                     for k in self._params.component_list))
                                    for j in self._params.component_list)
                        return value
                    self.B_dew = Expression(self._params.component_list,
                                            rule=rule_B_dew)

                    def rule_activity_coeff_dew(self, i):
                        return exp(self.A_dew[i] + self.B_dew[i])
                    self.activity_coeff_comp_dew = \
                        Expression(self._params.component_list,
                                   rule=rule_activity_coeff_dew)

                    return sum(self.mole_frac[i] *
                               self.pressure /
                               (self.activity_coeff_comp[i] *
                                self._p_sat_dewT[i])
                               for i in self._params.component_list) - 1 == 0
                else:
                    def rule_Gij_coeff_dew(self, i, j):
                        if i != j:
                            return (self.vol_mol_comp[i] /
                                    self.vol_mol_comp[j]) * \
                            exp(-self.tau[i, j])
                        else:
                            return 1

                    self.Gij_coeff_dew = \
                        Expression(self._params.component_list,
                                   self._params.component_list,
                                   rule=rule_Gij_coeff_dew)

                    def rule_A_dew(self, i):
                        value_1 = log(sum(self.mole_frac[j] *
                                          self.Gij_coeff_dew[j, i]
                                          for j in self._params.component_list))
                        return value_1
                    self.A_dew = Expression(self._params.component_list,
                                            rule=rule_A_dew)

                    def rule_B_dew(self, i):
                        value = sum((self.mole_frac[j] *
                                     self.Gij_coeff_dew[i, j] /
                                     sum(self.mole_frac[k] *
                                         self.Gij_coeff_dew[k, j]
                                         for k in self._params.component_list))
                                    for j in self._params.component_list)
                        return value
                    self.B_dew = Expression(self._params.component_list,
                                            rule=rule_B_dew)

                    def rule_activity_coeff_dew(self, i):
                        return exp(1 - self.A_dew[i] - self.B_dew[i])
                    self.activity_coeff_comp_dew = \
                        Expression(self._params.component_list,
                                   rule=rule_activity_coeff_dew)

                    return sum(self.mole_frac[i] *
                               self.pressure /
                               (self.activity_coeff_comp[i] *
                                self._p_sat_dewT[i])
                               for i in self._params.component_list) - 1 == 0
            self.eq_temperature_dew = Constraint(rule=rule_temp_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_dew)
            self.del_component(self._p_sat_dewT)
