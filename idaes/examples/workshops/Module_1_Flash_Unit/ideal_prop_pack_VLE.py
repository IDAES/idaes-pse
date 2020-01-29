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
Ideal property package with VLE calucations. Correlations to compute
Cp_comp, h_comp and vapor pressure are obtained from "The properties of gases
and liquids by Robert C. Reid" and "Perry's Chemical Engineers Handbook by
Robert H. Perry". SI units.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, log, NonNegativeReals, \
    value, Var, Set, Param, sqrt, log10
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, In

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
from idaes.core.util.exceptions import BurntToast, ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.constants import Constants as const

# Some more inforation about this module
__author__ = "Jaffer Ghouse"
__version__ = "0.0.2"


# Set up logger
_log = logging.getLogger(__name__)


class IdealParameterData(PhysicalParameterBlock):
    """
    Property Parameter Block Class
    Contains parameters and indexing sets associated with properties for
    BTX system.
    """
    # Config block for the _IdealStateBlock
    CONFIG = PhysicalParameterBlock.CONFIG()

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
        '''
        Callable method for Block construction.
        '''
        super(IdealParameterData, self).build()

        self.state_block_class = IdealStateBlock

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

    def initialize(blk, state_args=None, state_vars_fixed=False,
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

                         flow_mol : value at which to initialize component
                                    flows
                         pressure : value at which to initialize pressure
                         temperature : value at which to initialize temperature
                         mole_frac_comp: value at which to initialize the
                                     component mixture mole fraction
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

        # Deactivate the constraints specific for outlet block i.e.
        # when defined state is False
        for k in blk.keys():
            if blk[k].config.defined_state is False:
                blk[k].sum_mole_frac_out.deactivate()

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
            if hasattr(blk[k], "eq_temperature_bubble"):
                calculate_variable_from_constraint(
                        blk[k].temperature_bubble,
                        blk[k].eq_temperature_bubble)

            if hasattr(blk[k], "eq_temperature_dew"):
                calculate_variable_from_constraint(blk[k].temperature_dew,
                                                   blk[k].eq_temperature_dew)

            if hasattr(blk[k], "eq_pressure_bubble"):
                calculate_variable_from_constraint(blk[k].pressure_bubble,
                                                   blk[k].eq_pressure_bubble)

            if hasattr(blk[k], "eq_pressure_dew"):
                calculate_variable_from_constraint(blk[k].pressure_dew,
                                                   blk[k].eq_pressure_dew)

        if outlvl > 0:
            _log.info("Dew and bubble points initialization for "
                      "{} completed".format(blk.name))

        # ---------------------------------------------------------------------
        # If flash, initialize T1 and Teq
        for k in blk.keys():
            if blk[k].config.has_phase_equilibrium:
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
            if blk[k]._params.config.valid_phase == "Liq":
                blk[k].flow_mol_phase['Liq'].value = \
                    blk[k].flow_mol.value

                for j in blk[k]._params.component_list:
                    blk[k].mole_frac_phase_comp['Liq', j].value = \
                        blk[k].mole_frac_comp[j].value

            elif blk[k]._params.config.valid_phase == "Vap":
                blk[k].flow_mol_phase['Vap'].value = \
                    blk[k].flow_mol.value

                for j in blk[k]._params.component_list:
                    blk[k].mole_frac_phase_comp['Vap', j].value = \
                        blk[k].mole_frac_comp[j].value

            else:
                # Seems to work best with default values for phase flows
                for j in blk[k]._params.component_list:
                    blk[k].mole_frac_phase_comp['Vap', j].value = \
                        blk[k].mole_frac_comp[j].value
                    blk[k].mole_frac_phase_comp['Liq', j].value = \
                        blk[k].mole_frac_comp[j].value

                    calculate_variable_from_constraint(
                            blk[k].pressure_sat[j],
                            blk[k].eq_pressure_sat[j])

        # ---------------------------------------------------------------------
        # Solve phase equilibrium constraints
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Deactivate all property constraints
                if c.local_name not in ("total_flow_balance",
                                        "component_flow_balances",
                                        "equilibrium_constraint",
                                        "sum_mole_frac",
                                        "_t1_constraint",
                                        "_teq_constraint",
                                        "eq_pressure_dew",
                                        "eq_pressure_bubble",
                                        "eq_temperature_dew",
                                        "eq_temperature_bubble",
                                        "eq_pressure_sat"):
                    c.deactivate()

        results = solve_indexed_blocks(opt, [blk], tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Phase state initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Phase state initialization for "
                             "{} failed".format(blk.name))

        # ---------------------------------------------------------------------
        # Initialize other properties
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Activate all constraints except sum_mole_frac_out
                if c.local_name not in ("sum_mole_frac_out"):
                    c.activate()

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Property initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Property initialization for "
                             "{} failed".format(blk.name))

        # ---------------------------------------------------------------------
        # Return state to initial conditions
        for k in blk.keys():
            if (blk[k].config.defined_state is False):
                blk[k].sum_mole_frac_out.activate()

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

        # Check for valid phase indicator and consistent flags
        if self.config.has_phase_equilibrium and \
                self._params.config.valid_phase in ['Vap', 'Liq']:
            raise ConfigurationError("Inconsistent inputs. Valid phase"
                                     " flag not set to VL for the state"
                                     " block but has_phase_equilibrium"
                                     " is set to True.")

        # Add state variables
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Component molar flowrate [mol/s]')
        self.mole_frac_comp = Var(
                self._params.component_list,
                bounds=(0, 1),
                initialize=1 / len(self._params.component_list),
                doc='Mixture mole fractions [-]')
        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            doc='State pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               doc='State temperature [K]')

        # Add supporting variables
        self.flow_mol_phase = Var(self._params.phase_list,
                                  initialize=0.5,
                                  doc='Phase molar flow rates [mol/s]')

        self.mole_frac_phase_comp = Var(
            self._params.phase_list,
            self._params.component_list,
            initialize=1 / len(self._params.component_list),
            bounds=(0, 1),
            doc='Phase mole fractions [-]')

        if not self.config.has_phase_equilibrium and \
                self._params.config.valid_phase == "Liq":
            self._make_liq_phase_eq()
        elif not self.config.has_phase_equilibrium and \
                self._params.config.valid_phase == "Vap":
            self._make_vap_phase_eq()
        elif (self.config.has_phase_equilibrium) or \
                (self._params.config.valid_phase ==
                    ('Liq', 'Vap')) or \
                (self._params.config.valid_phase ==
                    ('Vap', 'Liq')):
            self._make_flash_eq()
        else:
            raise BurntToast("{} found unexpected combination of valid_phases "
                             "and has_phase_equilibrium. Please contact the "
                             "IDAES developers with this bug."
                             .format(self.name))

    def _make_liq_phase_eq(self):
        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Liq'] == b.flow_mol
        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i] == b.mole_frac_phase_comp['Liq', i]
        self.component_flow_balances = Constraint(self._params.component_list,
                                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = Constraint(
                expr=1 == sum(self.mole_frac_comp[i]
                              for i in self._params.component_list))

    def _make_vap_phase_eq(self):
        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Vap'] == b.flow_mol
        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac_comp[i] == b.mole_frac_phase_comp['Vap', i]
        self.component_flow_balances = Constraint(self._params.component_list,
                                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = \
                Constraint(expr=1 == sum(self.mole_frac_comp[i]
                           for i in self._params.component_list))

    def _make_flash_eq(self):

        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Liq'] + \
                b.flow_mol_phase['Vap'] == b.flow_mol
        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac_comp[i] == \
                b.flow_mol_phase['Liq'] * b.mole_frac_phase_comp['Liq', i] + \
                b.flow_mol_phase['Vap'] * b.mole_frac_phase_comp['Vap', i]
        self.component_flow_balances = Constraint(self._params.component_list,
                                                  rule=rule_comp_mass_balance)

        def rule_mole_frac(b):
            return sum(b.mole_frac_phase_comp['Liq', i]
                       for i in b._params.component_list) -\
                sum(b.mole_frac_phase_comp['Vap', i]
                    for i in b._params.component_list) == 0
        self.sum_mole_frac = Constraint(rule=rule_mole_frac)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = \
                Constraint(expr=1 == sum(self.mole_frac_comp[i]
                           for i in self._params.component_list))

        if self.config.has_phase_equilibrium:
            # Definition of equilibrium temperature for smooth VLE
            self._teq = Var(initialize=self.temperature.value,
                            doc='Temperature for calculating '
                                'phase equilibrium')
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
                return b._t1 == 0.5 * \
                    (b.temperature + b.temperature_bubble +
                     sqrt((b.temperature - b.temperature_bubble)**2 +
                          b.eps_1**2))
            self._t1_constraint = Constraint(rule=rule_t1)

            # PSE paper Eqn 14
            # TODO : Add option for supercritical extension
            def rule_teq(b):
                return b._teq == 0.5 * (b._t1 + b.temperature_dew -
                                        sqrt((b._t1 - b.temperature_dew)**2 +
                                             b.eps_2**2))
            self._teq_constraint = Constraint(rule=rule_teq)

            def rule_tr_eq(b, i):
                return b._teq / b._params.temperature_crit[i]
            self._tr_eq = Expression(self._params.component_list,
                                     rule=rule_tr_eq,
                                     doc='Component reduced temperatures [-]')

            def rule_equilibrium(b, i):
                return b.fug_vap[i] == b.fug_liq[i]
            self.equilibrium_constraint = \
                Constraint(self._params.component_list, rule=rule_equilibrium)

# -----------------------------------------------------------------------------
# Property Methods
    def _dens_mol_phase(self):
        self.dens_mol_phase = Var(self._params.phase_list,
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
        self.enth_mol_phase_comp = Var(self._params.phase_list,
                                       self._params.component_list,
                                       doc="Phase-component molar specific "
                                           "enthalpies [J/mol]")

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
            return self.flow_mol_phase[p] * self.mole_frac_phase_comp[p, j]
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

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def default_energy_balance_type(self):
        return EnergyBalanceType.enthalpyTotal

    def define_state_vars(self):
        """Define state vars."""
        return {"flow_mol": self.flow_mol,
                "mole_frac_comp": self.mole_frac_comp,
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
# Ideal-Ideal properties allow for the simplifications below
# Other methods require more complex equations with shadow compositions

# For future work, propose the following:
# Core class writes a set of constraints Phi_L_i == Phi_V_i
# Phi_L_i and Phi_V_i make calls to submethods which add shadow compositions
# as needed
    def _temperature_bubble(self):
        self.temperature_bubble = Var(initialize=298.15,
                                      doc="Bubble point temperature (K)")

        def rule_psat_bubble(b, j):
            return 1e5*10**(b._params.pressure_sat_coeff[j, 'A'] -
                            b._params.pressure_sat_coeff[j, 'B'] /
                            (b.temperature_bubble +
                             b._params.pressure_sat_coeff[j, 'C']))
        try:
            # Try to build expression
            self._p_sat_bubbleT = Expression(self._params.component_list,
                                             rule=rule_psat_bubble)

            def rule_temp_bubble(b):
                return sum(b._p_sat_bubbleT[i] * b.mole_frac_comp[i]
                           for i in b._params.component_list) - \
                    b.pressure == 0
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
                                        for i in b._params.component_list) \
                    - 1 == 0
            self.eq_temperature_dew = Constraint(rule=rule_temp_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_dew)
            self.del_component(self._p_sat_dewT)

    def _pressure_bubble(self):
        self.pressure_bubble = Var(initialize=298.15,
                                   doc="Bubble point pressure (Pa)")

        def rule_psat_bubble(b, j):
            return 1e5*10**(b._params.pressure_sat_coeff[j, 'A'] -
                            b._params.pressure_sat_coeff[j, 'B'] /
                            (b.temperature +
                             b._params.pressure_sat_coeff[j, 'C']))

        try:
            # Try to build expression
            self._p_sat_bubbleP = Expression(self._params.component_list,
                                             rule=rule_psat_bubble)

            def rule_pressure_bubble(b):
                return sum(b._p_sat_bubbleP[i] * b.mole_frac_comp[i]
                           for i in b._params.component_list) \
                    - b.pressure_bubble == 0
            self.eq_pressure_bubble = Constraint(rule=rule_pressure_bubble)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.pressure_bubble)
            self.del_component(self._p_sat_bubbleP)

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
                        for i in b._params.component_list) \
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
        return b.dens_mol_phase['Liq'] == 1e3 * sum(
            b.mole_frac_phase_comp['Liq', j] *
            b._params.dens_liq_params[j, '1'] /
            b._params.dens_liq_params[j, '2'] **
            (1 + (1 - b.temperature /
                  b._params.dens_liq_params[j, '3']) **
             b._params.dens_liq_params[j, '4'])
            for j in b._params.component_list)

    def _fug_liq(self):
        def fug_liq_rule(b, i):
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
                    b.mole_frac_phase_comp['Liq', j] * b.pressure /
                    b._params.pressure_ref))

# -----------------------------------------------------------------------------
# Vapour phase properties
    def _dens_mol_vap(b):
        return b.pressure == (b.dens_mol_phase['Vap'] *
                              const.gas_constant *
                              b.temperature)

    def _fug_vap(self):
        def fug_vap_rule(b, i):
            return b.mole_frac_phase_comp['Vap', i] * b.pressure
        self.fug_vap = Expression(self._params.component_list,
                                  rule=fug_vap_rule)

    def _dh_vap(self):

        # heat of vaporization
        add_object_reference(self, "dh_vap", self._params.dh_vap)

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
            b.ds_vap[j] + ((b._params.cp_ig['Vap', j, '5'] / 4) *
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
                    b.mole_frac_phase_comp['Vap', j] * b.pressure /
                    b._params.pressure_ref))
