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

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, Expression, log, NonNegativeReals,\
    value, Var, exp, Set, Param, sqrt
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
from idaes.core.util.misc import add_object_reference
from idaes.core.util.exceptions import ConfigurationError
from idaes.ui.report import degrees_of_freedom

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
             'mole_frac': {'method': None, 'units': 'none'},
             'temperature': {'method': None, 'units': 'K'},
             'pressure': {'method': None, 'units': 'Pa'},
             'flow_mol_phase': {'method': None, 'units': 'mol/s'},
             'dens_mol_phase': {'method': '_dens_mol_phase',
                                'units': 'mol/m^3'},
             'pressure_sat': {'method': '_pressure_sat', 'units': 'Pa'},
             'mole_frac_phase': {'method': '_mole_frac_phase',
                                 'units': 'no unit'},
             'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp',
                                     'units': 'J/mol'},
             'enth_mol_phase': {'method': '_enth_mol_phase',
                                'units': 'J/mol'},
             'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp',
                                     'units': 'J/mol'},
             'entr_mol_phase': {'method': '_entr_mol_phase',
                                'units': 'J/mol'},
             'temperature_bub': {'method': '_temperature_bub',
                                 'units': 'K'},
             'temperature_dew': {'method': '_temperature_dew',
                                 'units': 'K'},
             'pressure_bub': {'method': '_pressure_bub',
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

    def initialize(blk, flow_mol=None, mole_frac=None,
                   temperature=None, pressure=None, state_vars_fixed=False,
                   hold_state=False, outlvl=1,
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

        _log.info('Starting {} initialisation'.format(blk.name))

        # Deactivate the constraints specific for outlet block i.e.
        # when defined state is False
        for k in blk.keys():
            if blk[k].config.defined_state is False:
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

                for j in blk[k].component_list_ref:
                    if blk[k].mole_frac[j].fixed is True:
                        Xflag[k, j] = True
                    else:
                        Xflag[k, j] = False
                        if mole_frac is None:
                            blk[k].mole_frac[j].fix(1 / len(blk[k].
                                                    component_list_ref))
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
                        blk[k].temperature.fix(325)
                    else:
                        blk[k].temperature.fix(temperature)

            # ---------------------------------------------------------------------
            # If input block, return flags, else release state
            flags = {"Fflag": Fflag,
                     "Xflag": Xflag,
                     "Pflag": Pflag,
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

        if optarg is None:
            sopt = {'tol': 1e-8}
        else:
            sopt = optarg

        opt = SolverFactory('ipopt')
        opt.options = sopt

        # ---------------------------------------------------------------------
        for k in blk.keys():

            blk[k].eq_total.deactivate()
            blk[k].eq_comp.deactivate()
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k]._params.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k]._params.config.valid_phase ==
                        ('Vap', 'Liq')):
                blk[k].equilibrium_constraint.deactivate()
                blk[k].eq_sum_mol_frac.deactivate()
                try:
                    blk[k].eq_h_liq.deactivate()
                except AttributeError:
                    pass
                try:
                    blk[k].eq_h_vap.deactivate()
                except AttributeError:
                    pass
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k]._params.config.valid_phase == "Liq":
                try:
                    blk[k].eq_h_liq.deactivate()
                except AttributeError:
                    pass
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k]._params.config.valid_phase == "Vap":
                try:
                    blk[k].eq_h_vap.deactivate()
                except AttributeError:
                    pass

        if (blk[k].config.has_phase_equilibrium) or \
                (blk[k]._params.config.valid_phase ==
                    ('Liq', 'Vap')) or \
                (blk[k]._params.config.valid_phase ==
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

        for k in blk.keys():
            blk[k].eq_total.activate()
            blk[k].eq_comp.activate()
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k]._params.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k]._params.config.valid_phase ==
                        ('Vap', 'Liq')):
                blk[k].equilibrium_constraint.activate()
                blk[k].eq_sum_mol_frac.activate()

        results = solve_indexed_blocks(opt, [blk], tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Initialisation step 2 for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Initialisation step 2 for "
                             "{} failed".format(blk.name))

        for k in blk.keys():
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k]._params.config.valid_phase == "Liq":
                try:
                    blk[k].eq_h_liq.activate()
                except AttributeError:
                    pass
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k]._params.config.valid_phase == "Vap":
                try:
                    blk[k].eq_h_vap.activate()
                except AttributeError:
                    pass
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k]._params.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k]._params.config.valid_phase ==
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
                _log.info("Initialisation step 3 for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Initialisation step 3 for "
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
            _log.info("Initialisation completed for {}".format(blk.name))

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
            for j in blk[k].component_list_ref:
                if flags['Xflag'][k, j] is False:
                    blk[k].mole_frac[j].unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()

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

        # List of valid phases in property package
        add_object_reference(self, "phase_list_ref",
                             self._params.phase_list)

        # Component list - a list of component identifiers
        add_object_reference(self, "component_list_ref",
                             self._params.component_list)

        # Add state variables
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Component molar flowrate [mol/s]')
        self.mole_frac = Var(self.component_list_ref,
                             bounds=(0, 1),
                             initialize=1 / len(self.component_list_ref),
                             doc='Mixture mole fractions [-]')
        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            doc='State pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               doc='State temperature [K]')

        # Add supporting variables
        self.flow_mol_phase = Var(self.phase_list_ref,
                                  initialize=0.5,
                                  doc='Phase molar flow rates [mol/s]')

        self.mole_frac_phase = Var(self.phase_list_ref,
                                   self.component_list_ref,
                                   initialize=1 / len(self.component_list_ref),
                                   bounds=(0, 1),
                                   doc='Phase mole fractions [-]')

        if not self.config.has_phase_equilibrium and \
                self._params.config.valid_phase == "Liq":
            self._make_liq_phase_eq()

        if (self.config.has_phase_equilibrium) or \
                (self._params.config.valid_phase ==
                    ('Liq', 'Vap')) or \
                (self._params.config.valid_phase ==
                    ('Vap', 'Liq')):
            self._make_flash_eq()

        if not self.config.has_phase_equilibrium and \
                self._params.config.valid_phase == "Vap":
            self._make_vap_phase_eq()

    def _make_liq_phase_eq(self):
        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Liq'] == b.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac[i] == \
                b.flow_mol_phase['Liq'] * b.mole_frac_phase['Liq', i]
        self.eq_comp = Constraint(self.component_list_ref,
                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = Constraint(expr=sum(self.mole_frac[i]
                                              for i in self.component_list_ref)
                                              == 1)

    def _make_vap_phase_eq(self):
        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Vap'] == b.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac[i] == \
                b.flow_mol_phase['Vap'] * b.mole_frac_phase['Vap', i]
        self.eq_comp = Constraint(self.component_list_ref,
                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = Constraint(expr=sum(self.mole_frac[i]
                                              for i in self.component_list_ref)
                                              == 1)

    def _make_flash_eq(self):
        # List of Reaction Indicies
        add_object_reference(self, "phase_equilibrium_idx_ref",
                             self._params.phase_equilibrium_idx)

        # Reaction Stoichiometry
        add_object_reference(self, "phase_equilibrium_list_ref",
                             self._params.phase_equilibrium_list)

        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Liq'] + \
                b.flow_mol_phase['Vap'] == b.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac[i] == \
                b.flow_mol_phase['Liq'] * b.mole_frac_phase['Liq', i] + \
                b.flow_mol_phase['Vap'] * b.mole_frac_phase['Vap', i]
        self.eq_comp = Constraint(self.component_list_ref,
                                  rule=rule_comp_mass_balance)

        def rule_mole_frac(b):
            return sum(b.mole_frac_phase['Liq', i]
                       for i in b.component_list_ref) -\
                sum(b.mole_frac_phase['Vap', i]
                    for i in b.component_list_ref) == 0
        self.eq_sum_mol_frac = Constraint(rule=rule_mole_frac)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = Constraint(expr=sum(self.mole_frac[i]
                                              for i in self.component_list_ref)
                                              == 1)

        # Definition of equilibrium temperature for smooth VLE
        self._teq = Var(initialize=self.temperature.value,
                        doc='Temperature for calculating phase equilibrium')
        self._t1 = Var(initialize=self.temperature.value,
                       doc='Intermediate temperature for calculating Teq')

        self.eps_1 = Param(default=0.01,
                           doc='Smoothing parameter for Teq')
        self.eps_2 = Param(default=0.0005,
                           doc='Smoothing parameter for Teq')

        # PSE paper Eqn 13
        def rule_t1(b):
            return b._t1 == 0.5*(b.temperature + b.temperature_bub +
                                 sqrt((b.temperature-b.temperature_bub)**2 +
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
                self.component_list_ref,
                rule=rule_tr_eq,
                doc='Component reduced temperatures [-]')

        if self.config.has_phase_equilibrium:
            def rule_equilibrium(b, i):
                return b.fug_vap[i] == b.fug_liq[i]
            self.equilibrium_constraint = Constraint(
                    self.component_list_ref, rule=rule_equilibrium)

# -----------------------------------------------------------------------------
# Property Methods
    def _dens_mol_phase(self):
        self.dens_mol_phase = Var(self.phase_list_ref,
                                  doc="Molar density [mol/m^3]")

        def rule_dens_mol_phase(b, p):
            if p == 'Vap':
                return b._dens_mol_vap()
            else:
                return b._dens_mol_liq()
        self.eq_dens_mol_phase = Constraint(self.phase_list_ref,
                                            rule=rule_dens_mol_phase)

    def _enth_mol_phase_comp(self):
        self.enth_mol_phase_comp = Var(
                self.phase_list_ref,
                self.component_list_ref,
                doc='Phase-component molar specific enthalpies [J/mol]')

        def rule_enth_mol_phase_comp(b, p, j):
            if p == 'Vap':
                return b._enth_mol_comp_vap(j)
            else:
                return b._enth_mol_comp_liq(j)
        self.eq_enth_mol_phase_comp = Constraint(
                self.phase_list_ref,
                self.component_list_ref,
                rule=rule_enth_mol_phase_comp)

    def _enth_mol_phase(self):
        self.enth_mol_phase = Var(
                self.phase_list_ref,
                doc='Phase molar specific enthalpies [J/mol]')

        def rule_enth_mol_phase(b, p):
            return b.enth_mol_phase[p] == sum(b.enth_mol_phase_comp[p, i] *
                                              b.mole_frac_phase[p, i]
                                              for i in b.component_list_ref)
        self.eq_enth_mol_phase = Constraint(self.phase_list_ref,
                                            rule=rule_enth_mol_phase)

    def _entr_mol_phase_comp(self):
        self.entr_mol_phase_comp = Var(
                self.phase_list_ref,
                self.component_list_ref,
                doc='Phase-component molar specific entropies [J/mol.K]')

        def rule_entr_mol_phase_comp(b, p, j):
            if p == 'Vap':
                return b._entr_mol_comp_vap(j)
            else:
                return b._entr_mol_comp_liq(j)
        self.eq_entr_mol_phase_comp = Constraint(
                self.phase_list_ref,
                self.component_list_ref,
                rule=rule_entr_mol_phase_comp)

    def _entr_mol_phase(self):
        self.entr_mol_phase = Var(
                self.phase_list_ref,
                doc='Phase molar specific enthropies [J/mol.K]')

        def rule_entr_mol_phase(b, p):
            return b.entr_mol_phase[p] == sum(b.entr_mol_phase_comp[p, i] *
                                              b.mole_frac_phase[p, i]
                                              for i in b.component_list_ref)
        self.eq_entr_mol_phase = Constraint(self.phase_list_ref,
                                            rule=rule_entr_mol_phase)

# -----------------------------------------------------------------------------
# General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if j in self.component_list_ref:
            return self.flow_mol_phase[p] * self.mole_frac_phase[p, j]
        else:
            return 0

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.flow_mol_phase[p] * self.enth_mol_phase[p]

    def get_material_density_terms(self, p, j):
        """Create material density terms."""
        if j in self.component_list_ref:
            return self.dens_mol_phase[p] * self.mole_frac_phase[p, j]
        else:
            return 0

    def get_enthalpy_density_terms(self, p):
        """Create enthalpy density terms."""
        return self.dens_mol_phase[p] * self.enth_mol_phase[p]

    def get_material_flow_basis(b):
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

        if hasattr(self, "eq_bubble_temp"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.temperature_bub,
                                           self.eq_bubble_temp)

        return self.temperature_bub.value

        if clear_components is True:
            self.del_component(self.eq_bubble_temp)
            self.del_component(self._p_sat_bubbleT)
            self.del_component(self.temperature_bub)

    def calculate_dew_point_temperature(self, clear_components=True):
        """"To compute the dew point temperature of the mixture."""

        if hasattr(self, "eq_dew_temp"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.temperature_dew,
                                           self.eq_dew_temp)

        return self.temperature_dew.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_dew_temp)
            self.del_component(self._p_sat_dewT)
            self.del_component(self.temperature_dew)

    def calculate_bubble_point_pressure(self, clear_components=True):
        """"To compute the bubble point pressure of the mixture."""

        if hasattr(self, "eq_bubble_pressure"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.pressure_bub,
                                           self.eq_bubble_pressure)

        return self.pressure_bub.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_bubble_pressure)
            self.del_component(self._p_sat_bubbleP)
            self.del_component(self.pressure_bub)

    def calculate_dew_point_pressure(self, clear_components=True):
        """"To compute the dew point pressure of the mixture."""

        if hasattr(self, "eq_dew_pressure"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.pressure_dew,
                                           self.eq_dew_pressure)

        return self.pressure_dew.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_dew_pressure)
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
    def _temperature_bub(self):
        self.temperature_bub = Var(initialize=298.15,
                                   doc="Bubble point temperature (K)")

        def rule_psat_bubble(b, j):
            return b._params.pressure_crit[j] * \
                exp((b._params.pressure_sat_coeff[j, 'A'] *
                     (1 - b.temperature_bub /
                      b._params.temperature_crit[j]) +
                     b._params.pressure_sat_coeff[j, 'B'] *
                     (1 - b.temperature_bub /
                      b._params.temperature_crit[j])**1.5 +
                     b._params.pressure_sat_coeff[j, 'C'] *
                     (1 - b.temperature_bub /
                      b._params.temperature_crit[j])**3 +
                     b._params.pressure_sat_coeff[j, 'D'] *
                     (1 - b.temperature_bub /
                      b._params.temperature_crit[j])**6) /
                    (1 - (1 - b.temperature_bub /
                          b._params.temperature_crit[j])))
        try:
            # Try to build expression
            self._p_sat_bubbleT = Expression(self.component_list_ref,
                                             rule=rule_psat_bubble)

            def rule_temp_bubble(b):
                return sum(b._p_sat_bubbleT[i] * b.mole_frac[i]
                           for i in b.component_list_ref) - \
                    b.pressure == 0
            self.eq_bubble_temp = Constraint(rule=rule_temp_bubble)

        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_bub)
            self.del_component(self._p_sat_bubbleT)

    def _temperature_dew(self):

        self.temperature_dew = Var(initialize=298.15,
                                   doc="Dew point temperature (K)")

        def rule_psat_dew(b, j):
            return b._params.pressure_crit[j] * \
                exp((b._params.pressure_sat_coeff[j, 'A'] *
                     (1 - b.temperature_dew /
                      b._params.temperature_crit[j]) +
                     b._params.pressure_sat_coeff[j, 'B'] *
                     (1 - b.temperature_dew /
                      b._params.temperature_crit[j])**1.5 +
                     b._params.pressure_sat_coeff[j, 'C'] *
                     (1 - b.temperature_dew /
                      b._params.temperature_crit[j])**3 +
                     b._params.pressure_sat_coeff[j, 'D'] *
                     (1 - b.temperature_dew /
                      b._params.temperature_crit[j])**6) /
                    (1 - (1 - b.temperature_dew /
                          b._params.temperature_crit[j])))

        try:
            # Try to build expression
            self._p_sat_dewT = Expression(self.component_list_ref,
                                          rule=rule_psat_dew)

            def rule_temp_dew(b):
                return b.pressure * sum(b.mole_frac[i] /
                                        b._p_sat_dewT[i]
                                        for i in b.component_list_ref) \
                    - 1 == 0
            self.eq_dew_temp = Constraint(rule=rule_temp_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_dew)
            self.del_component(self._p_sat_dewT)

    def _pressure_bub(self):
        self.pressure_bub = Var(initialize=298.15,
                                doc="Bubble point pressure (Pa)")

        def rule_psat_bubble(b, j):
            return b._params.pressure_crit[j] * \
                exp((b._params.pressure_sat_coeff[j, 'A'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j]) +
                     b._params.pressure_sat_coeff[j, 'B'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j])**1.5 +
                     b._params.pressure_sat_coeff[j, 'C'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j])**3 +
                     b._params.pressure_sat_coeff[j, 'D'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j])**6) /
                    (1 - (1 - b.temperature /
                          b._params.temperature_crit[j])))

        try:
            # Try to build expression
            self._p_sat_bubbleP = Expression(self.component_list_ref,
                                             rule=rule_psat_bubble)

            def rule_pressure_bubble(b):
                return sum(b._p_sat_bubbleP[i] * b.mole_frac[i]
                           for i in b.component_list_ref) \
                    - b.pressure_bub == 0
            self.eq_bubble_pressure = Constraint(rule=rule_pressure_bubble)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.pressure_bub)
            self.del_component(self._p_sat_bubbleP)

    def _pressure_dew(self):
        self.pressure_dew = Var(initialize=298.15,
                                doc="Dew point pressure (Pa)")

        def rule_psat_dew(b, j):
            return b._params.pressure_crit[j] * \
                exp((b._params.pressure_sat_coeff[j, 'A'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j]) +
                     b._params.pressure_sat_coeff[j, 'B'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j])**1.5 +
                     b._params.pressure_sat_coeff[j, 'C'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j])**3 +
                     b._params.pressure_sat_coeff[j, 'D'] *
                     (1 - b.temperature /
                      b._params.temperature_crit[j])**6) /
                    (1 - (1 - b.temperature /
                          b._params.temperature_crit[j])))

        try:
            # Try to build expression
            self._p_sat_dewP = Expression(self.component_list_ref,
                                          rule=rule_psat_dew)

            def rule_pressure_dew(b):
                return b.pressure_dew * \
                    sum(b.mole_frac[i] / b._p_sat_dewP[i]
                        for i in b.component_list_ref) \
                    - 1 == 0
            self.eq_dew_pressure = Constraint(rule=rule_pressure_dew)
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
                b.mole_frac_phase['Liq', j] *
                b._params.dens_liq_params[j, '1'] /
                b._params.dens_liq_params[j, '2'] **
                (1 + (1-b.temperature /
                      b._params.dens_liq_params[j, '3']) **
                 b._params.dens_liq_params[j, '4'])
                for j in b.component_list_ref)

    def _fug_liq(self):
        def fug_liq_rule(b, i):
            return b.pressure_sat[i] * b.mole_frac_phase['Liq', i]
        self.fug_liq = Expression(self.component_list_ref, rule=fug_liq_rule)

    def _pressure_sat(self):
        self.pressure_sat = Var(self.component_list_ref,
                                initialize=101325,
                                doc="Vapor pressure [Pa]")

        def rule_P_sat(b, j):
            return (b._tr_eq[j]) * \
                log(b.pressure_sat[j] / b._params.pressure_crit[j]) == \
                (b._params.pressure_sat_coeff[j, 'A']*(1-b._tr_eq[j]) +
                 b._params.pressure_sat_coeff[j, 'B']*(1-b._tr_eq[j])**1.5 +
                 b._params.pressure_sat_coeff[j, 'C']*(1-b._tr_eq[j])**3 +
                 b._params.pressure_sat_coeff[j, 'D']*(1-b._tr_eq[j])**6)
        self.eq_P_sat = Constraint(self.component_list_ref, rule=rule_P_sat)

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
                b._params.gas_const *
                log(b.mole_frac_phase['Liq', j]*b.pressure/b._params.pressure_ref))

# -----------------------------------------------------------------------------
# Vapour phase properties
    def _dens_mol_vap(b):
        return b.pressure == (b.dens_mol_phase['Vap'] *
                              b._params.gas_const *
                              b.temperature)

    def _fug_vap(self):
        def fug_vap_rule(b, i):
            return b.mole_frac_phase['Vap', i] * b.pressure
        self.fug_vap = Expression(self.component_list_ref, rule=fug_vap_rule)

    def _dh_vap(self):
        # heat of vaporization
        add_object_reference(self, "dh_vap",
                             self._params.dh_vap)

    def _ds_vap(self):
        # entropy of vaporization = dh_Vap/T_boil
        # TODO : something more rigorous would be nice
        self.ds_vap = Var(self.component_list_ref,
                          initialize=86,
                          doc="Entropy of vaporization [J/mol.K]")

        def rule_ds_vap(b, j):
            return b.dh_vap[j] == (b.ds_vap[j] *
                                   b._params.temperature_boil[j])
        self.eq_ds_vap = Constraint(self.component_list_ref,
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
                b._params.gas_const *
                log(b.mole_frac_phase['Vap', j]*b.pressure/b._params.pressure_ref))
