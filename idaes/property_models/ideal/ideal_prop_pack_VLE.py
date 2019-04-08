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
    value, Var, exp, Set
from pyomo.opt import SolverFactory, TerminationCondition
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
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
                              'units': 'Pa'}})

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
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Vap', 'Liq')):
                blk[k].eq_Keq.deactivate()
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
                    blk[k].config.parameters.config.valid_phase == "Liq":
                try:
                    blk[k].eq_h_liq.deactivate()
                except AttributeError:
                    pass
            if not blk[k].config.has_phase_equilibrium and \
                    blk[k].config.parameters.config.valid_phase == "Vap":
                try:
                    blk[k].eq_h_vap.deactivate()
                except AttributeError:
                    pass

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

        for k in blk.keys():
            blk[k].eq_total.activate()
            blk[k].eq_comp.activate()
            if (blk[k].config.has_phase_equilibrium) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Liq', 'Vap')) or \
                    (blk[k].config.parameters.config.valid_phase ==
                        ('Vap', 'Liq')):
                blk[k].eq_Keq.activate()
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
                self.config.parameters.config.valid_phase in ['Vap', 'Liq']:
            raise ConfigurationError("Inconsistent inputs. Valid phase"
                                     " flag not set to VL for the state"
                                     " block but has_phase_equilibrium"
                                     " is set to True.")

        self._make_params()
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
            self._make_flash_eq()

        if not self.config.has_phase_equilibrium and \
                self.config.parameters.config.valid_phase == "Vap":
            self._make_vap_phase_eq()

    def _make_params(self):
        """Make references to the necessary parameters."""
        # List of valid phases in property package
        add_object_reference(self, "phase_list_ref",
                             self.config.parameters.phase_list)

        # Component list - a list of component identifiers
        add_object_reference(self, "component_list_ref",
                             self.config.parameters.component_list)

        # Thermodynamic reference state
        add_object_reference(self, "pressure_ref_ref",
                             self.config.parameters.pressure_reference)
        add_object_reference(self, "temperature_ref_ref",
                             self.config.parameters.temperature_reference)

        # Gas Constant
        add_object_reference(self, "gas_const_ref",
                             self.config.parameters.gas_const)

        # Critical Properties
        add_object_reference(self, "pressure_critical_ref",
                             self.config.parameters.pressure_critical)
        add_object_reference(self, "temperature_critical_ref",
                             self.config.parameters.temperature_critical)

        # Molecular weights
        add_object_reference(self, "mw_comp_ref",
                             self.config.parameters.mw_comp)

        # Specific Enthalpy Coefficients
        add_object_reference(self, "CpIG_ref",
                             self.config.parameters.CpIG)

        # Vapor pressure coeeficients
        add_object_reference(self, "pressure_sat_coeff_ref",
                             self.config.parameters.pressure_sat_coeff)

        # heat of vaporization
        add_object_reference(self, "dh_vap_ref",
                             self.config.parameters.dh_vap)

    def _make_state_vars(self):
        """List the necessary state variable objects."""
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Component molar flowrate [mol/s]')
        self.mole_frac = Var(self.component_list_ref,
                             bounds=(0, 1),
                             initialize=1 / len(self.component_list_ref))
        self.pressure = Var(initialize=101325,
                            domain=NonNegativeReals,
                            doc='State pressure [Pa]')
        self.temperature = Var(initialize=298.15,
                               domain=NonNegativeReals,
                               doc='State temperature [K]')

    def _make_vars(self):
        self.flow_mol_phase = Var(self.phase_list_ref,
                                  initialize=0.5)

        self.mole_frac_phase = Var(self.phase_list_ref,
                                   self.component_list_ref,
                                   initialize=1 / len(self.component_list_ref),
                                   bounds=(0, 1))

    def _make_liq_phase_eq(self):
        def rule_total_mass_balance(self):
            return self.flow_mol_phase['Liq'] == self.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(self, i):
            return self.flow_mol * self.mole_frac[i] == \
                self.flow_mol_phase['Liq'] * self.mole_frac_phase['Liq', i]
        self.eq_comp = Constraint(self.component_list_ref,
                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = Constraint(expr=sum(self.mole_frac[i]
                                              for i in self.component_list_ref)
                                              == 1)

    def _make_vap_phase_eq(self):
        def rule_total_mass_balance(self):
            return self.flow_mol_phase['Vap'] == self.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(self, i):
            return self.flow_mol * self.mole_frac[i] == \
                self.flow_mol_phase['Vap'] * self.mole_frac_phase['Vap', i]
        self.eq_comp = Constraint(self.component_list_ref,
                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = Constraint(expr=sum(self.mole_frac[i]
                                              for i in self.component_list_ref)
                                              == 1)

    def _make_flash_eq(self):

        self.pressure_sat = Var(self.component_list_ref,
                                initialize=101325,
                                doc="vapor pressure ")

        def rule_total_mass_balance(self):
            return self.flow_mol_phase['Liq'] + \
                self.flow_mol_phase['Vap'] == self.flow_mol
        self.eq_total = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(self, i):
            return self.flow_mol * self.mole_frac[i] == \
                self.flow_mol_phase['Liq'] * self.mole_frac_phase['Liq', i] + \
                self.flow_mol_phase['Vap'] * self.mole_frac_phase['Vap', i]
        self.eq_comp = Constraint(self.component_list_ref,
                                  rule=rule_comp_mass_balance)

        def rule_mole_frac(self):
            return sum(self.mole_frac_phase['Liq', i]
                       for i in self.component_list_ref) -\
                sum(self.mole_frac_phase['Vap', i]
                    for i in self.component_list_ref) == 0
        self.eq_sum_mol_frac = Constraint(rule=rule_mole_frac)

        if self.config.defined_state is False:
            # applied at outlet only
            self.eq_mol_frac_out = Constraint(expr=sum(self.mole_frac[i]
                                              for i in self.component_list_ref)
                                              == 1)

        if self.config.has_phase_equilibrium:
            def rule_Keq(self, i):
                return self.mole_frac_phase['Vap', i] * self.pressure == \
                    self.pressure_sat[i] * self.mole_frac_phase['Liq', i]
            self.eq_Keq = Constraint(self.component_list_ref, rule=rule_Keq)

        def rule_temp_var_x(self, i):
            return 1 - self.temperature / self.temperature_critical_ref[i]
        self.x = Expression(self.component_list_ref, rule=rule_temp_var_x)

        def rule_P_vap(self, j):
            return (1 - self.x[j]) * \
                log(self.pressure_sat[j] / self.pressure_critical_ref[j]) == \
                (self.pressure_sat_coeff_ref[j, 'A'] * self.x[j] +
                 self.pressure_sat_coeff_ref[j, 'B'] * self.x[j]**1.5 +
                 self.pressure_sat_coeff_ref[j, 'C'] * self.x[j]**3 +
                 self.pressure_sat_coeff_ref[j, 'D'] * self.x[j]**6)
        self.eq_P_vap = Constraint(self.component_list_ref, rule=rule_P_vap)

    def _density_mol(self):
        self.density_mol = Var(self.phase_list_ref, doc="Molar density")

        def density_mol_calculation(self, p):
            if p == "Vap":
                return self.pressure == (self.density_mol[p] *
                                         self.gas_const_ref *
                                         self.temperature)
            elif p == "Liq":  # TODO: Add a correlation to compute liq density
                return self.density_mol[p] == 11.1E3  # mol/m3
        try:
            # Try to build constraint
            self.density_mol_calculation = Constraint(
                self.phase_list_ref, rule=density_mol_calculation)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.density_mol)
            self.del_component(self.density_mol_calculation)
            raise

    def _enthalpy_comp_liq(self):
        # Liquid phase comp enthalpy
        self.enthalpy_comp_liq = Var(self.component_list_ref, initialize=10000)

        def rule_hl_ig_pc(b, j):
            return self.enthalpy_comp_liq[j] * 1E3 == \
                ((self.CpIG_ref['Liq', j, '5'] / 5) *
                    (self.temperature**5 - self.temperature_ref_ref**5)
                    + (self.CpIG_ref['Liq', j, '4'] / 4) *
                      (self.temperature**4 - self.temperature_ref_ref**4)
                    + (self.CpIG_ref['Liq', j, '3'] / 3) *
                      (self.temperature**3 - self.temperature_ref_ref**3)
                    + (self.CpIG_ref['Liq', j, '2'] / 2) *
                      (self.temperature**2 - self.temperature_ref_ref**2)
                    + self.CpIG_ref['Liq', j, '1'] *
                      (self.temperature - self.temperature_ref_ref))
        self.eq_hl_ig_pc = Constraint(self.component_list_ref,
                                      rule=rule_hl_ig_pc)

    def _enthalpy_liq(self):
        # Liquid phase enthalpy
        self.enthalpy_liq = Var()

        def rule_hliq(self):
            return self.enthalpy_liq == sum(self.enthalpy_comp_liq[i] *
                                            self.mole_frac_phase['Liq', i]
                                            for i in self.component_list_ref)
        self.eq_h_liq = Constraint(rule=rule_hliq)

    def _enthalpy_comp_vap(self):
        # Vapor comp enthalpy
        self.enthalpy_comp_vap = Var(self.component_list_ref, initialize=40000)

        def rule_hv_ig_pc(b, j):
            return self.enthalpy_comp_vap[j] == self.dh_vap_ref[j] + \
                ((self.CpIG_ref['Vap', j, '5'] / 5) *
                    (self.temperature**5 - self.temperature_ref_ref**5)
                    + (self.CpIG_ref['Vap', j, '4'] / 4) *
                      (self.temperature**4 - self.temperature_ref_ref**4)
                    + (self.CpIG_ref['Vap', j, '3'] / 3) *
                      (self.temperature**3 - self.temperature_ref_ref**3)
                    + (self.CpIG_ref['Vap', j, '2'] / 2) *
                      (self.temperature**2 - self.temperature_ref_ref**2)
                    + self.CpIG_ref['Vap', j, '1'] *
                      (self.temperature - self.temperature_ref_ref))
        self.eq_hv_ig_pc = Constraint(self.component_list_ref,
                                      rule=rule_hv_ig_pc)

    def _enthalpy_vap(self):
        # Vapor phase enthalpy
        self.enthalpy_vap = Var()

        def rule_hvap(self):
            return self.enthalpy_vap == sum(self.enthalpy_comp_vap[i] *
                                            self.mole_frac_phase['Vap', i]
                                            for i in self.component_list_ref)
        self.eq_h_vap = Constraint(rule=rule_hvap)

    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if (p == "Vap") and (j in self.component_list_ref):
            return self.flow_mol_phase['Vap'] * self.mole_frac_phase['Vap', j]
        elif (p == "Liq") and (j in self.component_list_ref):
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
            if j in self.component_list_ref:
                return self.density_mol[p] * self.mole_frac_phase['Liq', j]
            else:
                return 0
        elif p == "Vap":
            if j in self.component_list_ref:
                return self.density_mol[p] * self.mole_frac_phase['Vap', j]
            else:
                return 0

    def get_enthalpy_density_terms(self, p):
        """Create enthalpy density terms."""
        if p == "Liq":
            return self.density_mol[p] * self.enthalpy_liq
        elif p == "Vap":
            return self.density_mol[p] * self.enthalpy_vap

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

    def _temperature_bubble(self):

        self.temperature_bubble = Var(initialize=298.15,
                                      doc="Bubble point temperature (K)")

        def rule_psat_bubble(m, j):
            return self.pressure_critical_ref[j] * \
                exp((self.pressure_sat_coeff_ref[j, 'A'] *
                    (1 - self.temperature_bubble /
                    self.temperature_critical_ref[j]) +
                    self.pressure_sat_coeff_ref[j, 'B'] *
                    (1 - self.temperature_bubble /
                    self.temperature_critical_ref[j])**1.5 +
                    self.pressure_sat_coeff_ref[j, 'C'] *
                    (1 - self.temperature_bubble /
                    self.temperature_critical_ref[j])**3 +
                    self.pressure_sat_coeff_ref[j, 'D'] *
                    (1 - self.temperature_bubble /
                    self.temperature_critical_ref[j])**6) /
                    (1 - (1 - self.temperature_bubble /
                          self.temperature_critical_ref[j])))
        try:
            # Try to build expression
            self._p_sat_bubbleT = Expression(self.component_list_ref,
                                             rule=rule_psat_bubble)

            def rule_temp_bubble(self):
                return sum(self._p_sat_bubbleT[i] * self.mole_frac[i]
                           for i in self.component_list_ref) - \
                    self.pressure == 0
            self.eq_bubble_temp = Constraint(rule=rule_temp_bubble)

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
            return self.pressure_critical_ref[j] * \
                exp((self.pressure_sat_coeff_ref[j, 'A'] *
                    (1 - self.temperature_dew /
                    self.temperature_critical_ref[j]) +
                    self.pressure_sat_coeff_ref[j, 'B'] *
                    (1 - self.temperature_dew /
                    self.temperature_critical_ref[j])**1.5 +
                    self.pressure_sat_coeff_ref[j, 'C'] *
                    (1 - self.temperature_dew /
                    self.temperature_critical_ref[j])**3 +
                    self.pressure_sat_coeff_ref[j, 'D'] *
                    (1 - self.temperature_dew /
                    self.temperature_critical_ref[j])**6) /
                    (1 - (1 - self.temperature_dew /
                          self.temperature_critical_ref[j])))

        try:
            # Try to build expression
            self._p_sat_dewT = Expression(self.component_list_ref,
                                          rule=rule_psat_dew)

            def rule_temp_dew(self):
                return self.pressure * sum(self.mole_frac[i] /
                                           self._p_sat_dewT[i]
                                           for i in self.component_list_ref) \
                    - 1 == 0
            self.eq_dew_temp = Constraint(rule=rule_temp_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.temperature_dew)
            self.del_component(self._p_sat_dewT)

    def _pressure_bubble(self):
        self.pressure_bubble = Var(initialize=298.15,
                                   doc="Bubble point pressure (Pa)")

        def rule_psat_bubble(m, j):
            return self.pressure_critical_ref[j] * \
                exp((self.pressure_sat_coeff_ref[j, 'A'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j]) +
                    self.pressure_sat_coeff_ref[j, 'B'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j])**1.5 +
                    self.pressure_sat_coeff_ref[j, 'C'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j])**3 +
                    self.pressure_sat_coeff_ref[j, 'D'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j])**6) /
                    (1 - (1 - self.temperature /
                          self.temperature_critical_ref[j])))

        try:
            # Try to build expression
            self._p_sat_bubbleP = Expression(self.component_list_ref,
                                             rule=rule_psat_bubble)

            def rule_pressure_bubble(self):
                return sum(self._p_sat_bubbleP[i] * self.mole_frac[i]
                           for i in self.component_list_ref) \
                    - self.pressure_bubble == 0
            self.eq_bubble_pressure = Constraint(rule=rule_pressure_bubble)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.pressure_bubble)
            self.del_component(self._p_sat_bubbleP)

    def _pressure_dew(self):
        self.pressure_dew = Var(initialize=298.15,
                                doc="Dew point pressure (Pa)")

        def rule_psat_dew(m, j):
            return self.pressure_critical_ref[j] * \
                exp((self.pressure_sat_coeff_ref[j, 'A'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j]) +
                    self.pressure_sat_coeff_ref[j, 'B'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j])**1.5 +
                    self.pressure_sat_coeff_ref[j, 'C'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j])**3 +
                    self.pressure_sat_coeff_ref[j, 'D'] *
                    (1 - self.temperature /
                    self.temperature_critical_ref[j])**6) /
                    (1 - (1 - self.temperature /
                          self.temperature_critical_ref[j])))

        try:
            # Try to build expression
            self._p_sat_dewP = Expression(self.component_list_ref,
                                          rule=rule_psat_dew)

            def rule_pressure_dew(self):
                return self.pressure_dew * \
                    sum(self.mole_frac[i] / self._p_sat_dewP[i]
                        for i in self.component_list_ref) \
                    - 1 == 0
            self.eq_dew_pressure = Constraint(rule=rule_pressure_dew)
        except AttributeError:
            # If expression fails, clean up so that DAE can try again later
            # Deleting only var/expression as expression construction will fail
            # first; if it passes then constraint construction will not fail.
            self.del_component(self.pressure_dew)
            self.del_component(self._p_sat_dewP)

    # Property package utility functions
    def calculate_bubble_point_temperature(self, clear_components=True):
        """"To compute the bubble point temperature of the mixture."""

        if hasattr(self, "eq_bubble_temp"):
            # Do not delete components if the block already has the components
            clear_components = False

        calculate_variable_from_constraint(self.temperature_bubble,
                                           self.eq_bubble_temp)

        return self.temperature_bubble.value

        if clear_components is True:
            self.del_component(self.eq_bubble_temp)
            self.del_component(self._p_sat_bubbleT)
            self.del_component(self.temperature_bubble)

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

        calculate_variable_from_constraint(self.pressure_bubble,
                                           self.eq_bubble_pressure)

        return self.pressure_bubble.value

        # Delete the var/constraint created in this method that are part of the
        # IdealStateBlock if the user desires
        if clear_components is True:
            self.del_component(self.eq_bubble_pressure)
            self.del_component(self._p_sat_bubbleP)
            self.del_component(self.pressure_bubble)

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
