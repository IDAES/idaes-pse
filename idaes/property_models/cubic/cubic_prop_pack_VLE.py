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
General Cubic Equation of State property package with VLE calucations.
Correlations to compute Cp_comp, h_comp and vapor pressure are obtained from
"The properties of gases and liquids by Robert C. Reid" and "Perry's Chemical
Engineers Handbook by Robert H. Perry". SI units.
"""

# Chages the divide behavior to not do integer division
from __future__ import division

# Import Python libraries
import math
import logging
import os
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           exp,
                           Expression,
                           ExternalFunction,
                           log,
                           NonNegativeReals,
                           Set,
                           SolverFactory,
                           sqrt,
                           TerminationCondition,
                           Param,
                           PositiveReals,
                           value,
                           Var)
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util.exceptions import BurntToast
from idaes.core.util.model_statistics import degrees_of_freedom


# Set up logger
_log = logging.getLogger(__name__)
_so = os.path.join(os.path.dirname(__file__), "../cubic_eos/cubic_roots.so")


def cubic_roots_available():
    """Make sure the compiled cubic root functions are available. Yes, in
    Windows the .so extention is still used.
    """
    return os.path.isfile(_so)


class CubicEoS(Enum):
    PR = 0
    SRK = 1


EoS_param = {
        CubicEoS.PR: {'u': 2, 'w': -1, 'omegaA': 0.45724, 'coeff_b': 0.07780},
        CubicEoS.SRK: {'u': 1, 'w': 0, 'omegaA': 0.42748, 'coeff_b': 0.08664}
        }


class CubicParameterData(PhysicalParameterBlock):
    """
    General Property Parameter Block Class
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
        super(CubicParameterData, self).build()

        self.state_block_class = CubicStateBlock

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
             'dens_mass_phase': {'method': '_dens_mass_phase',
                                 'units': 'kg/m^3'},
             'pressure_sat': {'method': '_pressure_sat', 'units': 'Pa'},
             'mole_frac_phase': {'method': '_mole_frac_phase',
                                 'units': 'no unit'},
             'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp',
                                     'units': 'J/mol'},
             'enth_mol_phase': {'method': '_enth_mol_phase',
                                'units': 'J/mol'},
             'enth_mol': {'method': '_enth_mol', 'units': 'J/mol'},
             'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp',
                                     'units': 'J/mol'},
             'entr_mol_phase': {'method': '_entr_mol_phase',
                                'units': 'J/mol'},
             'entr_mol': {'method': '_entr_mol', 'units': 'J/mol.K'},
             'fug_phase': {'method': '_fug_phase', 'units': 'Pa'},
             'fug_coeff_phase': {'method': '_fug_coeff_phase', 'units': '-'},
             'gibbs_mol_phase': {'method': '_gibbs_mol_phase',
                                 'units': 'J/mol'},
             'mw': {'method': '_mw', 'units': 'kg/mol'},
             'mw_phase': {'method': '_mw_phase', 'units': 'kg/mol'},
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


class _CubicStateBlock(StateBlock):
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
            flow_mol : value at which to initialize molar flow (default=None)
            mole_frac : dict of values to use when initializing mole fractions
                        (default = None)
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
                blk[k].sum_mole_frac_out.deactivate()

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
        # If present, initialize bubble and dew point calculations
        # Antoine equation
        def antoine_P(b, j, T):
            return 1e5*10**(b._params.antoine[j, '1'] -
                            b._params.antoine[j, '2'] /
                            (T + b._params.antoine[j, '3']))

        # Bubble temperature initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_tbub"):
                Tbub0 = 0
                for j in blk[k]._params.component_list:
                    Tbub0 += value(
                            blk[k].mole_frac[j] *
                            (blk[k]._params.antoine[j, '2'] /
                             (blk[k]._params.antoine[j, '1'] -
                              math.log10(value(blk[k].pressure*1e-5))) -
                             blk[k]._params.antoine[j, '3']))

                err = 1
                counter = 0

                while err > 1e-2 and counter < 100:
                    f = value(sum(antoine_P(blk[k], j, Tbub0) *
                                  blk[k].mole_frac[j]
                                  for j in blk[k]._params.component_list) -
                              blk[k].pressure)
                    df = value(sum(
                            blk[k].mole_frac[j] *
                            blk[k]._params.antoine[j, '2'] *
                            math.log(10)*antoine_P(blk[k], j, Tbub0) /
                            (Tbub0 + blk[k]._params.antoine[j, '3'])**2
                            for j in blk[k]._params.component_list))

                    Tbub1 = Tbub0 - f/df

                    err = abs(Tbub1 - Tbub0)
                    Tbub0 = Tbub1
                    counter += 1

                blk[k].temperature_bubble.value = Tbub0

                for j in blk[k]._params.component_list:
                    blk[k]._mole_frac_tbub[j].value = value(
                            blk[k].mole_frac[j]*blk[k].pressure /
                            antoine_P(blk[k], j, Tbub0))

        # Dew temperature initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_tdew"):
                Tdew0 = 0
                for j in blk[k]._params.component_list:
                    Tdew0 += value(
                            blk[k].mole_frac[j] *
                            (blk[k]._params.antoine[j, '2'] /
                             (blk[k]._params.antoine[j, '1'] -
                              math.log10(value(blk[k].pressure*1e-5))) -
                             blk[k]._params.antoine[j, '3']))

                err = 1
                counter = 0

                while err > 1e-2 and counter < 100:
                    f = value(blk[k].pressure *
                              sum(blk[k].mole_frac[j] /
                                  antoine_P(blk[k], j, Tdew0)
                                  for j in blk[k]._params.component_list) - 1)
                    df = -value(blk[k].pressure*math.log(10) *
                                sum(blk[k].mole_frac[j] *
                                    blk[k]._params.antoine[j, '2'] /
                                    ((Tdew0 +
                                      blk[k]._params.antoine[j, '3'])**2 *
                                    antoine_P(blk[k], j, Tdew0))
                                    for j in blk[k]._params.component_list))

                    Tdew1 = Tdew0 - f/df

                    err = abs(Tdew1 - Tdew0)
                    Tdew0 = Tdew1
                    counter += 1

                blk[k].temperature_dew.value = Tdew0

                for j in blk[k]._params.component_list:
                    blk[k]._mole_frac_tdew[j].value = value(
                            blk[k].mole_frac[j]*blk[k].pressure /
                            antoine_P(blk[k], j, Tdew0))

#        # Bubble pressure initialization
#        for k in blk.keys():
#            if hasattr(blk[k], "_mole_frac_pbub"):
#                blk[k].pressure_bubble.value = value(
#                        sum(blk[k].mole_frac[j] *
#                            antoine_P(blk[k], j, blk[k].temperature)
#                            for j in blk[k]._params.component_list))
#
#                for j in blk[k]._params.component_list:
#                    blk[k]._mole_frac_pbub[j].value = value(
#                            blk[k].mole_frac[j] *
#                            antoine_P(blk[k], j, blk[k].temperature) /
#                            blk[k].pressure_bubble)
#
#                blk[k].pressure_bubble.display()
#                blk[k]._mole_frac_pbub.display()
#
#        # Dew pressure initialization
#        for k in blk.keys():
#            if hasattr(blk[k], "_mole_frac_pdew"):
#                blk[k].pressure_dew.value = value(
#                        sum(1/(blk[k].mole_frac[j] /
#                               antoine_P(blk[k], j, blk[k].temperature))
#                            for j in blk[k]._params.component_list))
#
#                for j in blk[k]._params.component_list:
#                    blk[k]._mole_frac_pdew[j].value = value(
#                            blk[k].mole_frac[j]*blk[k].pressure_bubble /
#                            antoine_P(blk[k], j, blk[k].temperature))

        if outlvl > 0:
            _log.info("Dew and bubble points initialization for "
                      "{} completed".format(blk.name))

        # ---------------------------------------------------------------------
        # If flash, initialize T1 and Teq
        for k in blk.keys():
            if ((blk[k]._params.config.valid_phase == ('Liq', 'Vap')) or
                    (blk[k]._params.config.valid_phase == ('Vap', 'Liq'))):
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
                    blk[k].mole_frac_phase['Liq', j].value = \
                        blk[k].mole_frac[j].value

            elif blk[k]._params.config.valid_phase == "Vap":
                blk[k].flow_mol_phase['Vap'].value = \
                    blk[k].flow_mol.value

                for j in blk[k]._params.component_list:
                    blk[k].mole_frac_phase['Vap', j].value = \
                        blk[k].mole_frac[j].value

            else:
                # Seems to work best with default values for phase flows
                for j in blk[k]._params.component_list:
                    blk[k].mole_frac_phase['Vap', j].value = \
                        blk[k].mole_frac[j].value
                    blk[k].mole_frac_phase['Liq', j].value = \
                        blk[k].mole_frac[j].value

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
            for j in blk[k]._params.component_list:
                if flags['Xflag'][k, j] is False:
                    blk[k].mole_frac[j].unfix()
            if flags['Pflag'][k] is False:
                blk[k].pressure.unfix()
            if flags['Tflag'][k] is False:
                blk[k].temperature.unfix()

        if outlvl > 0:
            if outlvl > 0:
                _log.info('{} states released.'.format(blk.name))

@declare_process_block_class("CubicStateBlock",
                             block_class=_CubicStateBlock)
class CubicStateBlockData(StateBlockData):
    """An general property package for cubic equations of state with VLE."""

    def build(self):
        """Callable method for Block construction."""
        super(CubicStateBlockData, self).build()

        # Add state variables
        self.flow_mol = Var(initialize=1.0,
                            domain=NonNegativeReals,
                            doc='Component molar flowrate [mol/s]')
        self.mole_frac = Var(self._params.component_list,
                             bounds=(0, None),
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
                                  domain=NonNegativeReals,
                                  doc='Phase molar flow rates [mol/s]')

        self.mole_frac_phase = Var(
            self._params.phase_list,
            self._params.component_list,
            initialize=1 / len(self._params.component_list),
            bounds=(0, None),
            doc='Phase mole fractions [-]')

        if self._params.config.valid_phase == "Liq":
            self._make_liq_phase_eq()
        elif self._params.config.valid_phase == "Vap":
            self._make_vap_phase_eq()
        elif ((self._params.config.valid_phase == ('Liq', 'Vap')) or
                (self._params.config.valid_phase == ('Vap', 'Liq'))):
            self._make_flash_eq()
        else:
            raise BurntToast("{} found unexpected value for valid_phases. "
                             "Please contact the "
                             "IDAES developers with this bug."
                             .format(self.name))

    def _make_liq_phase_eq(self):
        # Add supporting equations for Cubic EoS
        self.common_cubic()

        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Liq'] == b.flow_mol
        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac[i] == b.mole_frac_phase['Liq', i]
        self.component_flow_balances = Constraint(self._params.component_list,
                                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = Constraint(
                expr=1 == sum(self.mole_frac[i]
                              for i in self._params.component_list))

    def _make_vap_phase_eq(self):
        # Add supporting equations for Cubic EoS
        self.common_cubic()

        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Vap'] == b.flow_mol
        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.mole_frac[i] == b.mole_frac_phase['Vap', i]
        self.component_flow_balances = Constraint(self._params.component_list,
                                                  rule=rule_comp_mass_balance)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = \
                Constraint(expr=1 == sum(self.mole_frac[i]
                           for i in self._params.component_list))

    def _make_flash_eq(self):

        def rule_total_mass_balance(b):
            return b.flow_mol_phase['Liq'] + \
                b.flow_mol_phase['Vap'] == b.flow_mol
        self.total_flow_balance = Constraint(rule=rule_total_mass_balance)

        def rule_comp_mass_balance(b, i):
            return b.flow_mol * b.mole_frac[i] == \
                b.flow_mol_phase['Liq'] * b.mole_frac_phase['Liq', i] + \
                b.flow_mol_phase['Vap'] * b.mole_frac_phase['Vap', i]
        self.component_flow_balances = Constraint(self._params.component_list,
                                                  rule=rule_comp_mass_balance)

        def rule_mole_frac(b):
            return sum(b.mole_frac_phase['Liq', i]
                       for i in b._params.component_list) -\
                sum(b.mole_frac_phase['Vap', i]
                    for i in b._params.component_list) == 0
        self.sum_mole_frac = Constraint(rule=rule_mole_frac)

        if self.config.defined_state is False:
            # applied at outlet only
            self.sum_mole_frac_out = \
                Constraint(expr=1 == sum(self.mole_frac[i]
                           for i in self._params.component_list))

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

        # Add supporting equations for Cubic EoS
        self.common_cubic()

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

        def rule_fug_phase_eq(b, p, j):
            if p == 'Vap':
                return b._fug_vap_eq(j)
            else:
                return b._fug_liq_eq(j)
        self._fug_phase_eq = Expression(self._params.phase_list,
                                        self._params.component_list,
                                        rule=rule_fug_phase_eq)

        def rule_equilibrium(b, i):
            return b._fug_phase_eq["Vap", i] == b._fug_phase_eq["Liq", i]
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

    def _dens_mass_phase(self):
        self.dens_mass_phase = Var(self._params.phase_list,
                                   doc="Mass density [kg/m^3]")

        def rule_dens_mass_phase(b, p):
            if p == 'Vap':
                return b._dens_mass_vap()
            else:
                return b._dens_mass_liq()
        self.eq_dens_mass_phase = Constraint(self._params.phase_list,
                                             rule=rule_dens_mass_phase)

#    def _enth_mol_phase_comp(self):
#        self.enth_mol_phase_comp = Var(self._params.phase_list,
#                                       self._params.component_list,
#                                       doc="Phase-component molar specific "
#                                           "enthalpies [J/mol]")
#
#        def rule_enth_mol_phase_comp(b, p, j):
#            if p == 'Vap':
#                return b._enth_mol_comp_vap(j)
#            else:
#                return b._enth_mol_comp_liq(j)
#        self.eq_enth_mol_phase_comp = Constraint(
#            self._params.phase_list,
#            self._params.component_list,
#            rule=rule_enth_mol_phase_comp)
#
    def _enth_mol_phase(self):
        self.enth_mol_phase = Var(
            self._params.phase_list,
            doc='Phase molar specific enthalpies [J/mol]')

        def rule_enth_mol_phase(b, p):
            if p == "Vap":
                return b.enth_mol_phase[p] == b._enth_mol_vap()
            else:
                return b.enth_mol_phase[p] == b._enth_mol_liq()
        self.eq_enth_mol_phase = Constraint(self._params.phase_list,
                                            rule=rule_enth_mol_phase)

#    def _entr_mol_phase_comp(self):
#        self.entr_mol_phase_comp = Var(
#            self._params.phase_list,
#            self._params.component_list,
#            doc='Phase-component molar specific entropies [J/mol.K]')
#
#        def rule_entr_mol_phase_comp(b, p, j):
#            if p == 'Vap':
#                return b._entr_mol_comp_vap(j)
#            else:
#                return b._entr_mol_comp_liq(j)
#        self.eq_entr_mol_phase_comp = Constraint(
#            self._params.phase_list,
#            self._params.component_list,
#            rule=rule_entr_mol_phase_comp)
#
    def _enth_mol(self):
        self.enth_mol = Var(
            doc='Mixture molar specific enthalpies [J/mol]')

        def rule_enth_mol(b):
            return b.enth_mol*b.flow_mol == sum(
                    b.flow_mol_phase[p]*b.enth_mol_phase[p]
                    for p in b._params.phase_list)
        self.eq_enth_mol = Constraint(rule=rule_enth_mol)

    def _entr_mol(self):
        self.entr_mol = Var(
            doc='Mixture molar specific entropies [J/mol.K]')

        def rule_entr_mol(b):
            return b.entr_mol*b.flow_mol == sum(
                    b.flow_mol_phase[p]*b.entr_mol_phase[p]
                    for p in b._params.phase_list)
        self.eq_entr_mol = Constraint(rule=rule_entr_mol)

    def _entr_mol_phase(self):
        self.entr_mol_phase = Var(
            self._params.phase_list,
            doc='Phase molar specific entropies [J/mol.K]')

        def rule_entr_mol_phase(b, p):
            if p == "Vap":
                return b.entr_mol_phase[p] == b._entr_mol_vap()
            else:
                return b.entr_mol_phase[p] == b._entr_mol_liq()
        self.eq_entr_mol_phase = Constraint(self._params.phase_list,
                                            rule=rule_entr_mol_phase)

    def _fug_phase(self):
        def rule_fug_phase(b, p, j):
            if p == 'Vap':
                return b._fug_vap(j)
            else:
                return b._fug_liq(j)
        self.fug_phase = Expression(self._params.phase_list,
                                    self._params.component_list,
                                    rule=rule_fug_phase)

    def _fug_coeff_phase(self):
        def rule_fug_coeff_phase(b, p, j):
            if p == 'Vap':
                return b._fug_coeff_vap(j)
            else:
                return b._fug_coeff_liq(j)
        self.fug_coeff_phase = Expression(self._params.phase_list,
                                          self._params.component_list,
                                          rule=rule_fug_coeff_phase)

    def _gibbs_mol_phase(self):
        self.gibbs_mol_phase = Var(
            self._params.phase_list,
            doc='Phase molar specific Gibbs energy [J/mol]')

        def rule_gibbs_mol_phase(b, p):
            return b.gibbs_mol_phase[p] == (
                    b.enth_mol_phase[p] - b.temperature*b.entr_mol_phase[p])
        self.eq_gibbs_mol_phase = Constraint(self._params.phase_list,
                                             rule=rule_gibbs_mol_phase)

    def _mw(self):
        def rule_mw(b):
            return sum(b.mw_phase[p] for p in b._params.phase_list)
        self.mw = Expression(rule=rule_mw)

    def _mw_phase(self):
        def rule_mw_phase(b, p):
            return sum(b.mole_frac_phase[p, j]*b._params.mw_comp[j]
                       for j in b._params.component_list)
        self.mw_phase = Expression(self._params.phase_list,
                                   rule=rule_mw_phase)

# -----------------------------------------------------------------------------
# General Methods
    def get_material_flow_terms(self, p, j):
        """Create material flow terms for control volume."""
        if j in self._params.component_list:
            return self.flow_mol_phase[p] * self.mole_frac_phase[p, j]
        else:
            return 0

    def get_enthalpy_flow_terms(self, p):
        """Create enthalpy flow terms."""
        return self.flow_mol_phase[p] * self.enth_mol_phase[p]

    def get_material_density_terms(self, p, j):
        """Create material density terms."""
        if j in self._params.component_list:
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

# -----------------------------------------------------------------------------
# Bubble and Dew Points
    def _temperature_bubble(self):
        self.temperature_bubble = Var(
                doc="Bubble point temperature (K)")

        self._mole_frac_tbub = Var(
                self._params.component_list,
                initialize=1/len(self._params.component_list),
                bounds=(0, None),
                doc="Vapor mole fractions at bubble point")

        self._sum_mole_frac_tbub = Constraint(
                expr=1 == sum(self._mole_frac_tbub[j]
                              for j in self._params.component_list))

        def rule_bubble_temp(b, j):
            return b.bubble_temp_liq(j) == b.bubble_temp_vap(j)
        self.eq_temperature_bubble = Constraint(self._params.component_list,
                                                rule=rule_bubble_temp)

    def _temperature_dew(self):
        self.temperature_dew = Var(
                doc="Dew point temperature (K)")

        self._mole_frac_tdew = Var(
                self._params.component_list,
                initialize=1/len(self._params.component_list),
                bounds=(0, None),
                doc="Liquid mole fractions at dew point")

        self._sum_mole_frac_tdew = Constraint(
                expr=1 == sum(self._mole_frac_tdew[j]
                              for j in self._params.component_list))

        def rule_dew_temp(b, j):
            return b.dew_temp_liq(j) == b.dew_temp_vap(j)
        self.eq_temperature_dew = Constraint(self._params.component_list,
                                             rule=rule_dew_temp)

    def _pressure_bubble(self):
        self.pressure_bubble = Var(
                domain=PositiveReals,
                doc="Bubble point pressure (Pa)")

        self._mole_frac_pbub = Var(
                self._params.component_list,
                initialize=1/len(self._params.component_list),
                bounds=(0, None),
                doc="Vapor mole fractions at bubble point")

        self._sum_mole_frac_pbub = Constraint(
                expr=1 == sum(self._mole_frac_pbub[j]
                              for j in self._params.component_list))

        def rule_bubble_pres(b, j):
            return b.bubble_pres_liq(j) == b.bubble_pres_vap(j)
        self.eq_pressure_bubble = Constraint(self._params.component_list,
                                             rule=rule_bubble_pres)

    def _pressure_dew(self):
        self.pressure_dew = Var(
                domain=PositiveReals,
                doc="Dew point pressure (Pa)")

        self._mole_frac_pdew = Var(
                self._params.component_list,
                initialize=1/len(self._params.component_list),
                bounds=(0, None),
                doc="Liquid mole fractions at dew point")

        self._sum_mole_frac_pdew = Constraint(
                expr=1 == sum(self._mole_frac_pdew[j]
                              for j in self._params.component_list))

        def rule_dew_press(b, j):
            return b.dew_press_liq(j) == b.dew_press_vap(j)
        self.eq_pressure_dew = Constraint(self._params.component_list,
                                          rule=rule_dew_press)

# -----------------------------------------------------------------------------
# Liquid phase properties
    def _vol_mol_liq(b):
        return b._vol_mol_cubic("Liq")

    def _dens_mol_liq(b):
        return b._dens_mol_cubic("Liq")

    def _dens_mass_liq(b):
        return b._dens_mass_cubic("Liq")

    def _fug_liq(b, j):
        return b._fug_cubic("Liq", j)

    def _fug_liq_eq(b, j):
        return b._fug_cubic_eq("Liq", j)

    def _fug_coeff_liq(b, j):
        return b._fug_coeff_cubic("Liq", j)

    def _enth_mol_liq(b):
        return b._enth_mol_cubic("Liq")

    def _enth_mol_liq_ig(b):
        return b._enth_mol_ig("Liq")

#    def _enth_mol_comp_liq(b, j):
#        return b.enth_mol_phase_comp["Liq", j] == b.dh_liq[j] + \
#            ((b._params.cp_ig["Liq", j, '5'] / 5) *
#             (b.temperature**5 - b._params.temperature_ref**5)
#             + (b._params.cp_ig["Liq", j, '4'] / 4) *
#               (b.temperature**4 - b._params.temperature_ref**4)
#             + (b._params.cp_ig["Liq", j, '3'] / 3) *
#               (b.temperature**3 - b._params.temperature_ref**3)
#             + (b._params.cp_ig["Liq", j, '2'] / 2) *
#               (b.temperature**2 - b._params.temperature_ref**2)
#             + b._params.cp_ig["Liq", j, '1'] *
#               (b.temperature - b._params.temperature_ref))
#
#    def _entr_mol_comp_liq(b, j):
#        return b.entr_mol_phase_comp["Liq", j] == (
#            b.ds_liq[j] + ((b._params.cp_ig["Liq", j, '5'] / 4) *
#                           (b.temperature**4 - b._params.temperature_ref**4)
#                           + (b._params.cp_ig["Liq", j, '4'] / 3) *
#                           (b.temperature**3 - b._params.temperature_ref**3)
#                           + (b._params.cp_ig["Liq", j, '3'] / 2) *
#                           (b.temperature**2 - b._params.temperature_ref**2)
#                           + b._params.cp_ig["Liq", j, '2'] *
#                           (b.temperature - b._params.temperature_ref)
#                           + b._params.cp_ig["Liq", j, '1'] *
#                           log(b.temperature / b._params.temperature_ref)) -
#            b._params.gas_const * log(b.mole_frac_phase["Liq", j] * b.pressure /
#                                      b._params.pressure_ref))

    def _entr_mol_liq(b):
        return b._entr_mol_cubic("Liq")

    def _entr_mol_liq_ig(b):
        return b._entr_mol_ig("Liq")

    def bubble_temp_liq(b, j):
        def a(k):
            return (b.omegaA*((b._params.gas_const *
                               b._params.temperature_crit[k])**2 /
                              b._params.pressure_crit[k]) *
                    ((1+b.fw[k]*(1-sqrt(b.temperature_bubble /
                                 b._params.temperature_crit[k])))**2))

        am = sum(sum(b.mole_frac[i]*b.mole_frac[j] *
                     sqrt(a(i)*a(j))*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b.mole_frac[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure/(b._params.gas_const*b.temperature_bubble)**2
        B = bm*b.pressure/(b._params.gas_const*b.temperature_bubble)

        delta = (2*sqrt(a(j))/am *
                 sum(b.mole_frac[i]*sqrt(a(i))*(1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_liq(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b.mole_frac[j]*phi

    def dew_temp_liq(b, j):
        def a(k):
            return (b.omegaA*((b._params.gas_const *
                               b._params.temperature_crit[k])**2 /
                              b._params.pressure_crit[k]) *
                    ((1+b.fw[k]*(1-sqrt(b.temperature_dew /
                                 b._params.temperature_crit[k])))**2))

        am = sum(sum(b._mole_frac_tdew[i]*b._mole_frac_tdew[j] *
                     sqrt(a(i)*a(j))*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b._mole_frac_tdew[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure/(b._params.gas_const*b.temperature_dew)**2
        B = bm*b.pressure/(b._params.gas_const*b.temperature_dew)

        delta = (2*sqrt(a(j))/am *
                 sum(b._mole_frac_tdew[i]*sqrt(a(i))*(1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_liq(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b._mole_frac_tdew[j]*phi

    def bubble_pres_liq(b, j):
        am = sum(sum(b.mole_frac[i]*b.mole_frac[j] *
                     sqrt(b.a[i]*b.a[j])*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b.mole_frac[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure_bubble/(b._params.gas_const*b.temperature)**2
        B = bm*b.pressure_bubble/(b._params.gas_const*b.temperature)

        delta = (2*sqrt(b.a[j])/am *
                 sum(b.mole_frac[i]*sqrt(b.a[i])*(1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_liq(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b.mole_frac[j]*phi

    def dew_press_liq(b, j):
        am = sum(sum(b._mole_frac_pdew[i]*b._mole_frac_pdew[j] *
                     sqrt(b.a[i]*b.a[j])*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b._mole_frac_pdew[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure_dew/(b._params.gas_const*b.temperature)**2
        B = bm*b.pressure_dew/(b._params.gas_const*b.temperature)

        delta = (2*sqrt(b.a[j])/am *
                 sum(b._mole_frac_pdew[i]*sqrt(b.a[i]) *
                     (1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_liq(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b._mole_frac_pdew[j]*phi

# -----------------------------------------------------------------------------
# Vapour phase properties
    def _vol_mol_vap(b):
        return b._vol_mol_cubic("Vap")

    def _dens_mol_vap(b):
        return b._dens_mol_cubic("Vap")

    def _dens_mass_vap(b):
        return b._dens_mass_cubic("Vap")

    def _fug_vap(b, j):
        return b._fug_cubic("Vap", j)

    def _fug_vap_eq(b, j):
        return b._fug_cubic_eq("Vap", j)

    def _fug_coeff_vap(b, j):
        return b._fug_coeff_cubic("Vap", j)

    def _enth_mol_vap(b):
        return b._enth_mol_cubic("Vap")

    def _enth_mol_vap_ig(b):
        return b._enth_mol_ig("Vap")

#    def _enth_mol_comp_vap(b, j):
#        return b.enth_mol_phase_comp['Vap', j] == b.dh_vap[j] + \
#            ((b._params.cp_ig['Vap', j, '5'] / 5) *
#             (b.temperature**5 - b._params.temperature_ref**5)
#             + (b._params.cp_ig['Vap', j, '4'] / 4) *
#               (b.temperature**4 - b._params.temperature_ref**4)
#             + (b._params.cp_ig['Vap', j, '3'] / 3) *
#               (b.temperature**3 - b._params.temperature_ref**3)
#             + (b._params.cp_ig['Vap', j, '2'] / 2) *
#               (b.temperature**2 - b._params.temperature_ref**2)
#             + b._params.cp_ig['Vap', j, '1'] *
#               (b.temperature - b._params.temperature_ref))
#
#    def _entr_mol_comp_vap(b, j):
#        return b.entr_mol_phase_comp['Vap', j] == (
#            b.ds_vap[j] + ((b._params.cp_ig['Vap', j, '5'] / 4) *
#                           (b.temperature**4 - b._params.temperature_ref**4)
#                           + (b._params.cp_ig['Vap', j, '4'] / 3) *
#                           (b.temperature**3 - b._params.temperature_ref**3)
#                           + (b._params.cp_ig['Vap', j, '3'] / 2) *
#                           (b.temperature**2 - b._params.temperature_ref**2)
#                           + b._params.cp_ig['Vap', j, '2'] *
#                           (b.temperature - b._params.temperature_ref)
#                           + b._params.cp_ig['Vap', j, '1'] *
#                           log(b.temperature / b._params.temperature_ref)) -
#            b._params.gas_const * log(b.mole_frac_phase['Vap', j] * b.pressure /
#                                      b._params.pressure_ref))

    def _entr_mol_vap(b):
        return b._entr_mol_cubic("Vap")

    def _entr_mol_vap_ig(b):
        return b._entr_mol_ig("Vap")

    def bubble_temp_vap(b, j):
        def a(k):
            return (b.omegaA*((b._params.gas_const *
                               b._params.temperature_crit[k])**2 /
                              b._params.pressure_crit[k]) *
                    ((1+b.fw[k]*(1-sqrt(b.temperature_bubble /
                                 b._params.temperature_crit[k])))**2))

        am = sum(sum(b._mole_frac_tbub[i]*b._mole_frac_tbub[j] *
                     sqrt(a(i)*a(j))*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b._mole_frac_tbub[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure/(b._params.gas_const*b.temperature_bubble)**2
        B = bm*b.pressure/(b._params.gas_const*b.temperature_bubble)

        delta = (2*sqrt(a(j))/am *
                 sum(b._mole_frac_tbub[i]*sqrt(a(i))*(1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_vap(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b._mole_frac_tbub[j]*phi

    def dew_temp_vap(b, j):
        def a(k):
            return (b.omegaA*((b._params.gas_const *
                               b._params.temperature_crit[k])**2 /
                              b._params.pressure_crit[k]) *
                    ((1+b.fw[k]*(1-sqrt(b.temperature_dew /
                                 b._params.temperature_crit[k])))**2))

        am = sum(sum(b.mole_frac[i]*b.mole_frac[j] *
                     sqrt(a(i)*a(j))*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b.mole_frac[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure/(b._params.gas_const*b.temperature_dew)**2
        B = bm*b.pressure/(b._params.gas_const*b.temperature_dew)

        delta = (2*sqrt(a(j))/am *
                 sum(b.mole_frac[i]*sqrt(a(i))*(1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_vap(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b.mole_frac[j]*phi

    def bubble_pres_vap(b, j):
        am = sum(sum(b._mole_frac_pbub[i]*b._mole_frac_pbub[j] *
                     sqrt(b.a[i]*b.a[j])*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b._mole_frac_pbub[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure_bubble/(b._params.gas_const*b.temperature)**2
        B = bm*b.pressure_bubble/(b._params.gas_const*b.temperature)

        delta = (2*sqrt(b.a[j])/am *
                 sum(b._mole_frac_pbub[i]*sqrt(b.a[i]) *
                     (1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_vap(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b._mole_frac_pbub[j]*phi

    def dew_press_vap(b, j):
        am = sum(sum(b.mole_frac[i]*b.mole_frac[j] *
                     sqrt(b.a[i]*b.a[j])*(1-b._params.kappa[i, j])
                     for j in b._params.component_list)
                 for i in b._params.component_list)
        bm = sum(b.mole_frac[i]*b.b[i] for i in b._params.component_list)

        A = am*b.pressure_dew/(b._params.gas_const*b.temperature)**2
        B = bm*b.pressure_dew/(b._params.gas_const*b.temperature)

        delta = (2*sqrt(b.a[j])/am *
                 sum(b.mole_frac[i]*sqrt(b.a[i])*(1-b._params.kappa[j, i])
                     for i in b._params.component_list))

        Z = b.proc_Z_vap(b._ext_func_param, A, B)

        phi = exp((b.b[j]/bm*(Z-1)*(B*b.EoS_p) -
                   log(Z-B)*(B*b.EoS_p) +
                   A*(b.b[j]/bm - delta) *
                   log((2*Z + B*(b.EoS_u + b.EoS_p)) /
                       (2*Z + B*(b.EoS_u-b.EoS_p))))/(B*b.EoS_p))

        return b.mole_frac[j]*phi

# -----------------------------------------------------------------------------
# Common Cubic Functions
    def common_cubic(blk):
        if hasattr(blk, "omegaA"):
            return

        blk.omegaA = EoS_param[blk._params.cubic_type]['omegaA']

        blk.EoS_Bc = EoS_param[blk._params.cubic_type]['coeff_b']
        blk.EoS_u = EoS_param[blk._params.cubic_type]['u']
        blk.EoS_w = EoS_param[blk._params.cubic_type]['w']
        blk.EoS_p = sqrt(blk.EoS_u**2 - 4*blk.EoS_w)

        # Create expressions for coefficients
        def func_fw(b, j):
            if b._params.cubic_type == CubicEoS.PR:
                return 0.37464 + 1.54226*b._params.omega[j] - \
                       0.26992*b._params.omega[j]**2
            elif b._params.cubic_type == CubicEoS.SRK:
                return 0.48 + 1.574*b._params.omega[j] - \
                       0.176*b._params.omega[j]**2
            else:
                raise BurntToast(
                        "{} received unrecognised cubic type. This should "
                        "never happen, so please contact the IDAES developers "
                        "with this bug.".format(b.name))
        blk.fw = Param(blk._params.component_list,
                       initialize=func_fw,
                       doc='EoS S factor')

        def func_b(b, j):
            return b.EoS_Bc*b._params.gas_const *\
                   b._params.temperature_crit[j]/b._params.pressure_crit[j]
        blk.b = Param(blk._params.component_list,
                      initialize=func_b,
                      doc='Component b coefficient')

        def func_a(b, j):
            return (b.omegaA*((b._params.gas_const *
                               b._params.temperature_crit[j])**2 /
                              b._params.pressure_crit[j]) *
                    ((1+b.fw[j]*(1-sqrt(b.temperature /
                                        b._params.temperature_crit[j])))**2))
        blk.a = Expression(blk._params.component_list,
                           rule=func_a,
                           doc='Component a coefficient')

        def func_a_eq(b, j):
            return (b.omegaA*((b._params.gas_const *
                               b._params.temperature_crit[j])**2 /
                              b._params.pressure_crit[j]) *
                    ((1+b.fw[j]*(1-sqrt(b._teq /
                                        b._params.temperature_crit[j])))**2))
        blk._a_eq = Expression(blk._params.component_list,
                               rule=func_a_eq,
                               doc='Component a coefficient at Teq')

        def rule_am(b, p):
            return sum(sum(
                b.mole_frac_phase[p, i]*b.mole_frac_phase[p, j] *
                sqrt(b.a[i]*b.a[j])*(1-b._params.kappa[i, j])
                for j in b._params.component_list)
                for i in b._params.component_list)
        blk.am = Expression(blk._params.phase_list, rule=rule_am)

        def rule_am_eq(b, p):
            return sum(sum(
                b.mole_frac_phase[p, i]*b.mole_frac_phase[p, j] *
                sqrt(b._a_eq[i]*b._a_eq[j])*(1-b._params.kappa[i, j])
                for j in b._params.component_list)
                for i in b._params.component_list)
        blk._am_eq = Expression(blk._params.phase_list, rule=rule_am_eq)

        def rule_bm(b, p):
            return sum(b.mole_frac_phase[p, i]*b.b[i]
                       for i in b._params.component_list)
        blk.bm = Expression(blk._params.phase_list, rule=rule_bm)

        def rule_A(b, p):
            return (b.am[p]*b.pressure /
                    (b._params.gas_const*b.temperature)**2)
        blk.A = Expression(blk._params.phase_list, rule=rule_A)

        def rule_B(b, p):
            return (b.bm[p]*b.pressure /
                    (b._params.gas_const*b.temperature))
        blk.B = Expression(blk._params.phase_list, rule=rule_B)

        def rule_A_eq(b, p):
            return (b._am_eq[p]*b.pressure /
                    (b._params.gas_const*b._teq)**2)
        blk._A_eq = Expression(blk._params.phase_list, rule=rule_A_eq)

        def rule_B_eq(b, p):
            return (b.bm[p]*b.pressure /
                    (b._params.gas_const*b._teq))
        blk._B_eq = Expression(blk._params.phase_list, rule=rule_B_eq)

        blk.proc_Z_liq = ExternalFunction(library=_so,
                                          function="ceos_z_liq")
        blk.proc_Z_vap = ExternalFunction(library=_so,
                                          function="ceos_z_vap")
        blk.proc_Z_liq_x = ExternalFunction(library=_so,
                                            function="ceos_z_liq_extend")
        blk.proc_Z_vap_x = ExternalFunction(library=_so,
                                            function="ceos_z_vap_extend")

        def rule_delta(b, p, i):
            return (2*sqrt(blk.a[i])/b.am[p] *
                    sum(b.mole_frac_phase[p, j]*sqrt(blk.a[j]) *
                        (1-b._params.kappa[i, j])
                        for j in b._params.component_list))
        blk.delta = Expression(blk._params.phase_list,
                               blk._params.component_list,
                               rule=rule_delta)

        def rule_delta_eq(b, p, i):
            return (2*sqrt(blk._a_eq[i])/b._am_eq[p] *
                    sum(b.mole_frac_phase[p, j]*sqrt(blk._a_eq[j]) *
                        (1-b._params.kappa[i, j])
                        for j in b._params.component_list))
        blk._delta_eq = Expression(blk._params.phase_list,
                                   blk._params.component_list,
                                   rule=rule_delta_eq)

        def rule_dadT(b, p):
            return -((b._params.gas_const/2)*sqrt(b.omegaA) *
                     sum(sum(b.mole_frac_phase[p, i] *
                             b.mole_frac_phase[p, j] *
                             (1-b._params.kappa[i, j]) *
                             (b.fw[j]*sqrt(b.a[i] *
                                           b._params.temperature_crit[j] /
                                           b._params.pressure_crit[j]) +
                              b.fw[i]*sqrt(b.a[j] *
                                           b._params.temperature_crit[i] /
                                           b._params.pressure_crit[i]))
                             for j in b._params.component_list)
                         for i in b._params.component_list) /
                     sqrt(b.temperature))
        blk.dadT = Expression(blk._params.phase_list, rule=rule_dadT)

        blk._ext_func_param = Param(default=blk._params.cubic_type.value)

        def rule_compress_fact(b, p):
            if p == "Vap":
                return b.proc_Z_vap(b._ext_func_param, b.A[p], b.B[p])
            else:
                return b.proc_Z_liq(b._ext_func_param, b.A[p], b.B[p])
        blk.compress_fact = Expression(blk._params.phase_list,
                                       rule=rule_compress_fact)

        def rule_compress_fact_eq(b, p):
            if p == "Vap":
                return b.proc_Z_vap(b._ext_func_param, b._A_eq[p], b._B_eq[p])
            else:
                return b.proc_Z_liq(b._ext_func_param, b._A_eq[p], b._B_eq[p])
        blk._compress_fact_eq = Expression(blk._params.phase_list,
                                           rule=rule_compress_fact_eq)

    def _vol_mol_cubic(b, p):
        return (b.pressure*b.vol_mol_phase[p] ==
                b.compress_fact[p]*b._params.gas_const*b.temperature)

    def _dens_mol_cubic(b, p):
        return b.pressure == (b.dens_mol_phase[p]*b.compress_fact[p] *
                              b._params.gas_const*b.temperature)

    def _dens_mass_cubic(b, p):
        return b.dens_mass_phase[p] == b.dens_mol_phase[p]*b.mw_phase[p]

    def _fug_cubic(b, p, j):
        return b.mole_frac_phase[p, j]*b.pressure*b.fug_coeff_phase[p, j]

    def _fug_coeff_cubic(b, p, j):
        return exp((b.b[j]/b.bm[p]*(b.compress_fact[p]-1) *
                    (b.B[p]*b.EoS_p) -
                    log(b.compress_fact[p]-b.B[p]) *
                    (b.B[p]*b.EoS_p) +
                    b.A[p]*(b.b[j]/b.bm[p] - b.delta[p, j]) *
                    log((2*b.compress_fact[p] +
                         b.B[p]*(b.EoS_u + b.EoS_p)) /
                        (2*b.compress_fact[p] +
                         b.B[p]*(b.EoS_u-b.EoS_p)))) /
                   (b.B[p]*b.EoS_p))

    def _fug_cubic_eq(b, p, j):
        return (b.mole_frac_phase[p, j]*b.pressure *
                exp((b.b[j]/b.bm[p]*(b._compress_fact_eq[p]-1) *
                     (b._B_eq[p]*b.EoS_p) -
                     log(b._compress_fact_eq[p]-b._B_eq[p]) *
                     (b._B_eq[p]*b.EoS_p) +
                     b._A_eq[p]*(b.b[j]/b.bm[p] - b._delta_eq[p, j]) *
                     log((2*b._compress_fact_eq[p] +
                          b._B_eq[p]*(b.EoS_u + b.EoS_p)) /
                         (2*b._compress_fact_eq[p] +
                          b._B_eq[p]*(b.EoS_u-b.EoS_p)))) /
                    (b._B_eq[p]*b.EoS_p)))

    def _enth_mol_cubic(b, p):
        return (((b.temperature*b.dadT[p] - b.am[p]) *
                 log((2*b.compress_fact[p] + b.B[p]*(b.EoS_u + b.EoS_p)) /
                     (2*b.compress_fact[p] + b.B[p]*(b.EoS_u - b.EoS_p))) +
                 b._params.gas_const*b.temperature*(b.compress_fact[p]-1) *
                 b.bm[p]*b.EoS_p) / (b.bm[p]*b.EoS_p) + b._enth_mol_ig(p))

    def _enth_mol_ig(b, p):
        return sum(b.mole_frac_phase[p, j]*b._enth_mol_comp_ig(j)
                   for j in b._params.component_list)

    def _entr_mol_cubic(b, p):
        return ((b._params.gas_const*log((b.compress_fact[p]-b.B[p]) /
                                         b.compress_fact[p])*b.bm[p]*b.EoS_p +
                 b._params.gas_const*log(b.compress_fact[p] *
                                         b._params.pressure_ref/b.pressure) *
                 b.bm[p]*b.EoS_p +
                 b.dadT[p]*log((2*b.compress_fact[p] +
                                b.B[p]*(b.EoS_u + b.EoS_p)) /
                               (2*b.compress_fact[p] +
                                b.B[p]*(b.EoS_u-b.EoS_p)))) /
                (b.bm[p]*b.EoS_p) + b._entr_mol_ig(p))

    def _entr_mol_ig(b, p):
        return sum(b.mole_frac_phase[p, j]*b._entr_mol_comp_ig(j)
                   for j in b._params.component_list)

# -----------------------------------------------------------------------------
# Pure component properties
    # TODO : Implementation question. Should this be here as part of a list of
    # options, or part of the ParameterBlock as a user-specified method?
    def _enth_mol_comp_ig(b, j):
        return (
            (b._params.cp_ig[j, "5"]/5) *
            (b.temperature**5-b._params.temperature_ref**5) +
            (b._params.cp_ig[j, "4"]/4) *
            (b.temperature**4-b._params.temperature_ref**4) +
            (b._params.cp_ig[j, "3"]/3) *
            (b.temperature**3-b._params.temperature_ref**3) +
            (b._params.cp_ig[j, "2"]/2) *
            (b.temperature**2-b._params.temperature_ref**2) +
            b._params.cp_ig[j, "1"] *
            (b.temperature-b._params.temperature_ref))

    # TODO : Implementation question. Should this be here as part of a list of
    # options, or part of the ParameterBlock as a user-specified method?
    def _entr_mol_comp_ig(b, j):
        return ((b._params.cp_ig[j, '4']/3) *
                (b.temperature**3-b._params.temperature_ref**3) +
                (b._params.cp_ig[j, '3']/2) *
                (b.temperature**2-b._params.temperature_ref**2) +
                b._params.cp_ig[j, '2'] *
                (b.temperature-b._params.temperature_ref) +
                b._params.cp_ig[j, '1'] *
                log(b.temperature/b._params.temperature_ref))
