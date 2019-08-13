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
Framework for generic property packages
"""
# Import Python libraries
import math
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Expression,
                           SolverFactory,
                           TerminationCondition,
                           value,
                           Var)
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.util.initialization import solve_indexed_blocks
from idaes.core.util.model_statistics import degrees_of_freedom

# Import sub-model libraries
import idaes.property_models.core.generic.state_methods as state_methods


# Set up logger
_log = logging.getLogger(__name__)


# TODO: Need clean-up methods for all methods to work with Pyomo DAE
# TODO: Sub-libraries for common forms of derived properties (e.g. Gibbs energy, mixture and phase properties)

class GenericParameterData(PhysicalParameterBlock):
    """
    General Property Parameter Block Class
    """
    CONFIG = PhysicalParameterBlock.CONFIG()

    # General options
    CONFIG.declare("state_definition", ConfigValue(
        default=state_methods.FPTx,
        description="Choice of State Variable",
        doc="""Flag indicating the set of state variables to use for property
        package. Values should be a valid Python method which creates the
        required state variables. Default = state_methods.FPHx"""))
    CONFIG.declare("state_bounds", ConfigValue(
        domain=dict,
        description="Bounds for state variables",
        doc="""A dict containing bounds to use for state variables."""))

    CONFIG.declare("phase_equilibrium_formulation", ConfigValue(
        default=None,
        description="Phase equilibrium formulation to use",
        doc="""Flag indicating what formulation to use for calcuating phase
        equilibrium. Value should be a valid Python method or None. Default =
        None, indicating no phase equilibrium will occur."""))

    CONFIG.declare("bubble_temperature", ConfigValue(
        description="Method to use to calculate bubble temperature",
        doc="""Flag indicating what formulation to use for calculating bubble
        temperature. Value should be a valid Python method."""))
    CONFIG.declare("dew_temperature", ConfigValue(
        description="Method to use to calculate dew temperature",
        doc="""Flag indicating what formulation to use for calculating dew
        temperature. Value should be a valid Python method."""))
    CONFIG.declare("bubble_pressure", ConfigValue(
        description="Method to use to calculate bubble pressure",
        doc="""Flag indicating what formulation to use for calculating bubble
        pressure. Value should be a valid Python method."""))
    CONFIG.declare("dew_pressure", ConfigValue(
        description="Method to use to calculate dew pressure",
        doc="""Flag indicating what formulation to use for calculating dew
        pressure. Value should be a valid Python method."""))

    # Equation of state options
    CONFIG.declare("equation_of_state", ConfigValue(
        domain=dict,
        description="Equation of state for each phase",
        doc="""Flag containing a dict indicating the equation of state for
        each phase. Value should be a dict with keys for each valid phase and
        values being a valid Python module with the necessary methods for
        thedesired equation of state."""))

    # Pure component property options
    CONFIG.declare("dens_mol_liq", ConfigValue(
        description="Method to use to calculate liquid phase molar density",
        doc="""Flag indicating what method to use when calcuating liquid phase
        molar density."""))
    CONFIG.declare("enth_mol_liq", ConfigValue(
        description="Method to calculate liquid component molar enthalpies",
        doc="""Flag indicating what method to use when calculating liquid phase
        component molar enthalpies."""))
    CONFIG.declare("enth_mol_vap", ConfigValue(
        description="Method to calculate vapor component molar enthalpies",
        doc="""Flag indicating what method to use when calculating vapor phase
        component molar enthalpies."""))
    CONFIG.declare("entr_mol_liq", ConfigValue(
        description="Method to calculate liquid component molar entropies",
        doc="""Flag indicating what method to use when calculating liquid phase
        component molar entropies."""))
    CONFIG.declare("entr_mol_vap", ConfigValue(
        description="Method to calculate vapor component molar entropies",
        doc="""Flag indicating what method to use when calculating vapor phase
        component molar entropies."""))
    CONFIG.declare("pressure_sat", ConfigValue(
        description="Method to use to calculate saturation pressure",
        doc="""Flag indicating what method to use when calcuating saturation
        pressure. Value should be a valid Python method which takes two
        arguments: temperature and a component name."""))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(GenericParameterData, self).build()

        self.state_block_class = GenericStateBlock

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        # TODO : Need to fix to have methods for things that may or may not be
        # created by state var methods
        obj.add_properties(
            {'flow_mol': {'method': None, 'units': 'mol/s'},
             'mole_frac': {'method': None, 'units': 'none'},
             'mole_frac_phase': {'method': None, 'units': 'none'},
             'phase_frac': {'method': None, 'units': 'none'},
             'temperature': {'method': None, 'units': 'K'},
             'pressure': {'method': None, 'units': 'Pa'},
             'flow_mol_phase': {'method': None, 'units': 'mol/s'},
             'dens_mass': {'method': '_dens_mass', 'units': 'kg/m^3'},
             'dens_mass_phase': {'method': '_dens_mass_phase',
                                 'units': 'kg/m^3'},
             'dens_mol': {'method': '_dens_mol', 'units': 'mol/m^3'},
             'dens_mol_phase': {'method': '_dens_mol_phase',
                                'units': 'mol/m^3'},
             'enth_mol': {'method': '_enth_mol', 'units': 'J/mol'},
             'enth_mol_phase': {'method': '_enth_mol_phase', 'units': 'J/mol'},
             'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp',
                                     'units': 'J/mol'},
             'entr_mol': {'method': '_entr_mol', 'units': 'J/mol.K'},
             'entr_mol_phase': {'method': '_entr_mol_phase',
                                'units': 'J/mol.K'},
             'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp',
                                     'units': 'J/mol.K'},
             'fug': {'method': '_fug', 'units': 'Pa'},
             'fug_coeff': {'method': '_fug_coeff', 'units': '-'},
             'gibbs_mol': {'method': '_gibbs_mol', 'units': 'J/mol'},
             'gibbs_mol_phase': {'method': '_gibbs_mol_phase',
                                 'units': 'J/mol'},
             'gibbs_mol_phase_comp': {'method': '_gibbs_mol_phase_comp',
                                      'units': 'J/mol'},
             'mw': {'method': '_mw', 'units': 'kg/mol'},
             'mw_phase': {'method': '_mw_phase', 'units': 'kg/mol'},
             'pressure_bubble': {'method': '_pressure_bubble', 'units': 'Pa'},
             'pressure_dew': {'method': '_pressure_dew', 'units': 'Pa'},
             'pressure_sat': {'method': '_pressure_sat', 'units': 'Pa'},
             'temperature_bubble': {'method': '_temperature_bubble',
                                    'units': 'K'},
             'temperature_dew': {'method': '_temperature_dew', 'units': 'K'}})

        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _GenericStateBlock(StateBlock):
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

                    if f/df > 20:
                        Tbub1 = Tbub0 - 20
                    elif f/df < -20:
                        Tbub1 = Tbub0 + 20
                    else:
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

                    if f/df > 20:
                        Tdew1 = Tdew0 - 20
                    elif f/df < -20:
                        Tdew1 = Tdew0 + 20
                    else:
                        Tdew1 = Tdew0 - f/df

                    err = abs(Tdew1 - Tdew0)
                    Tdew0 = Tdew1
                    counter += 1

                blk[k].temperature_dew.value = Tdew0

                for j in blk[k]._params.component_list:
                    blk[k]._mole_frac_tdew[j].value = value(
                            blk[k].mole_frac[j]*blk[k].pressure /
                            antoine_P(blk[k], j, Tdew0))

        # Bubble pressure initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_pbub"):
                blk[k].pressure_bubble.value = value(
                        sum(blk[k].mole_frac[j] *
                            antoine_P(blk[k], j, blk[k].temperature)
                            for j in blk[k]._params.component_list))

                for j in blk[k]._params.component_list:
                    blk[k]._mole_frac_pbub[j].value = value(
                            blk[k].mole_frac[j] *
                            antoine_P(blk[k], j, blk[k].temperature) /
                            blk[k].pressure_bubble)

                blk[k].pressure_bubble.display()
                blk[k]._mole_frac_pbub.display()

        # Dew pressure initialization
        for k in blk.keys():
            if hasattr(blk[k], "_mole_frac_pdew"):
                blk[k].pressure_dew.value = value(
                        sum(1/(blk[k].mole_frac[j] /
                               antoine_P(blk[k], j, blk[k].temperature))
                            for j in blk[k]._params.component_list))

                for j in blk[k]._params.component_list:
                    blk[k]._mole_frac_pdew[j].value = value(
                            blk[k].mole_frac[j]*blk[k].pressure_bubble /
                            antoine_P(blk[k], j, blk[k].temperature))

        # Solve bubble and dew point constraints
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Deactivate all property constraints
                if c.local_name not in ("eq_pressure_dew",
                                        "eq_pressure_bubble",
                                        "eq_temperature_dew",
                                        "eq_temperature_bubble",
                                        "_sum_mole_frac_tbub",
                                        "_sum_mole_frac_tdew",
                                        "_sum_mole_frac_pbub",
                                        "_sum_mole_frac_pdew"):
                    c.deactivate()

        results = solve_indexed_blocks(opt, [blk], tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Dew and bubble points initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Dew and bubble points initialization for "
                             "{} failed".format(blk.name))

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
        # TODO : This will need to be generalised more when we move to a
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
                if blk[k].temperature.value > blk[k].temperature_dew.value:
                    # Pure vapour
                    blk[k].flow_mol_phase["Vap"].value = blk[k].flow_mol.value
                    blk[k].flow_mol_phase["Vap"].value = \
                        1e-5*blk[k].flow_mol.value

                    for j in blk[k]._params.component_list:
                        blk[k].mole_frac_phase['Vap', j].value = \
                            blk[k].mole_frac[j].value
                        blk[k].mole_frac_phase['Liq', j].value = \
                            blk[k]._mole_frac_tdew[j].value
                elif blk[k].temperature.value < \
                        blk[k].temperature_bubble.value:
                    # Pure liquid
                    blk[k].flow_mol_phase["Vap"].value = \
                        1e-5*blk[k].flow_mol.value
                    blk[k].flow_mol_phase["Vap"].value = blk[k].flow_mol.value

                    for j in blk[k]._params.component_list:
                        blk[k].mole_frac_phase['Vap', j].value = \
                            blk[k]._mole_frac_tbub[j].value
                        blk[k].mole_frac_phase['Liq', j].value = \
                            blk[k].mole_frac[j].value
                else:
                    # Two-phase
                    # TODO : Try to find some better guesses than this
                    blk[k].flow_mol_phase["Vap"].value = \
                        0.5*blk[k].flow_mol.value
                    blk[k].flow_mol_phase["Vap"].value = \
                        0.5*blk[k].flow_mol.value

                    for j in blk[k]._params.component_list:
                        blk[k].mole_frac_phase['Vap', j].value = \
                            blk[k].mole_frac[j].value
                        blk[k].mole_frac_phase['Liq', j].value = \
                            blk[k].mole_frac[j].value

        # ---------------------------------------------------------------------
        # Solve phase equilibrium constraints
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Activate equilibrium constraints
                if c.local_name in ("total_flow_balance",
                                    "component_flow_balances",
                                    "equilibrium_constraint",
                                    "sum_mole_frac",
                                    "_t1_constraint",
                                    "_teq_constraint"):
                    c.activate()

        results = solve_indexed_blocks(opt, [blk], tee=stee)

        if outlvl > 0:
            if results.solver.termination_condition \
                    == TerminationCondition.optimal:
                _log.info("Phase equilibrium initialization for "
                          "{} completed".format(blk.name))
            else:
                _log.warning("Phase equilibrium initialization for "
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

@declare_process_block_class("GenericStateBlock",
                             block_class=_GenericStateBlock)
class GenericStateBlockData(StateBlockData):
    """A modular, general purpose property package."""

    def build(self):
        """Callable method for Block construction."""
        super(GenericStateBlockData, self).build()

        # Add state vairables and assoicated methods
        self._params.config.state_definition(self)

        # Create common components for each property package
        for p in self._params.phase_list:
            self._params.config.equation_of_state[p].common(self)

        # Add phase equilibrium constraints if necessary
        if (self._params.config.phase_equilibrium_formulation is not None and
                (not self.config.defined_state or self.always_flash)):
            self._params.config.phase_equilibrium_formulation(self)

    # -------------------------------------------------------------------------
    # Bubble and Dew Points
    def _temperature_bubble(b):
        b.temperature_bubble = Var(
                doc="Bubble point temperature",
                bounds=(b.temperature.lb, b.temperature.ub))

        b._mole_frac_tbub = Var(
                b._params.component_list,
                initialize=1/len(b._params.component_list),
                bounds=(0, None),
                doc="Vapor mole fractions at bubble temperature")

        b._params.config.bubble_temperature(b)

    def _temperature_dew(b):
        b.temperature_dew = Var(
                doc="Dew point temperature",
                bounds=(b.temperature.lb, b.temperature.ub))

        b._mole_frac_tdew = Var(
                b._params.component_list,
                initialize=1/len(b._params.component_list),
                bounds=(0, None),
                doc="Liquid mole fractions at dew temperature")

        b._params.config.dew_temperature(b)

    def _pressure_bubble(b):
        b.pressure_bubble = Var(
                doc="Bubble point pressure",
                bounds=(b.pressure.lb, b.pressure.ub))

        b._mole_frac_pbub = Var(
                b._params.component_list,
                initialize=1/len(b._params.component_list),
                bounds=(0, None),
                doc="Vapor mole fractions at bubble pressure")

        b._params.config.bubble_pressure(b)

    def _pressure_dew(b):
        b.pressure_dew = Var(
                doc="Dew point pressure",
                bounds=(b.pressure.lb, b.pressure.ub))

        b._mole_frac_pdew = Var(
                b._params.component_list,
                initialize=1/len(b._params.component_list),
                bounds=(0, None),
                doc="Liquid mole fractions at dew pressure")

        b._params.config.dew_pressure(b)

    # -------------------------------------------------------------------------
    # Property Methods
    def _dens_mass(self):
        def rule_dens_mass(b):
            return sum(b.dens_mass_phase[p]*b.phase_frac[p]
                       for p in b._params.phase_list)
        self.dens_mass = Expression(
                doc="Mixture mass density",
                rule=rule_dens_mass)

    def _dens_mass_phase(self):
        def rule_dens_mass_phase(b, p):
            return b._params.config.equation_of_state[p].dens_mass(b, p)
        self.dens_mass_phase = Expression(
                self._params.phase_list,
                doc="Mass density of each phase",
                rule=rule_dens_mass_phase)

    def _dens_mol(self):
        def rule_dens_mol(b):
            return sum(b.dens_mol_phase[p]*b.phase_frac[p]
                       for p in b._params.phase_list)
        self.dens_mol = Expression(
                doc="Mixture molar density",
                rule=rule_dens_mol)

    def _dens_mol_phase(self):
        def rule_dens_mol_phase(b, p):
            return b._params.config.equation_of_state[p].dens_mol(b, p)
        self.dens_mol_phase = Expression(
                self._params.phase_list,
                doc="Molar density of each phase",
                rule=rule_dens_mol_phase)

    def _enth_mol(self):
        self.enth_mol = Var(doc="Mixture molar enthalpy")

        def rule_enth_mol(b):
            return b.enth_mol == sum(
                    sum(b.mole_frac_phase[p, j] *
                        b.enth_mol_phase_comp[p, j]
                        for j in b._params.component_list)
                    for p in b._params.phase_list)
        self.eq_enth_mol_phase = Constraint(rule=rule_enth_mol)

    def _enth_mol_phase(self):
        def rule_enth_mol_phase(b, p):
            return sum(b.mole_frac_phase[p, j] *
                       b.enth_mol_phase_comp[p, j]
                       for j in b._params.component_list)
        self.enth_mol_phase = Expression(self._params.phase_list,
                                         rule=rule_enth_mol_phase)

    def _enth_mol_phase_comp(self):
        def rule_enth_mol_phase_comp(b, p, j):
            return b._params.config.equation_of_state[p].enth_mol_comp(b, p, j)
        self.enth_mol_phase_comp = Expression(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_enth_mol_phase_comp)

    def _entr_mol(self):
        self.entr_mol = Var(doc="Mixture molar entropy")

        def rule_entr_mol(b):
            return b.entr_mol == sum(
                    sum(b.mole_frac_phase[p, j] *
                        b.entr_mol_phase_comp[p, j]
                        for j in b._params.component_list)
                    for p in b._params.phase_list)
        self.eq_entr_mol_phase = Constraint(rule=rule_entr_mol)

    def _entr_mol_phase(self):
        def rule_entr_mol_phase(b, p):
            return sum(b.mole_frac_phase[p, j] *
                       b.entr_mol_phase_comp[p, j]
                       for j in b._params.component_list)
        self.entr_mol_phase = Expression(self._params.phase_list,
                                         rule=rule_entr_mol_phase)

    def _entr_mol_phase_comp(self):
        def rule_entr_mol_phase_comp(b, p, j):
            return b._params.config.equation_of_state[p].entr_mol_comp(b, p, j)
        self.entr_mol_phase_comp = Expression(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_entr_mol_phase_comp)

    def _fug(self):
        def rule_fug(b, p, j):
            return b._params.config.equation_of_state[p].fugacity(b, p, j)
        self.fug = Expression(self._params.phase_list,
                              self._params.component_list,
                              rule=rule_fug)

    def _fug_coeff(self):
        def rule_fug_coeff(b, p, j):
            return b._params.config.equation_of_state[p].fug_coeff(b, p, j)
        self.fug_coeff = Expression(self._params.phase_list,
                                    self._params.component_list,
                                    rule=rule_fug_coeff)

    def _gibbs_mol(self):
        self.gibbs_mol = Var(doc="Mixture molar Gibbs energy")

        def rule_gibbs_mol(b):
            return b.gibbs_mol == sum(
                    sum(b.mole_frac_phase[p, j] *
                        b.gibbs_mol_phase_comp[p, j]
                        for j in b._params.component_list)
                    for p in b._params.phase_list)
        self.eq_gibbs_mol_phase = Constraint(rule=rule_gibbs_mol)

    def _gibbs_mol_phase(self):
        def rule_gibbs_mol_phase(b, p):
            return sum(b.mole_frac_phase[p, j] *
                       b.gibbs_mol_phase_comp[p, j]
                       for j in b._params.component_list)
        self.gibbs_mol_phase = Expression(self._params.phase_list,
                                          rule=rule_gibbs_mol_phase)

    def _gibbs_mol_phase_comp(self):
        def rule_gibbs_mol_phase_comp(b, p, j):
            return (b.enth_mol_phase_comp[p, j] -
                    (b.entr_mol_phase_comp[p, j] -
                     b._params.config.equation_of_state[p]
                     .entr_mol_comp_ref(b, p, j))*b.temperature)
        self.gibbs_mol_phase_comp = Expression(
            self._params.phase_list,
            self._params.component_list,
            rule=rule_gibbs_mol_phase_comp)

    def _mw(self):
        self.mw = Expression(
                doc="Average molecular weight",
                expr=sum(self.phase_frac[p] *
                         sum(self.mole_frac_phase[p, j]*self._params.mw_comp[j]
                             for j in self._params.component_list)
                         for p in self._params.phase_list))

    def _mw_phase(self):
        def rule_mw_phase(b, p):
            return sum(b.mole_frac_phase[p, j]*b._params.mw_comp[j]
                       for j in b._params.component_list)
        self.mw_phase = Expression(
                self._params.phase_list,
                doc="Average molecular weight of each phase",
                rule=rule_mw_phase)
