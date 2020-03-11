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
import types

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Expression,
                           Set,
                           SolverFactory,
                           value,
                           Var)
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_activated_constraints)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyPackageError)
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


# TODO: Need clean-up methods for all methods to work with Pyomo DAE
# TODO: Need way to dynamically determine units of measurement....
class GenericPropertyPackageError(PropertyPackageError):
    # Error message for when a property is called for but no option provided
    def __init__(self, block, prop):
        self.prop = prop
        self.block = block

    def __str__(self):
        return f"Generic Property Package instance {self.block} called for " \
               f"{self.prop}, but was not provided with a method " \
               f"for this property. Please add a method for this property " \
               f"in the property parameter configuration."


def get_method(self, config_arg):
    """
    Method to inspect configuration argument and return the user-defined
    construction method associated with it.

    This method checks whether the value provided is a method or a
    module. If the value is a module, it looks in the module for a method
    with the same name as the config argument and returns this. If the
    value is a method, the method is returned. If the value is neither a
    module or a method, an ConfigurationError is raised.

    Args:
        config_arg : the configuration argument to look up

    Returns:
        A callable method or a ConfigurationError
    """
    try:
        c_arg = getattr(self.params.config, config_arg)
    except AttributeError:
        raise AttributeError("{} Generic Property Package called for invalid "
                             "configuration option {}. Please contact the "
                             "developer of the property package."
                             .format(self.name, config_arg))

    if c_arg is None:
        raise GenericPropertyPackageError(self, c_arg)

    if isinstance(c_arg, types.ModuleType):
        return getattr(c_arg, config_arg)
    elif callable(c_arg):
        return c_arg
    else:
        raise ConfigurationError(
                "{} Generic Property Package received invalid value "
                "for argumnet {}. Value must be either a module or a "
                "method".format(self.name, config_arg))


class GenericParameterData(PhysicalParameterBlock):
    """
    General Property Parameter Block Class
    """
    CONFIG = PhysicalParameterBlock.CONFIG()

    # General options
    CONFIG.declare("component_list", ConfigValue(
        description="List of components in material",
        doc="""A list of names for the components of interest in the mixture.
        """))
    CONFIG.declare("phase_list", ConfigValue(
        description="List of phases of interest",
        doc="""A list of phases of interest in the mixture for the property
        package."""))
    CONFIG.declare("phase_component_list", ConfigValue(
        description="List of components in each phase",
        doc="""A dict of component_lists for each phase. Keys should correspond
        to members of phase_list, and each value should be a list of component
        names."""))

    CONFIG.declare("state_definition", ConfigValue(
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
        doc="""Flag indicating what formulation to use for calculating phase
        equilibrium. Value should be a valid Python method or None. Default =
        None, indicating no phase equilibrium will occur."""))
    CONFIG.declare("phase_equilibrium_dict", ConfigValue(
        default=None,
        description="Phase equilibrium reactions to be modeled",
        doc="""Dict describing what phase equilibrium reactions should be
        included in the property package. Keys should be names from
        component_list, and values should be a 2-tuple of phases from
        phase_list which should be in equilibrium."""))

    CONFIG.declare("temperature_bubble", ConfigValue(
        description="Method to use to calculate bubble temperature",
        doc="""Flag indicating what formulation to use for calculating bubble
        temperature. Value should be a valid Python method."""))
    CONFIG.declare("temperature_dew", ConfigValue(
        description="Method to use to calculate dew temperature",
        doc="""Flag indicating what formulation to use for calculating dew
        temperature. Value should be a valid Python method."""))
    CONFIG.declare("pressure_bubble", ConfigValue(
        description="Method to use to calculate bubble pressure",
        doc="""Flag indicating what formulation to use for calculating bubble
        pressure. Value should be a valid Python method."""))
    CONFIG.declare("pressure_dew", ConfigValue(
        description="Method to use to calculate dew pressure",
        doc="""Flag indicating what formulation to use for calculating dew
        pressure. Value should be a valid Python method."""))

    # Equation of state options
    CONFIG.declare("equation_of_state", ConfigValue(
        description="Equation of state for each phase",
        doc="""Flag containing a dict indicating the equation of state for
        each phase. Value should be a dict with keys for each valid phase and
        values being a valid Python module with the necessary methods for
        the desired equation of state."""))

    # Pure component property options
    CONFIG.declare("dens_mol_liq_comp", ConfigValue(
        description="Method to use to calculate liquid phase molar density",
        doc="""Flag indicating what method to use when calculating liquid phase
        molar density."""))
    CONFIG.declare("enth_mol_liq_comp", ConfigValue(
        description="Method to calculate liquid component molar enthalpies",
        doc="""Flag indicating what method to use when calculating liquid phase
        component molar enthalpies."""))
    CONFIG.declare("enth_mol_ig_comp", ConfigValue(
        description="Method to calculate ideal gas component molar enthalpies",
        doc="""Flag indicating what method to use when calculating ideal gas
        phase component molar enthalpies."""))
    CONFIG.declare("entr_mol_liq_comp", ConfigValue(
        description="Method to calculate liquid component molar entropies",
        doc="""Flag indicating what method to use when calculating liquid phase
        component molar entropies."""))
    CONFIG.declare("entr_mol_ig_comp", ConfigValue(
        description="Method to calculate ideal gas component molar entropies",
        doc="""Flag indicating what method to use when calculating ideal gas
        phase component molar entropies."""))
    CONFIG.declare("pressure_sat_comp", ConfigValue(
        description="Method to use to calculate saturation pressure",
        doc="""Flag indicating what method to use when calculating saturation
        pressure. Value should be a valid Python method which takes two
        arguments: temperature and a component name."""))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        # Call super.build() to initialize Block
        super(GenericParameterData, self).build()

        # Call configure method to set construction arguments
        self.configure()

        # Build core components
        self.state_block_class = GenericStateBlock

        if self.config.phase_component_list is not None:
            # Check if phase list provided and cross-validate
            if self.config.phase_list is not None:
                # Phase list provided, cross-validate
                for p in self.config.phase_component_list:
                    if p not in self.config.phase_list:
                        raise ConfigurationError(
                            "{} mismatch between phase_list and "
                            "phase_component_list. Phase {} appears in "
                            "phase_component_list but not phase_list."
                            .format(self.name, p))

                for p in self.config.phase_list:
                    if p not in self.config.phase_component_list:
                        raise ConfigurationError(
                            "{} mismatch between phase_list and "
                            "phase_component_list. Phase {} appears in "
                            "phase_list but not phase_component_list."
                            .format(self.name, p))

                # Build phase_list if cross-validation passes
                self.phase_list = Set(initialize=self.config.phase_list,
                                      ordered=True)
            else:
                # No phase_list provided, build from phase_component_list
                self.phase_list = Set(
                        initialize=[p for p in
                                    self.config.phase_component_list],
                        ordered=True)

            if self.config.component_list is not None:
                # Component list provided, cross-validate
                for p in self.config.phase_component_list:
                    for j in self.config.phase_component_list[p]:
                        if j not in self.config.component_list:
                            raise ConfigurationError(
                                "{} mismatch between component_list and "
                                "phase_component_list. Component {} appears in"
                                " phase_component_list but not component_list."
                                .format(self.name, j))
                for j in self.config.component_list:
                    xcheck = False
                    for p in self.config.phase_component_list:
                        if j in self.config.phase_component_list[p]:
                            xcheck = True
                            break
                    if not xcheck:
                        raise ConfigurationError(
                                "{} mismatch between component_list and "
                                "phase_component_list. Component {} appears in"
                                " component_list but not phase_component_list."
                                .format(self.name, j))

                # Build component_list if cross-validation passes
                self.component_list = Set(
                        initialize=self.config.component_list,
                        ordered=True)
            else:
                # No component_list provided, build from phase_component_list
                c_list = []
                for p in self.config.phase_component_list:
                    for j in self.config.phase_component_list[p]:
                        if j not in c_list:
                            c_list.append(j)
                self.component_list = Set(initialize=c_list, ordered=True)

            # All validation passed, build phase_component_set
            pc_set = []
            for p in self.config.phase_component_list:
                for j in self.config.phase_component_list[p]:
                    pc_set.append((p, j))
            self._phase_component_set = Set(initialize=pc_set, ordered=True)

        elif (self.config.phase_list is not None and
              self.config.component_list is not None):
            # Have phase and component lists
            # Assume all components in all phases
            self.phase_list = Set(initialize=self.config.phase_list,
                                  ordered=True)
            self.component_list = Set(initialize=self.config.component_list,
                                      ordered=True)

            # Create phase-component set
            pc_set = []
            for p in self.phase_list:
                for j in self.component_list:
                    pc_set.append((p, j))
            self._phase_component_set = Set(initialize=pc_set, ordered=True)
        else:
            # User has not provided sufficient information.
            if self.config.component_list is None:
                raise ConfigurationError(
                        "{} Generic Property Package was not provided with a "
                        "component_list or a phase_component_list. Users must "
                        "provide either at least one of these arguments"
                        .format(self.name))

            if self.config.phase_list is None:
                raise ConfigurationError(
                        "{} Generic Property Package was not provided with a "
                        "phase_list or a phase_component_list. Users must "
                        "provide either at least one of these arguments"
                        .format(self.name))

        # Validate state definition
        if self.config.state_definition is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "state_definition configuration argument. Please fix "
                    "your property parameter definition to include this "
                    "configuration argument.".format(self.name))

        # Validate equation of state
        if self.config.equation_of_state is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with an "
                    "equation_of_state configuration argument. Please fix "
                    "your property parameter definition to include this "
                    "configuration argument.".format(self.name))
        if not isinstance(self.config.equation_of_state, dict):
            raise ConfigurationError(
                    "{} Generic Property Package was provided with an invalid "
                    "equation_of_state configuration argument. Argument must "
                    "be a dict with phases as keys.".format(self.name))
        if len(self.config.equation_of_state) != len(self.phase_list):
            raise ConfigurationError(
                    "{} Generic Property Package was provided with an invalid "
                    "equation_of_state configuration argument. A value must "
                    "be present for each phase.".format(self.name))
        for p in self.config.equation_of_state:
            if p not in self.phase_list:
                raise ConfigurationError(
                    "{} Generic Property Package unrecognised phase {} in "
                    "equation_of_state configuration argument. Keys must be "
                    "valid phases.".format(self.name, p))

        # Validate that user provided either both a phase equilibrium
        # formulation and a dict of phase equilibria or neither
        if ((self.config.phase_equilibrium_formulation is not None) ^
                (self.config.phase_equilibrium_dict is not None)):
            raise ConfigurationError(
                    "{} Generic Property Package provided with only one of "
                    "phase_equilibrium_formulation and phase_equilibrium_dict."
                    " Either both of these arguments need to be provided or "
                    "neither.".format(self.name))

        # Validate and build phase equilibrium list
        if self.config.phase_equilibrium_dict is not None:
            if not isinstance(self.config.phase_equilibrium_dict, dict):
                raise ConfigurationError(
                    "{} Generic Property Package provided with invalid "
                    "phase_equilibrium_dict - value must be a dict. "
                    "Please see the documentation for the correct form."
                    .format(self.name))
            # Validate phase_equilibrium_dict
            for v in self.config.phase_equilibrium_dict.values():
                if not (isinstance(v, list) and len(v) == 2):
                    raise ConfigurationError(
                        "{} Generic Property Package provided with invalid "
                        "phase_equilibrium_dict, {}. Values in dict must be "
                        "lists containing 2 values.".format(self.name, v))
                if v[0] not in self.component_list:
                    raise ConfigurationError(
                        "{} Generic Property Package provided with invalid "
                        "phase_equilibrium_dict. First value in each list "
                        "must be a valid component, received {}."
                        .format(self.name, v[0]))
                if not (isinstance(v[1], tuple) and len(v[1]) == 2):
                    raise ConfigurationError(
                        "{} Generic Property Package provided with invalid "
                        "phase_equilibrium_dict. Second value in each list "
                        "must be a 2-tuple containing 2 valid phases, "
                        "received {}.".format(self.name, v[1]))
                for p in v[1]:
                    if p not in self.phase_list:
                        raise ConfigurationError(
                            "{} Generic Property Package provided with invalid"
                            " phase_equilibrium_dict. Unrecognised phase {} "
                            "in tuple {}".format(self.name, p, v[1]))

            self.phase_equilibrium_list = self.config.phase_equilibrium_dict

            pe_set = []
            for k in self.config.phase_equilibrium_dict.keys():
                pe_set.append(k)
            self.phase_equilibrium_idx = Set(initialize=pe_set,
                                             ordered=True)
        self.parameters()

    def configure(self):
        raise PropertyPackageError(
                "{} User defined property package failed to define a "
                "configure method. Please contact the developer of the "
                "property package with this error.".format(self.name))

    def parameters(self):
        raise PropertyPackageError(
                "{} User defined property package failed to define a "
                "parameters method. Please contact the developer of the "
                "property package with this error.".format(self.name))

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        # TODO : Need to fix to have methods for things that may or may not be
        # created by state var methods
        obj.add_properties(
            {'flow_mol': {'method': None, 'units': 'mol/s'},
             'mole_frac_comp': {'method': None, 'units': 'none'},
             'mole_frac_phase_comp': {'method': None, 'units': 'none'},
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
             'fug_phase_comp': {'method': '_fug_phase_comp', 'units': 'Pa'},
             'fug_coeff_phase_comp': {'method': '_fug_coeff_phase_comp',
                                      'units': '-'},
             'gibbs_mol': {'method': '_gibbs_mol', 'units': 'J/mol'},
             'gibbs_mol_phase': {'method': '_gibbs_mol_phase',
                                 'units': 'J/mol'},
             'gibbs_mol_phase_comp': {'method': '_gibbs_mol_phase_comp',
                                      'units': 'J/mol'},
             'mw': {'method': '_mw', 'units': 'kg/mol'},
             'mw_phase': {'method': '_mw_phase', 'units': 'kg/mol'},
             'pressure_bubble': {'method': '_pressure_bubble', 'units': 'Pa'},
             'pressure_dew': {'method': '_pressure_dew', 'units': 'Pa'},
             'pressure_sat_comp': {'method': '_pressure_sat_comp',
                                   'units': 'Pa'},
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

    def initialize(blk, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-8}):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : a dict of initial values for the state variables
                    defined by the property package.
            outlvl : sets output level of initialization routine
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
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                         should unfix any state variables fixed during
                         initialization (default=False).
                         - True - states variables are not unfixed, and
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
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="properties")

        init_log.info('Starting initialization')

        for k in blk.keys():
            # Deactivate the constraints specific for outlet block i.e.
            # when defined state is False
            if blk[k].config.defined_state is False:
                try:
                    blk[k].sum_mole_frac_out.deactivate()
                except AttributeError:
                    pass

        # Fix state variables if not already fixed
        if state_vars_fixed is False:
            flag_dict = fix_state_vars(blk, state_args)
            # Confirm DoF for sanity
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise BurntToast("Degrees of freedom were not zero "
                                     "after trying to fix state variables. "
                                     "Something broke in the generic property "
                                     "package code - please inform the IDAES "
                                     "developers.")
        else:
            # When state vars are fixed, check that DoF is 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception("State vars fixed but degrees of "
                                    "freedom for state block is not zero "
                                    "during initialization.")

        # Set solver options
        if optarg is None:
            sopt = {'tol': 1e-8}
        else:
            sopt = optarg

        opt = SolverFactory('ipopt')
        opt.options = sopt

        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations
        for k in blk.keys():
            # Bubble temperature initialization
            if hasattr(blk[k], "_mole_frac_tbub"):
                # Use lowest component critical temperature as starting point
                # Starting high and moving down generally works better,
                # as it under-predicts next step due to exponential form of
                # Psat.
                # Subtract 1 to avoid potential singularities at Tcrit
                Tbub0 = min(blk[k].params.temperature_crit_comp[j]
                            for j in blk[k].params.component_list) - 1

                err = 1
                counter = 0

                # Newton solver with step limiter to prevent overshoot
                # Tolerance only needs to be ~1e-1
                # Iteration limit of 30
                while err > 1e-1 and counter < 30:
                    f = value(sum(blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(blk[k], j, Tbub0) *
                                  blk[k].mole_frac_comp[j]
                                  for j in blk[k].params.component_list) -
                              blk[k].pressure)
                    df = value(sum(
                           blk[k].mole_frac_comp[j]*blk[k].params.config
                           .pressure_sat_comp.pressure_sat_comp_dT(
                                   blk[k], j, Tbub0)
                           for j in blk[k].params.component_list))

                    # Limit temperature step to avoid excessive overshoot
                    # Only limit positive steps due to non-linearity
                    if f/df < -50:
                        Tbub1 = Tbub0 + 50
                    else:
                        Tbub1 = Tbub0 - f/df

                    err = abs(Tbub1 - Tbub0)
                    Tbub0 = Tbub1
                    counter += 1

                blk[k].temperature_bubble.value = Tbub0

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_tbub[j].value = value(
                            blk[k].mole_frac_comp[j]*blk[k].pressure /
                            blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(blk[k], j, Tbub0))

            # Bubble temperature initialization
            if hasattr(blk[k], "_mole_frac_tdew"):
                if hasattr(blk[k], "_mole_frac_tbub"):
                    # If Tbub has been calculated above, use this as the
                    # starting point
                    Tdew0 = blk[k].temperature_bubble.value
                else:
                    # Otherwise, use lowest component critical temperature as
                    # starting point
                    # Subtract 1 to avoid potential singularities at Tcrit
                    Tdew0 = min(blk[k].params.temperature_crit_comp[j]
                                for j in blk[k].params.component_list) - 1

                err = 1
                counter = 0

                # Newton solver with step limiter to prevent overshoot
                # Tolerance only needs to be ~1e-1
                # Iteration limit of 30
                while err > 1e-1 and counter < 30:
                    f = value(blk[k].pressure *
                              sum(blk[k].mole_frac_comp[j] /
                                  blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(blk[k], j, Tdew0)
                                  for j in blk[k].params.component_list) - 1)
                    df = -value(
                            blk[k].pressure *
                            sum(blk[k].mole_frac_comp[j] /
                                blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(blk[k], j, Tdew0)**2 *
                                blk[k].params.config
                                .pressure_sat_comp.pressure_sat_comp_dT(
                                        blk[k], j, Tdew0)
                                for j in blk[k].params.component_list))

                    # Limit temperature step to avoid excessive overshoot
                    if f/df < -50:
                        Tdew1 = Tdew0 + 50
                    else:
                        Tdew1 = Tdew0 - f/df

                    err = abs(Tdew1 - Tdew0)
                    Tdew0 = Tdew1
                    counter += 1

                blk[k].temperature_dew.value = Tdew0

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_tdew[j].value = value(
                            blk[k].mole_frac_comp[j]*blk[k].pressure /
                            blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(blk[k], j, Tdew0))

            # Bubble pressure initialization
            if hasattr(blk[k], "_mole_frac_pbub"):
                blk[k].pressure_bubble.value = value(
                        sum(blk[k].mole_frac_comp[j] *
                            blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(
                                          blk[k], j, blk[k].temperature)
                            for j in blk[k].params.component_list))

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_pbub[j].value = value(
                        blk[k].mole_frac_comp[j] *
                        blk[k].params.config.pressure_sat_comp
                              .pressure_sat_comp(
                                      blk[k], j, blk[k].temperature) /
                        blk[k].pressure_bubble)

            # Dew pressure initialization
            if hasattr(blk[k], "_mole_frac_pdew"):
                blk[k].pressure_dew.value = value(
                        sum(1/(blk[k].mole_frac_comp[j] /
                               blk[k].params.config.pressure_sat_comp
                               .pressure_sat_comp(
                                       blk[k], j, blk[k].temperature))
                            for j in blk[k].params.component_list))

                for j in blk[k].params.component_list:
                    blk[k]._mole_frac_pdew[j].value = value(
                            blk[k].mole_frac_comp[j]*blk[k].pressure_bubble /
                            blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(
                                          blk[k], j, blk[k].temperature))

            # Solve bubble and dew point constraints
            for c in blk[k].component_objects(Constraint):
                # Deactivate all constraints not associated wtih bubble and dew
                # points
                if c.local_name not in ("eq_pressure_dew",
                                        "eq_pressure_bubble",
                                        "eq_temperature_dew",
                                        "eq_temperature_bubble",
                                        "_sum_mole_frac_tbub",
                                        "_sum_mole_frac_tdew",
                                        "_sum_mole_frac_pbub",
                                        "_sum_mole_frac_pdew"):
                    c.deactivate()

        # If StateBlock has active constraints (i.e. has bubble and/or dew
        # point calculations), solve the block to converge these
        n_cons = 0
        for k in blk:
            n_cons += number_activated_constraints(blk[k])
        if n_cons > 0:
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Dew and bubble point initialization: {}.".format(idaeslog.condition(res))
            )
        # ---------------------------------------------------------------------
        # If StateBlock is using a smooth VLE, calculate _T1 and _Teq
        eq_check = 0
        for k in blk.keys():
            if hasattr(blk[k], "_t1"):
                blk[k]._t1.value = max(blk[k].temperature.value,
                                       blk[k].temperature_bubble.value)
                blk[k]._teq.value = min(blk[k]._t1.value,
                                        blk[k].temperature_dew.value)

                eq_check += 1

        if eq_check > 0:
            init_log.info("Equilibrium temperature initialization completed.")

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        for k in blk.keys():
            blk[k].params.config.state_definition.state_initialization(blk[k])

        if outlvl > 0:
            init_log.info("State variable initialization completed.")

        # ---------------------------------------------------------------------
        if (blk[k].params.config.phase_equilibrium_formulation is not None and
                (not blk[k].config.defined_state or blk[k].always_flash)):
            blk[k].params.config.phase_equilibrium_formulation \
                .phase_equil_initialization(blk[k])

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt,[blk],tee=slc.tee)
            init_log.info(
                "Phase equilibrium initialization: {}.".format(
                    idaeslog.condition(res)
                )
            )

        # ---------------------------------------------------------------------
        # Initialize other properties
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Activate all constraints except flagged do_not_initialize
                if c.local_name not in (
                        blk[k].params.config
                        .state_definition.do_not_initialize):
                    c.activate()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
        init_log.info("Property initialization: {}.".format(
            idaeslog.condition(res))
        )

        # ---------------------------------------------------------------------
        # Return constraints to initial state
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                if c.local_name in (
                        blk[k].params.config
                        .state_definition.do_not_initialize):
                    c.activate()

        if state_vars_fixed is False:
            if hold_state is True:
                return flag_dict
            else:
                blk.release_state(flag_dict)

        init_log.info("Property package initialization: {}.".format(
            idaeslog.condition(res))
        )

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        '''
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of initialization routine
        '''
        revert_state_vars(blk, flags)
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info_high("State released.")

@declare_process_block_class("GenericStateBlock",
                             block_class=_GenericStateBlock)
class GenericStateBlockData(StateBlockData):
    """A modular, general purpose property package."""

    def build(self):
        """Callable method for Block construction."""
        super(GenericStateBlockData, self).build()

        # Add state variables and associated methods
        self.params.config.state_definition.define_state(self)

        # Create common components for each property package
        for p in self.params.phase_list:
            self.params.config.equation_of_state[p].common(self)

        # Add phase equilibrium constraints if necessary
        if (self.params.config.phase_equilibrium_formulation is not None and
                (not self.config.defined_state or self.always_flash)):
            self.params.config.phase_equilibrium_formulation.phase_equil(self)

    def components_in_phase(self, phase):
        """
        Generator method which yields components present in a given phase.

        Args:
            phase - phase for which to yield components

        Yields:
            components present in phase.
        """
        if self.params.config.phase_component_list is None:
            # All components in all phases
            for j in self.params.component_list:
                yield j
        else:
            # Return only components for indicated phase
            for j in self.params.config.phase_component_list[phase]:
                yield j

    # -------------------------------------------------------------------------
    # Bubble and Dew Points
    def _temperature_bubble(b):
        if b.params.config.temperature_bubble is None:
            raise GenericPropertyPackageError(b, "temperature_bubble")

        b.temperature_bubble = Var(
                doc="Bubble point temperature",
                bounds=(b.temperature.lb, b.temperature.ub))

        b._mole_frac_tbub = Var(
                b.params.component_list,
                initialize=1/len(b.params.component_list),
                bounds=(0, None),
                doc="Vapor mole fractions at bubble temperature")

        b.params.config.temperature_bubble(b)

    def _temperature_dew(b):
        if b.params.config.temperature_dew is None:
            raise GenericPropertyPackageError(b, "temperature_dew")

        b.temperature_dew = Var(
                doc="Dew point temperature",
                bounds=(b.temperature.lb, b.temperature.ub))

        b._mole_frac_tdew = Var(
                b.params.component_list,
                initialize=1/len(b.params.component_list),
                bounds=(0, None),
                doc="Liquid mole fractions at dew temperature")

        b.params.config.temperature_dew(b)

    def _pressure_bubble(b):
        if b.params.config.pressure_bubble is None:
            raise GenericPropertyPackageError(b, "pressure_bubble")

        b.pressure_bubble = Var(
                doc="Bubble point pressure",
                bounds=(b.pressure.lb, b.pressure.ub))

        b._mole_frac_pbub = Var(
                b.params.component_list,
                initialize=1/len(b.params.component_list),
                bounds=(0, None),
                doc="Vapor mole fractions at bubble pressure")

        b.params.config.pressure_bubble(b)

    def _pressure_dew(b):
        if b.params.config.pressure_dew is None:
            raise GenericPropertyPackageError(b, "pressure_dew")

        b.pressure_dew = Var(
                doc="Dew point pressure",
                bounds=(b.pressure.lb, b.pressure.ub))

        b._mole_frac_pdew = Var(
                b.params.component_list,
                initialize=1/len(b.params.component_list),
                bounds=(0, None),
                doc="Liquid mole fractions at dew pressure")

        b.params.config.pressure_dew(b)

    # -------------------------------------------------------------------------
    # Property Methods
    def _dens_mass(self):
        def rule_dens_mass(b):
            return sum(b.dens_mass_phase[p]*b.phase_frac[p]
                       for p in b.params.phase_list)
        self.dens_mass = Expression(
                doc="Mixture mass density",
                rule=rule_dens_mass)

    def _dens_mass_phase(self):
        def rule_dens_mass_phase(b, p):
            return b.params.config.equation_of_state[p].dens_mass_phase(b, p)
        self.dens_mass_phase = Expression(
                self.params.phase_list,
                doc="Mass density of each phase",
                rule=rule_dens_mass_phase)

    def _dens_mol(self):
        def rule_dens_mol(b):
            return sum(b.dens_mol_phase[p]*b.phase_frac[p]
                       for p in b.params.phase_list)
        self.dens_mol = Expression(
                doc="Mixture molar density",
                rule=rule_dens_mol)

    def _dens_mol_phase(self):
        def rule_dens_mol_phase(b, p):
            return b.params.config.equation_of_state[p].dens_mol_phase(b, p)
        self.dens_mol_phase = Expression(
                self.params.phase_list,
                doc="Molar density of each phase",
                rule=rule_dens_mol_phase)

    def _enth_mol(self):
        def rule_enth_mol(b):
            return sum(b.enth_mol_phase[p]*b.phase_frac[p]
                       for p in b.params.phase_list)
        self.enth_mol = Expression(rule=rule_enth_mol,
                                   doc="Mixture molar enthalpy")

    def _enth_mol_phase(self):
        def rule_enth_mol_phase(b, p):
            return b.params.config.equation_of_state[p].enth_mol_phase(b, p)
        self.enth_mol_phase = Expression(self.params.phase_list,
                                         rule=rule_enth_mol_phase)

    def _enth_mol_phase_comp(self):
        def rule_enth_mol_phase_comp(b, p, j):
            return b.params.config.equation_of_state[p] \
                    .enth_mol_phase_comp(b, p, j)
        self.enth_mol_phase_comp = Expression(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_enth_mol_phase_comp)

    def _entr_mol(self):
        def rule_entr_mol(b):
            return sum(b.entr_mol_phase[p]*b.phase_frac[p]
                       for p in b.params.phase_list)
        self.entr_mol = Expression(rule=rule_entr_mol,
                                   doc="Mixture molar entropy")

    def _entr_mol_phase(self):
        def rule_entr_mol_phase(b, p):
            return b.params.config.equation_of_state[p].entr_mol_phase(b, p)
        self.entr_mol_phase = Expression(self.params.phase_list,
                                         rule=rule_entr_mol_phase)

    def _entr_mol_phase_comp(self):
        def rule_entr_mol_phase_comp(b, p, j):
            return b.params.config.equation_of_state[p] \
                    .entr_mol_phase_comp(b, p, j)
        self.entr_mol_phase_comp = Expression(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_entr_mol_phase_comp)

    def _fug_phase_comp(self):
        def rule_fug_phase_comp(b, p, j):
            return b.params.config.equation_of_state[p] \
                .fug_phase_comp(b, p, j)
        self.fug_phase_comp = Expression(self.params.phase_list,
                                         self.params.component_list,
                                         rule=rule_fug_phase_comp)

    def _fug_coeff_phase_comp(self):
        def rule_fug_coeff_phase_comp(b, p, j):
            return b.params.config.equation_of_state[p] \
                .fug_coeff_phase_comp(b, p, j)
        self.fug_coeff_phase_comp = Expression(
                self.params.phase_list,
                self.params.component_list,
                rule=rule_fug_coeff_phase_comp)

    def _gibbs_mol(self):
        def rule_gibbs_mol(b):
            return sum(b.gibbs_mol_phase[p]*b.phase_frac[p]
                       for p in b.params.phase_list)
        self.gibbs_mol = Expression(rule=rule_gibbs_mol,
                                    doc="Mixture molar Gibbs energy")

    def _gibbs_mol_phase(self):
        def rule_gibbs_mol_phase(b, p):
            return b.params.config.equation_of_state[p].gibbs_mol_phase(b, p)
        self.gibbs_mol_phase = Expression(self.params.phase_list,
                                          rule=rule_gibbs_mol_phase)

    def _gibbs_mol_phase_comp(self):
        def rule_gibbs_mol_phase_comp(b, p, j):
            return b.params.config.equation_of_state[p] \
                    .gibbs_mol_phase_comp(b, p, j)
        self.gibbs_mol_phase_comp = Expression(
            self.params.phase_list,
            self.params.component_list,
            rule=rule_gibbs_mol_phase_comp)

    def _mw(self):
        self.mw = Expression(
                doc="Average molecular weight",
                expr=sum(self.phase_frac[p] *
                         sum(self.mole_frac_phase_comp[p, j] *
                             self.params.mw_comp[j]
                             for j in self.params.component_list)
                         for p in self.params.phase_list))

    def _mw_phase(self):
        def rule_mw_phase(b, p):
            return sum(b.mole_frac_phase_comp[p, j]*b.params.mw_comp[j]
                       for j in b.params.component_list)
        self.mw_phase = Expression(
                self.params.phase_list,
                doc="Average molecular weight of each phase",
                rule=rule_mw_phase)

    def _pressure_sat_comp(self):
        def rule_pressure_sat_comp(b, j):
            return get_method(b, "pressure_sat_comp")(b, j, b.temperature)
        self.pressure_sat_comp = Expression(
            self.params.component_list,
            rule=rule_pressure_sat_comp)
