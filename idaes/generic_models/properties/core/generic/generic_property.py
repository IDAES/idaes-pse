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
                           Param,
                           SolverFactory,
                           value,
                           Var)
from pyomo.common.config import ConfigValue

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock,
                        Component)
from idaes.core.phases import *
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


# TODO: Handle parameter data
# TODO: Set a default state definition
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


def get_method(self, config_arg, comp=None):
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
    if comp is None:
        source_block = self.params.config
    else:
        source_block = self.params.get_component(comp).config

    try:
        c_arg = getattr(source_block, config_arg)
    except AttributeError:
        raise AttributeError("{} Generic Property Package called for invalid "
                             "configuration option {}. Please contact the "
                             "developer of the property package."
                             .format(self.name, config_arg))

    if c_arg is None:
        raise GenericPropertyPackageError(self, config_arg)

    if isinstance(c_arg, types.ModuleType):
        return getattr(c_arg, config_arg)
    elif callable(c_arg):
        return c_arg
    else:
        raise ConfigurationError(
                "{} Generic Property Package received invalid value "
                "for argumnet {}. Value must be either a module or a "
                "method".format(self.name, config_arg))


def get_component_object(self, comp):
    """
    Utility method to get a component object from the property parameter block.
    This code is used frequently throughout the generic property pacakge
    libraries.

    Args:
        comp: name of the component object to be returned.

    Returns:
        Component: Component object with name comp.

    """
    return self.params.get_component(comp)


@declare_process_block_class("GenericParameterBlock")
class GenericParameterData(PhysicalParameterBlock):
    """
    General Property Parameter Block Class
    """
    CONFIG = PhysicalParameterBlock.CONFIG()

    # General options
    CONFIG.declare("components", ConfigValue(
        domain=dict,
        description="Dictionary of components in material",
        doc="""A dict of the components of interest in the mixture.
        Keys are component names and values are configuration arguments to
        be passed to Component on construction.
        """))
    CONFIG.declare("phases", ConfigValue(
        description="Dictionary of phases of interest",
        doc="""A dict of the phases of interest in the mixture.
        Keys are phases names and values are configuration arguments to
        be passed to Phase on construction.
        """))

    # TODO : Should we allow different state variables in each phase?
    CONFIG.declare("state_definition", ConfigValue(
        # default=FPhx,
        description="Choice of State Variables",
        doc="""Flag indicating the set of state variables to use for property
        package. Values should be a valid Python method which creates the
        required state variables."""))
    CONFIG.declare("state_bounds", ConfigValue(
        domain=dict,
        description="Bounds for state variables",
        doc="""A dict containing bounds to use for state variables."""))

    # Reference State
    CONFIG.declare("pressure_ref", ConfigValue(
        description="Pressure at reference state"))
    CONFIG.declare("temperature_ref", ConfigValue(
        description="Temperature at reference state"))

    # Phase equilibrium config arguments
    CONFIG.declare("phases_in_equilibrium", ConfigValue(
        default=None,
        domain=list,
        description="List of phase pairs which are in equilibrium",
        doc="""List of phase pairs for which equilibrium constraints should be
        constructed. Values should be a 2-tuples containing valid phase
        names. Default = None."""))
    CONFIG.declare("phase_equilibrium_formulation", ConfigValue(
        default=None,
        domain=dict,
        description="Formulation to use when calculating equilibrium state",
        doc="""Method to use for calculating phase equilibrium state and
        how to handle disappearing phases. Value should be a valid Python
        method or None. Default = None, indicating no phase equilibrium will
        occur."""))

    # Bubble and dew point methods
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

    def build(self):
        '''
        Callable method for Block construction.
        '''
        # Call super.build() to initialize Block
        super(GenericParameterData, self).build()

        # Call configure method to set construction arguments
        self.configure()

        # Set reference to StateBlock constructor class
        self.state_block_class = GenericStateBlock

        # Add Component objects
        if self.config.components is None:
            raise ConfigurationError(
                "{} was not provided with a components argument."
                .format(self.name))

        for c, d in self.config.components.items():
            self.add_component(c, Component(default=d))

        # Add Phase objects
        if self.config.phases is None:
            raise ConfigurationError(
                "{} was not provided with a phases argument."
                .format(self.name))

        for p, d in self.config.phases.items():
            tmp_dict = {}
            ptype = Phase
            for k, v in d.items():
                if k == "type":
                    ptype = v
                else:
                    tmp_dict[k] = v

            if ptype is Phase:
                _log.warning("{} phase {} was not assigned a type. "
                             "Using generic Phase object."
                             .format(self.name, p))
            self.add_component(str(p), ptype(default=tmp_dict))

        # Validate phase-component lists, and build _phase_component_set
        pc_set = []
        for p in self.phase_list:
            pc_list = self.get_phase(p).config.component_list
            if pc_list is None:
                # No phase-component list, assume all components in phase
                for j in self.component_list:
                    pc_set.append((p, j))
            else:
                # Validate that component names are valid and add to pc_set
                for j in pc_list:
                    if j not in self.component_list:
                        raise ConfigurationError(
                            "{} phase-component list for phase {} contained "
                            "component {} which is not in the master "
                            "component list".format(self.name, p, j))
                    pc_set.append((p, j))
        self._phase_component_set = Set(initialize=pc_set, ordered=True)

        # Validate state definition
        if self.config.state_definition is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "state_definition configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name))

        # Validate reference state and create Params
        if self.config.pressure_ref is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "pressure_ref configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name))
        else:
            self.pressure_ref = Param(
                initialize=self.config.pressure_ref,
                mutable=True)

        if self.config.temperature_ref is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "temperature_ref configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name))
        else:
            self.temperature_ref = Param(
                initialize=self.config.temperature_ref,
                mutable=True)

        # Validate equations of state
        for p in self.phase_list:
            if self.get_phase(p).config.equation_of_state is None:
                raise ConfigurationError(
                    "{} phase {} was not provided with an "
                    "equation_of_state configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name, p))

        # Validate and build phase equilibrium list
        if self.config.phases_in_equilibrium is not None:
            # List of interacting phases - assume all matching components
            # in phase pairs are in equilibrium
            pe_dict = {}
            pe_set = []
            counter = 1

            # Validate phase equilibrium formulation
            if self.config.phase_equilibrium_formulation is None:
                raise ConfigurationError(
                    "{} Generic Property Package provided with a "
                    "phases_in_equilibrium argument but no method was "
                    "specified for phase_equilibrium_formulation."
                    .format(self.name))
            pie_config = self.config.phase_equilibrium_formulation

            for pp in self.config.phases_in_equilibrium:
                if (pp not in pie_config.keys() and
                        (pp[1], pp[0]) not in pie_config.keys()):
                    raise ConfigurationError(
                        "{} Generic Property Package provided with a "
                        "phases_in_equilibrium argument but "
                        "phase_equilibrium_formulation was not specified "
                        "for all phase pairs."
                        .format(self.name))
                for j in self.component_list:
                    if ((pp[0], j) in self._phase_component_set
                            and (pp[1], j) in self._phase_component_set):
                        # Component j is in both phases, in equilibrium
                        pe_dict["PE"+str(counter)] = {j: (pp[0], pp[1])}
                        pe_set.append("PE"+str(counter))
                        counter += 1

                        # Validate that component has an equilibrium form
                        a = self.get_component(j).config.phase_equilibrium_form
                        if a is None:
                            raise ConfigurationError(
                                "{} Generic Property Package component {} is "
                                "in equilibrium but phase_equilibrium_form"
                                "was not specified."
                                .format(self.name, j))
                        elif (pp not in a.keys() and
                              (pp[1], pp[0]) not in a.keys()):
                            raise ConfigurationError(
                                "{} Generic Property Package component {} is "
                                "in equilibrium but phase_equilibrium_form "
                                "was not specified for all appropriate phase "
                                "pairs."
                                .format(self.name, j))

            # Construct phase_equilibrium_list and phase_equilibrium_idx
            self._pe_pairs = Set(initialize=self.config.phases_in_equilibrium,
                                 ordered=True)
            self.phase_equilibrium_list = pe_dict
            self.phase_equilibrium_idx = Set(initialize=pe_set,
                                             ordered=True)

        # Construct parameter components
        self.parameters()

        # For safety, fix all Vars in Component objects
        for c in self.component_list:
            cobj = self.get_component(c)
            for v in cobj.component_objects(Var):
                v.fix()

        self.config.state_definition.set_metadata(self)

    def configure(self):
        """
        Placeholder method to allow users to specify config arguments via a
        class. The user class should inherit from this one and implement a
        configure() method which sets the values of the desired config
        arguments.

        Args:
            None

        Returns:
            None
        """
        pass

    def parameters(self):
        """
        Placeholder method to allow users to specify parameters via a
        class. The user class should inherit from this one and implement a
        parameters() method which creates the reqruied components.

        Args:
            None

        Returns:
            None
        """
        pass

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        # TODO : Need to fix to have methods for things that may or may not be
        # created by state var methods
        obj.add_properties(
            {'flow_mol': {'method': None},
             'flow_mol_phase': {'method': None},
             'mole_frac_comp': {'method': None},
             'mole_frac_phase_comp': {'method': None},
             'phase_frac': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None},
             'dens_mass': {'method': '_dens_mass'},
             'dens_mass_phase': {'method': '_dens_mass_phase'},
             'dens_mol': {'method': '_dens_mol'},
             'dens_mol_phase': {'method': '_dens_mol_phase'},
             'enth_mol': {'method': '_enth_mol'},
             'enth_mol_phase': {'method': '_enth_mol_phase'},
             'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp'},
             'entr_mol': {'method': '_entr_mol'},
             'entr_mol_phase': {'method': '_entr_mol_phase'},
             'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp'},
             'fug_phase_comp': {'method': '_fug_phase_comp'},
             'fug_coeff_phase_comp': {'method': '_fug_coeff_phase_comp'},
             'gibbs_mol': {'method': '_gibbs_mol'},
             'gibbs_mol_phase': {'method': '_gibbs_mol_phase'},
             'gibbs_mol_phase_comp': {'method': '_gibbs_mol_phase_comp'},
             'mw': {'method': '_mw'},
             'mw_phase': {'method': '_mw_phase'},
             'pressure_bubble': {'method': '_pressure_bubble'},
             'pressure_dew': {'method': '_pressure_dew'},
             'pressure_sat_comp': {'method': '_pressure_sat_comp'},
             'temperature_bubble': {'method': '_temperature_bubble'},
             'temperature_dew': {'method': '_temperature_dew'}})

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
                Tbub0 = min(blk[k].params.get_component(j)
                            .temperature_crit.value
                            for j in blk[k].params.component_list) - 1

                err = 1
                counter = 0

                # Newton solver with step limiter to prevent overshoot
                # Tolerance only needs to be ~1e-1
                # Iteration limit of 30
                while err > 1e-1 and counter < 30:
                    f = value(sum(get_method(blk[k], "pressure_sat_comp", j)(
                                      blk[k],
                                      blk[k].params.get_component(j),
                                      Tbub0) *
                                  blk[k].mole_frac_comp[j]
                                  for j in blk[k].params.component_list) -
                              blk[k].pressure)
                    df = value(sum(
                           get_method(blk[k], "pressure_sat_comp", j)(
                                      blk[k],
                                      blk[k].params.get_component(j),
                                      Tbub0,
                                      dT=True)
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
                            get_method(blk[k], "pressure_sat_comp", j)(
                                       blk[k],
                                       blk[k].params.get_component(j),
                                       Tbub0))

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
                    Tdew0 = min(blk[k].params.get_component(j).temperature_crit
                                for j in blk[k].params.component_list) - 1

                err = 1
                counter = 0

                # Newton solver with step limiter to prevent overshoot
                # Tolerance only needs to be ~1e-1
                # Iteration limit of 30
                while err > 1e-1 and counter < 30:
                    f = value(blk[k].pressure *
                              sum(blk[k].mole_frac_comp[j] /
                                  get_method(blk[k], "pressure_sat_comp", j)(
                                       blk[k],
                                       blk[k].params.get_component(j),
                                       Tdew0)
                                  for j in blk[k].params.component_list) - 1)
                    df = -value(
                            blk[k].pressure *
                            sum(blk[k].mole_frac_comp[j] /
                                get_method(blk[k], "pressure_sat_comp", j)(
                                       blk[k],
                                       blk[k].params.get_component(j),
                                       Tdew0)**2 *
                                get_method(blk[k], "pressure_sat_comp", j)(
                                       blk[k],
                                       blk[k].params.get_component(j),
                                       Tdew0,
                                       dT=True)
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
                            get_method(blk[k], "pressure_sat_comp", j)(
                                       blk[k],
                                       blk[k].params.get_component(j),
                                       Tdew0))

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
                "Dew and bubble point initialization: {}."
                .format(idaeslog.condition(res))
            )
        # ---------------------------------------------------------------------
        # Calculate _teq
        for k in blk.keys():
            for pp in blk[k].params._pe_pairs:
                blk[k].params.config.phase_equilibrium_formulation[pp] \
                    .calculate_teq(blk[k], pp)

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
            for c in blk[k].component_objects(Constraint):
                # Activate common constraints
                if c.local_name in ("total_flow_balance",
                                    "component_flow_balances",
                                    "sum_mole_frac",
                                    "equilibrium_constraint"):
                    c.activate()
            for pp in blk[k].params._pe_pairs:
                # Activate formulation specific constraints
                blk[k].params.config.phase_equilibrium_formulation[pp] \
                    .phase_equil_initialization(blk[k], pp)

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
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
    def build(self):
        super(GenericStateBlockData, self).build()

        # Add state variables and associated methods
        self.params.config.state_definition.define_state(self)

        # Create common components for each property package
        for p in self.params.phase_list:
            self.params.get_phase(p).config.equation_of_state.common(self)

        # Add phase equilibrium constraints if necessary
        if (self.params.config.phases_in_equilibrium is not None and
                (not self.config.defined_state or self.always_flash)):

            # Add equilibrium temperature variable
            self._teq = Var(
                self.params._pe_pairs,
                initialize=value(self.temperature),
                doc='Temperature for calculating phase equilibrium')

            pe_form_config = self.params.config.phase_equilibrium_formulation
            for pp in self.params._pe_pairs:
                pe_form_config[pp].phase_equil(self, pp)

            def rule_equilibrium(b, phase1, phase2, j):
                config = b.params.get_component(j).config
                try:
                    e_mthd = config.phase_equilibrium_form[(phase1, phase2)]
                except KeyError:
                    e_mthd = config.phase_equilibrium_form[(phase2, phase1)]
                if e_mthd is None:
                    raise GenericPropertyPackageError(b,
                                                      "phase_equilibrium_form")
                return e_mthd(self, phase1, phase2, j)
            self.equilibrium_constraint = Constraint(
                self.params._pe_pairs,
                self.params.component_list,
                rule=rule_equilibrium)

    def components_in_phase(self, phase):
        """
        Generator method which yields components present in a given phase.

        Args:
            phase - phase for which to yield components

        Yields:
            components present in phase.
        """
        if self.params.get_phase(phase).config.component_list is None:
            # All components in all phases
            for j in self.params.component_list:
                yield j
        else:
            # Return only components for indicated phase
            for j in self.params.get_phase(phase).config.component_list:
                yield j

    # -------------------------------------------------------------------------
    # Bubble and Dew Points
    def _temperature_bubble(b):
        if b.params.config.temperature_bubble is None:
            raise GenericPropertyPackageError(b, "temperature_bubble")

        try:
            b.temperature_bubble = Var(
                    doc="Bubble point temperature",
                    bounds=(b.temperature.lb, b.temperature.ub))

            b._mole_frac_tbub = Var(
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Vapor mole fractions at bubble temperature")

            b.params.config.temperature_bubble(b)
        except AttributeError:
            b.del_component(b.temperature_bubble)
            b.del_component(b.mole_frac_tbub)
            raise

    def _temperature_dew(b):
        if b.params.config.temperature_dew is None:
            raise GenericPropertyPackageError(b, "temperature_dew")

        try:
            b.temperature_dew = Var(
                    doc="Dew point temperature",
                    bounds=(b.temperature.lb, b.temperature.ub))

            b._mole_frac_tdew = Var(
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Liquid mole fractions at dew temperature")

            b.params.config.temperature_dew(b)
        except AttributeError:
            b.del_component(b.temperature_dew)
            b.del_component(b.mole_frac_tdew)
            raise

    def _pressure_bubble(b):
        if b.params.config.pressure_bubble is None:
            raise GenericPropertyPackageError(b, "pressure_bubble")

        try:
            b.pressure_bubble = Var(
                    doc="Bubble point pressure",
                    bounds=(b.pressure.lb, b.pressure.ub))

            b._mole_frac_pbub = Var(
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Vapor mole fractions at bubble pressure")

            b.params.config.pressure_bubble(b)
        except AttributeError:
            b.del_component(b.pressure_bubble)
            b.del_component(b.mole_frac_pbub)
            raise

    def _pressure_dew(b):
        if b.params.config.pressure_dew is None:
            raise GenericPropertyPackageError(b, "pressure_dew")

        try:
            b.pressure_dew = Var(
                    doc="Dew point pressure",
                    bounds=(b.pressure.lb, b.pressure.ub))

            b._mole_frac_pdew = Var(
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Liquid mole fractions at dew pressure")

            b.params.config.pressure_dew(b)
        except AttributeError:
            b.del_component(b.pressure_dew)
            b.del_component(b.mole_frac_pdew)
            raise

    # -------------------------------------------------------------------------
    # Property Methods
    def _dens_mass(self):
        try:
            def rule_dens_mass(b):
                return sum(b.dens_mass_phase[p]*b.phase_frac[p]
                           for p in b.params.phase_list)
            self.dens_mass = Expression(
                    doc="Mixture mass density",
                    rule=rule_dens_mass)
        except AttributeError:
            self.del_component(self.dens_mass)
            raise

    def _dens_mass_phase(self):
        try:
            def rule_dens_mass_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.dens_mass_phase(b, p)
            self.dens_mass_phase = Expression(
                    self.params.phase_list,
                    doc="Mass density of each phase",
                    rule=rule_dens_mass_phase)
        except AttributeError:
            self.del_component(self.dens_mass_phass)
            raise

    def _dens_mol(self):
        try:
            def rule_dens_mol(b):
                return sum(b.dens_mol_phase[p]*b.phase_frac[p]
                           for p in b.params.phase_list)
            self.dens_mol = Expression(
                    doc="Mixture molar density",
                    rule=rule_dens_mol)
        except AttributeError:
            self.del_component(self.dens_mol)
            raise

    def _dens_mol_phase(self):
        try:
            def rule_dens_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.dens_mol_phase(b, p)
            self.dens_mol_phase = Expression(
                    self.params.phase_list,
                    doc="Molar density of each phase",
                    rule=rule_dens_mol_phase)
        except AttributeError:
            self.del_component(self.dens_mol_phase)
            raise

    def _enth_mol(self):
        try:
            def rule_enth_mol(b):
                return sum(b.enth_mol_phase[p]*b.phase_frac[p]
                           for p in b.params.phase_list)
            self.enth_mol = Expression(rule=rule_enth_mol,
                                       doc="Mixture molar enthalpy")
        except AttributeError:
            self.del_component(self.enth_mol)
            raise

    def _enth_mol_phase(self):
        try:
            def rule_enth_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.enth_mol_phase(b, p)
            self.enth_mol_phase = Expression(self.params.phase_list,
                                             rule=rule_enth_mol_phase)
        except AttributeError:
            self.del_component(self.enth_mol_phase)
            raise

    def _enth_mol_phase_comp(self):
        try:
            def rule_enth_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.enth_mol_phase_comp(b, p, j)
            self.enth_mol_phase_comp = Expression(
                self.params.phase_list,
                self.params.component_list,
                rule=rule_enth_mol_phase_comp)
        except AttributeError:
            self.del_component(self.enth_mol_phase_comp)
            raise

    def _entr_mol(self):
        try:
            def rule_entr_mol(b):
                return sum(b.entr_mol_phase[p]*b.phase_frac[p]
                           for p in b.params.phase_list)
            self.entr_mol = Expression(rule=rule_entr_mol,
                                       doc="Mixture molar entropy")
        except AttributeError:
            self.del_component(self.entr_mol)
            raise

    def _entr_mol_phase(self):
        try:
            def rule_entr_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.entr_mol_phase(b, p)
            self.entr_mol_phase = Expression(self.params.phase_list,
                                             rule=rule_entr_mol_phase)
        except AttributeError:
            self.del_component(self.entr_mol_phase)
            raise

    def _entr_mol_phase_comp(self):
        try:
            def rule_entr_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.entr_mol_phase_comp(b, p, j)
            self.entr_mol_phase_comp = Expression(
                self.params.phase_list,
                self.params.component_list,
                rule=rule_entr_mol_phase_comp)
        except AttributeError:
            self.del_component(self.entr_mol_phase_comp)
            raise

    def _fug_phase_comp(self):
        try:
            def rule_fug_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.fug_phase_comp(b, p, j)
            self.fug_phase_comp = Expression(self.params.phase_list,
                                             self.params.component_list,
                                             rule=rule_fug_phase_comp)
        except AttributeError:
            self.del_component(self.fug_phase_comp)
            raise

    def _fug_coeff_phase_comp(self):
        try:
            def rule_fug_coeff_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.fug_coeff_phase_comp(b, p, j)
            self.fug_coeff_phase_comp = Expression(
                    self.params.phase_list,
                    self.params.component_list,
                    rule=rule_fug_coeff_phase_comp)
        except AttributeError:
            self.del_component(self.fug_coeff_phase_comp)
            raise

    def _gibbs_mol(self):
        try:
            def rule_gibbs_mol(b):
                return sum(b.gibbs_mol_phase[p]*b.phase_frac[p]
                           for p in b.params.phase_list)
            self.gibbs_mol = Expression(rule=rule_gibbs_mol,
                                        doc="Mixture molar Gibbs energy")
        except AttributeError:
            self.del_component(self.gibbs_mol)
            raise

    def _gibbs_mol_phase(self):
        try:
            def rule_gibbs_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.gibbs_mol_phase(b, p)
            self.gibbs_mol_phase = Expression(self.params.phase_list,
                                              rule=rule_gibbs_mol_phase)
        except AttributeError:
            self.del_component(self.gibbs_mol_phase)
            raise

    def _gibbs_mol_phase_comp(self):
        try:
            def rule_gibbs_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.gibbs_mol_phase_comp(b, p, j)
            self.gibbs_mol_phase_comp = Expression(
                self.params.phase_list,
                self.params.component_list,
                rule=rule_gibbs_mol_phase_comp)
        except AttributeError:
            self.del_component(self.gibbs_mol_phase_comp)
            raise

    def _mw(self):
        try:
            self.mw = Expression(
                    doc="Average molecular weight",
                    expr=sum(self.phase_frac[p] *
                             sum(self.mole_frac_phase_comp[p, j] *
                                 self.params.get_component(j).mw_comp
                                 for j in self.params.component_list)
                             for p in self.params.phase_list))
        except AttributeError:
            self.del_component(self.mw)
            raise

    def _mw_phase(self):
        try:
            def rule_mw_phase(b, p):
                return sum(b.mole_frac_phase_comp[p, j] *
                           b.params.get_component(j).mw_comp
                           for j in b.params.component_list)
            self.mw_phase = Expression(
                    self.params.phase_list,
                    doc="Average molecular weight of each phase",
                    rule=rule_mw_phase)
        except AttributeError:
            self.del_component(self.mw_phase)
            raise

    def _pressure_sat_comp(self):
        try:
            def rule_pressure_sat_comp(b, j):
                cobj = b.params.get_Component(j)
                return get_method(b, "pressure_sat_comp")(
                    b, cobj, b.temperature)
            self.pressure_sat_comp = Expression(
                self.params.component_list,
                rule=rule_pressure_sat_comp)
        except AttributeError:
            self.del_component(self.pressure_sat_comp)
            raise
