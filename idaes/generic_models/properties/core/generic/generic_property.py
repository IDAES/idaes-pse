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


# TODO: Set a default state definition
# TODO: Probably should set an initial value for state variables
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
    CONFIG.declare("phase_equilibrium_dict", ConfigValue(
        default=None,
        domain=dict,
        description="Phase equilibrium reactions to be modeled",
        doc="""Dict describing what phase equilibrium reactions should
        be included in the property calculations. Keys should be unique
        identifiers for each phase equilibrium reaction, and values should be
        3-tuples where the first value is a valid name from the component_list,
        and the following two values should be phases from the  phase_list."""
        ))
    CONFIG.declare("phase_equilibrium_formulation", ConfigValue(
        default=None,
        description="Formulation to use when calculating equilibrium state",
        doc="""Method to use for calculating phase equilibrium state and
        how to handle disappearing phases. Value should be a valid Python
        method or None. Default = None, indicating no phase equilibrium will
        occur."""))

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
        # Check that only one of phases_in_equilibrium and
        # phase_equilibrium_dict was provided
        pe_dict = {}
        pe_set = []
        if (self.config.phases_in_equilibrium is not None and
                self.config.phase_equilibrium_dict is not None):
            raise ConfigurationError(
                "{} Generic Property Package was provided with both a "
                "phases_in_equilibrium and a phase_equilibrium_dict "
                "argument. Users should provide only one of these."
                .format(self.name))
        else:
            pe_dict = {}
            pe_set = []

            if self.config.phases_in_equilibrium is not None:
                # List of interacting phases - assume all matching components
                # in phase pairs are in equilibrium
                counter = 1
                for pp in self.config.phases_in_equilibrium:
                    for j in self.component_list:
                        if ((pp[0], j) in self._phase_component_set
                                and (pp[1], j) in self._phase_component_set):
                            # Component j is in both phases, in equilibrium
                            pe_dict["PE"+str(counter)] = {j: (pp[0], pp[1])}
                            pe_set.append("PE"+str(counter))
                            counter += 1
            elif self.config.phase_equilibrium_dict is not None:
                # Provided a phase_equilibrium_dict - validate names
                for i, v in self.config.phase_equilibrium_dict.items():
                    if v[0] not in self.component_list:
                        raise ConfigurationError(
                            "{} phase_equilibrium_dict contained entry ({}) "
                            "with unrecognised component name {}."
                            .format(self.name, i, v[0]))
                    if v[1] not in self.phase_list:
                        raise ConfigurationError(
                            "{} phase_equilibrium_dict contained entry ({}) "
                            "with unrecognised phase name {}."
                            .format(self.name, i, v[1]))
                    if v[2] not in self.phase_list:
                        raise ConfigurationError(
                            "{} phase_equilibrium_dict contained entry ({}) "
                            "with unrecognised phase name {}."
                            .format(self.name, i, v[2]))
                pe_dict[i] = {v[0]: (v[1], v[2])}
                pe_set.append(i)

            # Construct phase_equilibrium_list and phase_equilibrium_idx
            self.phase_equilibrium_list = pe_dict
            self.phase_equilibrium_idx = Set(initialize=pe_set,
                                             ordered=True)

        # Validate phase equilibrium formulation if required
        if (self.config.phases_in_equilibrium is not None or
                self.config.phase_equilibrium_dict is not None):
            if self.config.phase_equilibrium_formulation is None:
                raise ConfigurationError(
                    "{} Generic Property Package provided with a "
                    "phases_in_equilibrium or phase_equilibrium_dict argument,"
                    " but no method was specified for "
                    "phase_equilibrium_formulation.".format(self.name))

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

    @classmethod
    def define_metadata(cls, obj):
        """Define properties supported and units."""
        # TODO : Need to fix to have methods for things that may or may not be
        # created by state var methods
    #     obj.add_properties(
    #         {'flow_mol': {'method': None, 'units': 'mol/s'},
    #          'mole_frac_comp': {'method': None, 'units': 'none'},
    #          'mole_frac_phase_comp': {'method': None, 'units': 'none'},
    #          'phase_frac': {'method': None, 'units': 'none'},
    #          'temperature': {'method': None, 'units': 'K'},
    #          'pressure': {'method': None, 'units': 'Pa'},
    #          'flow_mol_phase': {'method': None, 'units': 'mol/s'},
    #          'dens_mass': {'method': '_dens_mass', 'units': 'kg/m^3'},
    #          'dens_mass_phase': {'method': '_dens_mass_phase',
    #                              'units': 'kg/m^3'},
    #          'dens_mol': {'method': '_dens_mol', 'units': 'mol/m^3'},
    #          'dens_mol_phase': {'method': '_dens_mol_phase',
    #                             'units': 'mol/m^3'},
    #          'enth_mol': {'method': '_enth_mol', 'units': 'J/mol'},
    #          'enth_mol_phase': {'method': '_enth_mol_phase', 'units': 'J/mol'},
    #          'enth_mol_phase_comp': {'method': '_enth_mol_phase_comp',
    #                                  'units': 'J/mol'},
    #          'entr_mol': {'method': '_entr_mol', 'units': 'J/mol.K'},
    #          'entr_mol_phase': {'method': '_entr_mol_phase',
    #                             'units': 'J/mol.K'},
    #          'entr_mol_phase_comp': {'method': '_entr_mol_phase_comp',
    #                                  'units': 'J/mol.K'},
    #          'fug_phase_comp': {'method': '_fug_phase_comp', 'units': 'Pa'},
    #          'fug_coeff_phase_comp': {'method': '_fug_coeff_phase_comp',
    #                                   'units': '-'},
    #          'gibbs_mol': {'method': '_gibbs_mol', 'units': 'J/mol'},
    #          'gibbs_mol_phase': {'method': '_gibbs_mol_phase',
    #                              'units': 'J/mol'},
    #          'gibbs_mol_phase_comp': {'method': '_gibbs_mol_phase_comp',
    #                                   'units': 'J/mol'},
    #          'mw': {'method': '_mw', 'units': 'kg/mol'},
    #          'mw_phase': {'method': '_mw_phase', 'units': 'kg/mol'},
    #          'pressure_bubble': {'method': '_pressure_bubble', 'units': 'Pa'},
    #          'pressure_dew': {'method': '_pressure_dew', 'units': 'Pa'},
    #          'pressure_sat_comp': {'method': '_pressure_sat_comp',
    #                                'units': 'Pa'},
    #          'temperature_bubble': {'method': '_temperature_bubble',
    #                                 'units': 'K'},
    #          'temperature_dew': {'method': '_temperature_dew', 'units': 'K'}})

        obj.add_default_units({'time': 's',
                                'length': 'm',
                                'mass': 'g',
                                'amount': 'mol',
                                'temperature': 'K',
                                'energy': 'J',
                                'holdup': 'mol'})


class _GenericStateBlock(StateBlock):
    def initialize(blk, state_args={}, state_vars_fixed=False,
                   hold_state=False, outlvl=idaeslog.NOTSET,
                   solver='ipopt', optarg={'tol': 1e-8}):
        pass

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        pass

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
