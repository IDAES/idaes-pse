##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
                           Var,
                           units as pyunits)
from pyomo.common.config import ConfigValue, In
from pyomo.core.base.units_container import _PyomoUnit

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlockData,
                        StateBlock)
from idaes.core.components import Component, __all_components__
from idaes.core.phases import Phase, AqueousPhase, __all_phases__
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            solve_indexed_blocks)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_activated_constraints)
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError)
import idaes.logger as idaeslog

from idaes.generic_models.properties.core.generic.utility import (
    get_method, GenericPropertyPackageError)
from idaes.generic_models.properties.core.phase_equil.bubble_dew import \
    LogBubbleDew


# Set up logger
_log = idaeslog.getLogger(__name__)


def set_param_value(b, param, units):
    # We cannot use the standard method in core.util.misc as here the parameter
    # data is directly attached to the config block, rather than in a parameter
    # data entry.
    param_obj = getattr(b, param)
    config = getattr(b.config, param)
    if isinstance(config, tuple):
        param_obj.value = pyunits.convert_value(
            config[0], from_units=config[1], to_units=units)
    else:
        _log.debug("{} no units provided for parameter {} - assuming default "
                   "units".format(b.name, param))
        param_obj.value = config


# TODO: Set a default state definition
# TODO: Need way to dynamically determine units of measurement....
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
    CONFIG.declare("phase_equilibrium_state", ConfigValue(
        default=None,
        domain=dict,
        description="Formulation to use when calculating equilibrium state",
        doc="""Method to use for calculating phase equilibrium state and
        how to handle disappearing phases. Value should be a valid Python
        method or None. Default = None, indicating no phase equilibrium will
        occur."""))

    # Bubble and dew point methods
    CONFIG.declare("bubble_dew_method", ConfigValue(
        default=LogBubbleDew,
        description="Method to use to calculate bubble and dew points",
        doc="""Flag indicating what formulation to use for calculating bubble
        and dew points. Value should be a valid Python class."""))

    # General parameter data dict
    CONFIG.declare("parameter_data", ConfigValue(
        default={},
        domain=dict,
        description="Dict containing initialization data for parameters"))

    # Base units of measurement
    CONFIG.declare("base_units", ConfigValue(
        default={},
        domain=dict,
        description="Base units for property package",
        doc="Dict containing definition of base units of measurement to use "
        "with property package."))

    def build(self):
        '''
        Callable method for Block construction.
        '''
        # Call super.build() to initialize Block
        super(GenericParameterData, self).build()

        # Validate and set base units of measurement
        self.get_metadata().add_default_units(self.config.base_units)
        units_meta = self.get_metadata().default_units

        for key, unit in self.config.base_units.items():
            if key in ['time', 'length', 'mass', 'amount', 'temperature',
                       "current", "luminous intensity"]:
                if not isinstance(unit, _PyomoUnit):
                    raise ConfigurationError(
                        "{} recieved unexpected units for quantity {}: {}. "
                        "Units must be instances of a Pyomo unit object."
                        .format(self.name, key, unit))
            else:
                raise ConfigurationError(
                    "{} defined units for an unexpected quantity {}. "
                    "Generic property packages only support units for the 7 "
                    "base SI quantities.".format(self.name, key))

        # Check that main 5 base units are assigned
        for k in ['time', 'length', 'mass', 'amount', 'temperature']:
            if not isinstance(units_meta[k], _PyomoUnit):
                raise ConfigurationError(
                    "{} units for quantity {} were not assigned. "
                    "Please make sure to provide units for all base units "
                    "when configuring the property package."
                    .format(self.name, k))

        # Call configure method to set construction arguments
        self.configure()

        # Build core components
        self._state_block_class = GenericStateBlock

        # Add Phase objects
        if self.config.phases is None:
            raise ConfigurationError(
                "{} was not provided with a phases argument."
                .format(self.name))

        # Add a flag indicating whether this is an electrolyte system or not
        self._electrolyte = False

        for p, d in self.config.phases.items():
            ptype = d.pop("type", None)

            if ptype is None:
                _log.warning("{} phase {} was not assigned a type. "
                             "Using generic Phase object."
                             .format(self.name, p))
                ptype = Phase
            elif ptype not in __all_phases__:
                raise TypeError(
                    "{} phase {} was assigned unrecognised type {}."
                    .format(self.name, p, str(ptype)))
            elif ptype is AqueousPhase:
                # If there is an aqueous phase, set _electrolyte = True
                self._electrolyte = True

            self.add_component(str(p), ptype(default=d))

        # Check if we need to create electrolyte component lists
        if self._electrolyte:
            self.add_component(
                "anion_set",
                Set(ordered=True,
                    doc="Set of anions present in aqueous phase"))
            self.add_component(
                "cation_set",
                Set(ordered=True,
                    doc="Set of cations present in aqueous phase"))
            self.add_component(
                "solvent_set",
                Set(ordered=True,
                    doc="Set of solvent species in aqueous phase"))
            self.add_component(
                "solute_set",
                Set(ordered=True,
                    doc="Set of molecular solutes in aqueous phase"))
            self.add_component(
                "_apparent_set",
                Set(ordered=True,
                    doc="Set of apparent-only species in aqueous phase"))
            self.add_component(
                "_non_aqueous_set",
                Set(ordered=True,
                    doc="Set of components not present in aqueous phase"))

        # Add Component objects
        if self.config.components is None:
            raise ConfigurationError(
                "{} was not provided with a components argument."
                .format(self.name))

        for c, d in self.config.components.items():
            ctype = d.pop("type", None)
            d["_electrolyte"] = self._electrolyte

            if ctype is None:
                _log.warning("{} component {} was not assigned a type. "
                             "Using generic Component object."
                             .format(self.name, c))
                ctype = Component
            elif ctype not in __all_components__:
                raise TypeError(
                    "{} component {} was assigned unrecognised type {}."
                    .format(self.name, c, str(ctype)))

            self.add_component(c, ctype(default=d))

        # If this is an electrolyte system, we now need to build the actual
        # component lists
        if self._electrolyte:
            true_species = []
            apparent_species = []
            all_species = []

            for j in self.anion_set:
                true_species.append(j)
                all_species.append(j)
            for j in self.cation_set:
                true_species.append(j)
                all_species.append(j)
            for j in self.solvent_set:
                true_species.append(j)
                apparent_species.append(j)
                all_species.append(j)
            for j in self.solute_set:
                true_species.append(j)
                apparent_species.append(j)
                all_species.append(j)
            for j in self._apparent_set:
                apparent_species.append(j)
                all_species.append(j)
            for j in self._non_aqueous_set:
                true_species.append(j)
                apparent_species.append(j)
                all_species.append(j)

            self.add_component(
                "component_list",
                Set(initialize=all_species,
                    ordered=True,
                    doc="Master set of all components in mixture"))
            self.add_component(
                "true_species_set",
                Set(initialize=true_species,
                    ordered=True,
                    doc="Set of true components in mixture"))
            self.add_component(
                "apparent_species_set",
                Set(initialize=apparent_species,
                    ordered=True,
                    doc="Set of apparent components in mixture"))

        # Validate phase-component lists, and build _phase_component_set
        if not self._electrolyte:
            pc_set = []
            for p in self.phase_list:
                pobj = self.get_phase(p)
                pc_list = self.get_phase(p).config.component_list
                if pc_list is None:
                    # No phase-component list, look at components to determine
                    # which are valid in current phase
                    for j in self.component_list:
                        if self.get_component(j)._is_phase_valid(pobj):
                            # If compoennt says phase is valid, add to set
                            pc_set.append((p, j))
                else:
                    # Validate that component names are valid and add to pc_set
                    for j in pc_list:
                        if j not in self.component_list:
                            # Unrecognised component
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {} which is not in the "
                                "master component list"
                                .format(self.name, p, j))
                        # Check that phase is valid for component
                        if not self.get_component(j)._is_phase_valid(pobj):
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {}, however this "
                                "component is not valid for the given "
                                "PhaseType".format(self.name, p, j))
                        pc_set.append((p, j))
            self._phase_component_set = Set(initialize=pc_set, ordered=True)
        else:
            pc_set_appr = []
            pc_set_true = []
            for p in self.phase_list:
                pobj = self.get_phase(p)
                pc_list = self.get_phase(p).config.component_list
                if pc_list is None:
                    # No phase-component list, look at components to determine
                    # which are valid in current phase
                    for j in self.true_species_set:
                        if self.get_component(j)._is_phase_valid(pobj):
                            # If compoennt says phase is valid, add to set
                            pc_set_true.append((p, j))
                    for j in self.apparent_species_set:
                        if self.get_component(j)._is_phase_valid(pobj):
                            # If compoennt says phase is valid, add to set
                            pc_set_appr.append((p, j))
                else:
                    # Validate that component names are valid and add to pc_set
                    for j in pc_list:
                        if (j not in self.true_species_set and
                                j not in self.true_species_set):
                            # Unrecognised component
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {} which is not in the "
                                "master component list"
                                .format(self.name, p, j))
                        # Check that phase is valid for component
                        if not self.get_component(j)._is_phase_valid(pobj):
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {}, however this "
                                "component is not valid for the given "
                                "PhaseType".format(self.name, p, j))
                        if j not in self.true_species_set:
                            pc_set_true.append((p, j))
                        if j not in self.apparent_species_set:
                            pc_set_appr.append((p, j))
            self.true_phase_component_set = Set(initialize=pc_set_true,
                                                ordered=True)
            self.apparent_phase_component_set = Set(initialize=pc_set_appr,
                                                    ordered=True)

        # Validate and construct elemental composition objects as appropriate
        element_comp = {}
        for c in self.component_list:
            cobj = self.get_component(c)
            e_comp = cobj.config.elemental_composition

            if e_comp is None:
                # Do nothing
                continue
            else:
                for k, v in e_comp.items():
                    if not isinstance(v, int):
                        raise ConfigurationError(
                            "{} values in elemental_composition must be "
                            "integers (not floats): {}: {}."
                            .format(self.name, k, str(v)))
                element_comp[c] = e_comp

        if len(element_comp) == 0:
            # No elemental compositions defined, don't define components
            pass
        elif len(element_comp) != len(self.component_list):
            # Not all components defined elemental compositions
            raise ConfigurationError(
                "{} not all Components declared an elemental_composition "
                "argument. Either all Components must declare this, or none."
                .format(self.name))
        else:
            # Add elemental composition components
            self.element_list = Set(ordered=True)

            # Iterate through all componets and collect composing elements
            # Add these to element_list
            for ec in element_comp.values():
                for e in ec.keys():
                    if e not in self.element_list:
                        self.element_list.add(e)

            self.element_comp = {}
            for c in self.component_list:
                cobj = self.get_component(c)

                self.element_comp[c] = {}
                for e in self.element_list:

                    if e not in cobj.config.elemental_composition:
                        self.element_comp[c][e] = 0
                    else:
                        self.element_comp[c][e] = \
                            cobj.config.elemental_composition[e]

        # Validate state definition
        if self.config.state_definition is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "state_definition configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name))

        units = self.get_metadata().derived_units

        # Validate reference state and create Params
        if self.config.pressure_ref is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "pressure_ref configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name))
        else:
            self.pressure_ref = Param(
                mutable=True,
                units=units["pressure"])
            set_param_value(self, "pressure_ref", units["pressure"])

        if self.config.temperature_ref is None:
            raise ConfigurationError(
                    "{} Generic Property Package was not provided with a "
                    "temperature_ref configuration argument. Please fix "
                    "your property parameter definition to include this."
                    .format(self.name))
        else:
            self.temperature_ref = Param(
                mutable=True,
                units=units["temperature"])
            set_param_value(self, "temperature_ref", units["temperature"])

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
            if self.config.phase_equilibrium_state is None:
                raise ConfigurationError(
                    "{} Generic Property Package provided with a "
                    "phases_in_equilibrium argument but no method was "
                    "specified for phase_equilibrium_state."
                    .format(self.name))
            pie_config = self.config.phase_equilibrium_state

            for pp in self.config.phases_in_equilibrium:
                if (pp not in pie_config.keys() and
                        (pp[1], pp[0]) not in pie_config.keys()):
                    raise ConfigurationError(
                        "{} Generic Property Package provided with a "
                        "phases_in_equilibrium argument but "
                        "phase_equilibrium_state was not specified "
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
                                "in equilibrium but phase_equilibrium_form "
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

        # Construct parameters
        for c in self.component_list:
            cobj = self.get_component(c)
            for a, v in cobj.config.items():
                if isinstance(v, types.ModuleType):
                    c_arg = getattr(v, a)
                else:
                    c_arg = v
                if hasattr(c_arg, "build_parameters"):
                    try:
                        c_arg.build_parameters(cobj)
                    except KeyError:
                        raise ConfigurationError(
                            "{} values were not defined for parameter {} in "
                            "component {}. Please check the parameter_data "
                            "argument to ensure values are provided."
                            .format(self.name, a, c))

            # Validate and construct Henry parameters (indexed by phase)
            if cobj.config.henry_component is not None:
                for p, meth in cobj.config.henry_component.items():
                    # First validate that p is a phase
                    if p not in self.phase_list:
                        raise ConfigurationError(
                            "{} component {} was marked as a Henry's Law "
                            "component in phase {}, but this is not a valid "
                            "phase name.".format(self.name, c, p))
                    elif not self.get_phase(p).is_liquid_phase():
                        raise ConfigurationError(
                            "{} component {} was marked as a Henry's Law "
                            "component in phase {}, but this is not a Liquid "
                            "phase.".format(self.name, c, p))
                    else:
                        meth.build_parameters(cobj, p)

        for p in self.phase_list:
            pobj = self.get_phase(p)
            pobj.config.equation_of_state.build_parameters(pobj)

        # Call custom user parameter method
        self.parameters()

        # For safety, fix all Vars in Component objects
        for v in self.component_objects(Var, descend_into=True):
            for i in v:
                if v[i].value is None:
                    raise ConfigurationError(
                        "{} parameter {} was not assigned"
                        " a value. Please check your configuration "
                        "arguments.".format(self.name, v.local_name))
                v[i].fix()

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
        parameters() method which creates the required components.

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
             'flow_mol_phase_comp': {'method': None},
             'flow_mol_comp': {'method': None},
             'mole_frac_comp': {'method': None},
             'mole_frac_phase_comp': {'method': None},
             'phase_frac': {'method': None},
             'temperature': {'method': None},
             'pressure': {'method': None},
             'compress_fact_phase': {'method': '_compress_fact_phase'},
             'dens_mass': {'method': '_dens_mass'},
             'dens_mass_phase': {'method': '_dens_mass_phase'},
             'dens_mol': {'method': '_dens_mol'},
             'dens_mol_phase': {'method': '_dens_mol_phase'},
             'cp_mol': {'method': '_cp_mol'},
             'cp_mol_phase': {'method': '_cp_mol_phase'},
             'cp_mol_phase_comp': {'method': '_cp_mol_phase_comp'},
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


class _GenericStateBlock(StateBlock):
    """
    This Class contains methods which should be applied to Property Blocks as a
    whole, rather than individual elements of indexed Property Blocks.
    """

    def _return_component_list(self):
        # Overload the _return_component_list method to handle electrolyte
        # systems where we have two component lists to choose from
        params = self._get_parameter_block()
        if not params._electrolyte:
            return params.component_list
        elif self._block_data_config_default["species_basis"] == "true":
            return params.true_species_set
        elif self._block_data_config_default["species_basis"] == "apparent":
            return params.apparent_species_set
        else:
            raise BurntToast(
                "{} unrecognized value for configuration argument "
                "'species_basis'; this should never happen. Please contact "
                "the IDAES developers with this bug.".format(self.name))

    def _return_phase_component_set(self):
        # Overload the _return_phase_component_set method to handle electrolyte
        # systems where we have two component lists to choose from
        params = self._get_parameter_block()
        if not params._electrolyte:
            return params._phase_component_set
        elif self._block_data_config_default["species_basis"] == "true":
            return params.true_phase_component_set
        elif self._block_data_config_default["species_basis"] == "apparent":
            return params.apparent_phase_component_set
        else:
            raise BurntToast(
                "{} unrecognized value for configuration argument "
                "'species_basis'; this should never happen. Please contact "
                "the IDAES developers with this bug.".format(self.name))

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
                if blk[k].always_flash:
                    # If not always flash, DoF is probably less than zero
                    # We will handle this elsewhere
                    dof = degrees_of_freedom(blk[k])
                    if dof != 0:
                        raise BurntToast(
                            "Degrees of freedom were not zero [{}] "
                            "after trying to fix state variables. "
                            "Something broke in the generic property "
                            "package code - please inform the IDAES "
                            "developers.".format(dof))
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
            T_units = blk[k].params.get_metadata().default_units["temperature"]
            # Bubble temperature initialization
            if hasattr(blk[k], "_mole_frac_tbub"):
                for pp in blk[k].params._pe_pairs:
                    valid_comps = _valid_VL_component_list(blk[k], pp)

                    if valid_comps == []:
                        continue

                    # Use lowest component temperature_crit as starting point
                    # Starting high and moving down generally works better,
                    # as it under-predicts next step due to exponential form of
                    # Psat.
                    # Subtract 1 to avoid potential singularities at Tcrit
                    Tbub0 = min(blk[k].params.get_component(j)
                                .temperature_crit.value
                                for j in valid_comps) - 1

                    err = 1
                    counter = 0

                    # Newton solver with step limiter to prevent overshoot
                    # Tolerance only needs to be ~1e-1
                    # Iteration limit of 30
                    while err > 1e-1 and counter < 30:
                        f = value(sum(
                            get_method(blk[k], "pressure_sat_comp", j)(
                                    blk[k],
                                    blk[k].params.get_component(j),
                                    Tbub0*T_units) *
                            blk[k].mole_frac_comp[j]
                            for j in valid_comps) -
                            blk[k].pressure)
                        df = value(sum(
                               get_method(blk[k], "pressure_sat_comp", j)(
                                          blk[k],
                                          blk[k].params.get_component(j),
                                          Tbub0*T_units,
                                          dT=True)
                               for j in valid_comps))

                        # Limit temperature step to avoid excessive overshoot
                        if f/df < -50:
                            Tbub1 = Tbub0 + 50
                        elif f/df > 50:
                            Tbub1 = Tbub0 - 50
                        else:
                            Tbub1 = Tbub0 - f/df

                        err = abs(Tbub1 - Tbub0)
                        Tbub0 = Tbub1
                        counter += 1

                    blk[k].temperature_bubble[pp].value = Tbub0

                    for j in valid_comps:
                        blk[k]._mole_frac_tbub[pp, j].value = value(
                                blk[k].mole_frac_comp[j] *
                                get_method(blk[k], "pressure_sat_comp", j)(
                                           blk[k],
                                           blk[k].params.get_component(j),
                                           Tbub0*T_units) /
                                blk[k].pressure)

            # Dew temperature initialization
            if hasattr(blk[k], "_mole_frac_tdew"):
                for pp in blk[k].params._pe_pairs:
                    valid_comps = _valid_VL_component_list(blk[k], pp)

                    if valid_comps == []:
                        continue

                    if hasattr(blk[k], "_mole_frac_tbub"):
                        # If Tbub has been calculated above, use this as the
                        # starting point
                        Tdew0 = blk[k].temperature_bubble[pp].value
                    else:
                        # Otherwise, use lowest component critical temperature
                        # as starting point
                        # Subtract 1 to avoid potential singularities at Tcrit
                        Tdew0 = min(
                            blk[k].params.get_component(j).
                            temperature_crit.value
                            for j in valid_comps) - 1

                    err = 1
                    counter = 0

                    # Newton solver with step limiter to prevent overshoot
                    # Tolerance only needs to be ~1e-1
                    # Iteration limit of 30
                    while err > 1e-1 and counter < 30:
                        f = value(
                            blk[k].pressure *
                            sum(blk[k].mole_frac_comp[j] /
                                get_method(blk[k], "pressure_sat_comp", j)(
                                           blk[k],
                                           blk[k].params.get_component(j),
                                           Tdew0*T_units)
                                for j in valid_comps) - 1)
                        df = -value(
                                blk[k].pressure *
                                sum(blk[k].mole_frac_comp[j] /
                                    get_method(blk[k], "pressure_sat_comp", j)(
                                           blk[k],
                                           blk[k].params.get_component(j),
                                           Tdew0*T_units)**2 *
                                    get_method(blk[k], "pressure_sat_comp", j)(
                                           blk[k],
                                           blk[k].params.get_component(j),
                                           Tdew0*T_units,
                                           dT=True)
                                    for j in valid_comps))

                        # Limit temperature step to avoid excessive overshoot
                        if f/df < -50:
                            Tdew1 = Tdew0 + 50
                        elif f/df > 50:
                            Tdew1 = Tdew0 - 50
                        else:
                            Tdew1 = Tdew0 - f/df

                        err = abs(Tdew1 - Tdew0)
                        Tdew0 = Tdew1
                        counter += 1

                    blk[k].temperature_dew[pp].value = Tdew0

                    for j in valid_comps:
                        blk[k]._mole_frac_tdew[pp, j].value = value(
                                blk[k].mole_frac_comp[j]*blk[k].pressure /
                                get_method(blk[k], "pressure_sat_comp", j)(
                                           blk[k],
                                           blk[k].params.get_component(j),
                                           Tdew0*T_units))

            # Bubble pressure initialization
            if hasattr(blk[k], "_mole_frac_pbub"):
                for pp in blk[k].params._pe_pairs:
                    valid_comps = _valid_VL_component_list(blk[k], pp)

                    if valid_comps == []:
                        continue

                    blk[k].pressure_bubble[pp].value = value(
                            sum(blk[k].mole_frac_comp[j] *
                                blk[k].params.config.pressure_sat_comp
                                      .pressure_sat_comp(
                                              blk[k],
                                              j,
                                              blk[k].temperature*T_units)
                                for j in valid_comps))

                    for j in valid_comps:
                        blk[k]._mole_frac_pbub[pp, j].value = value(
                            blk[k].mole_frac_comp[j] *
                            blk[k].params.config.pressure_sat_comp
                                  .pressure_sat_comp(
                                          blk[k],
                                          j,
                                          blk[k].temperature*T_units) /
                            blk[k].pressure_bubble)

            # Dew pressure initialization
            if hasattr(blk[k], "_mole_frac_pdew"):
                for pp in blk[k].params._pe_pairs:
                    valid_comps = _valid_VL_component_list(blk[k], pp)

                    if valid_comps == []:
                        continue

                    blk[k].pressure_dew[pp].value = value(
                            sum(1/(blk[k].mole_frac_comp[j] /
                                   blk[k].params.config.pressure_sat_comp
                                   .pressure_sat_comp(
                                           blk[k],
                                           j,
                                           blk[k].temperature*T_units))
                                for j in valid_comps))

                    for j in valid_comps:
                        blk[k]._mole_frac_pdew[pp, j].value = value(
                            blk[k].mole_frac_comp[j]*blk[k].pressure_bubble /
                            blk[k].params.config.pressure_sat_comp
                            .pressure_sat_comp(blk[k],
                                               j,
                                               blk[k].temperature*T_units))

            # Solve bubble and dew point constraints
            for c in blk[k].component_objects(Constraint):
                # Deactivate all constraints not associated wtih bubble and dew
                # points
                if c.local_name not in ("eq_pressure_dew",
                                        "eq_pressure_bubble",
                                        "eq_temperature_dew",
                                        "eq_temperature_bubble",
                                        "eq_mole_frac_tbub",
                                        "eq_mole_frac_tdew",
                                        "eq_mole_frac_pbub",
                                        "eq_mole_frac_pdew",
                                        "mole_frac_comp_eq"):
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
        # Calculate _teq if required
        if (blk[k].params.config.phases_in_equilibrium is not None and
                (not blk[k].config.defined_state or blk[k].always_flash)):
            for k in blk.keys():
                for pp in blk[k].params._pe_pairs:
                    blk[k].params.config.phase_equilibrium_state[pp] \
                        .calculate_teq(blk[k], pp)

            init_log.info("Equilibrium temperature initialization completed.")

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        for k in blk.keys():
            blk[k].params.config.state_definition.state_initialization(blk[k])

            # If state block has phase equilibrium, use the average of all
            # _teq's as an initial guess for T
            if (blk[k].params.config.phases_in_equilibrium is not None and
                    isinstance(blk[k].temperature, Var) and
                    not blk[k].temperature.fixed):
                blk[k].temperature.value = value(
                    sum(blk[k]._teq[i] for i in blk[k].params._pe_pairs) /
                    len(blk[k].params._pe_pairs))

        if outlvl > 0:
            init_log.info("State variable initialization completed.")

        # ---------------------------------------------------------------------
        n_cons = 0
        skip = False
        for k in blk.keys():
            if (blk[k].params.config.phase_equilibrium_state is not None and
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
                    blk[k].params.config.phase_equilibrium_state[pp] \
                        .phase_equil_initialization(blk[k], pp)

            n_cons += number_activated_constraints(blk[k])
            if degrees_of_freedom(blk[k]) < 0:
                # Skip solve if DoF < 0 - this is probably due to a
                # phase-component flow state with flash
                skip = True

        if n_cons > 0 and not skip:
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

        n_cons = 0
        skip = False
        for k in blk:
            if degrees_of_freedom(blk[k]) < 0:
                # Skip solve if DoF < 0 - this is probably due to a
                # phase-component flow state with flash
                skip = True
            n_cons += number_activated_constraints(blk[k])
        if n_cons > 0 and not skip:
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
    CONFIG = StateBlockData.CONFIG()

    CONFIG.declare("species_basis", ConfigValue(
        default="true",
        domain=In(["true", "apparent"]),
        doc="True or apparent component basis",
        description="Argument idicating whetehr the true or apparent species "
        "set should be used for indexing components"))

    def build(self):
        super(GenericStateBlockData, self).build()

        # Add state variables and associated methods
        self.params.config.state_definition.define_state(self)

        # Add equilibrium temperature variable if required
        if (self.params.config.phases_in_equilibrium is not None and
                (not self.config.defined_state or self.always_flash)):

            t_units = self.params.get_metadata().default_units["temperature"]
            self._teq = Var(
                self.params._pe_pairs,
                initialize=value(self.temperature),
                doc='Temperature for calculating phase equilibrium',
                units=t_units)

        # Create common components for each property package
        for p in self.params.phase_list:
            pobj = self.params.get_phase(p)
            pobj.config.equation_of_state.common(self, pobj)

        # Add phase equilibrium constraints if necessary
        if (self.params.config.phases_in_equilibrium is not None and
                (not self.config.defined_state or self.always_flash)):

            pe_form_config = self.params.config.phase_equilibrium_state
            for pp in self.params._pe_pairs:
                pe_form_config[pp].phase_equil(self, pp)

            def rule_equilibrium(b, phase1, phase2, j):
                if ((phase1, j) not in b.params._phase_component_set or
                        (phase2, j) not in b.params._phase_component_set):
                    return Constraint.Skip
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
        for j in self.params.component_list:
            if (phase, j) in self.params._phase_component_set:
                yield j

    # -------------------------------------------------------------------------
    # Bubble and Dew Points
    def _temperature_bubble(b):
        if b.params.config.bubble_dew_method is None:
            raise GenericPropertyPackageError(b, "temperature_bubble")

        t_units = b.params.get_metadata().default_units["temperature"]
        try:
            b.temperature_bubble = Var(
                    b.params._pe_pairs,
                    doc="Bubble point temperature",
                    bounds=(b.temperature.lb, b.temperature.ub),
                    units=t_units)

            b._mole_frac_tbub = Var(
                    b.params._pe_pairs,
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Vapor mole fractions at bubble temperature",
                    units=None)

            b.params.config.bubble_dew_method.temperature_bubble(b)
        except AttributeError:
            b.del_component(b.temperature_bubble)
            b.del_component(b._mole_frac_tbub)
            raise

    def _temperature_dew(b):
        if b.params.config.bubble_dew_method is None:
            raise GenericPropertyPackageError(b, "temperature_dew")

        t_units = b.params.get_metadata().default_units["temperature"]
        try:
            b.temperature_dew = Var(
                    b.params._pe_pairs,
                    doc="Dew point temperature",
                    bounds=(b.temperature.lb, b.temperature.ub),
                    units=t_units)

            b._mole_frac_tdew = Var(
                    b.params._pe_pairs,
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Liquid mole fractions at dew temperature",
                    units=None)

            b.params.config.bubble_dew_method.temperature_dew(b)
        except AttributeError:
            b.del_component(b.temperature_dew)
            b.del_component(b._mole_frac_tdew)
            raise

    def _pressure_bubble(b):
        if b.params.config.bubble_dew_method is None:
            raise GenericPropertyPackageError(b, "pressure_bubble")

        units_meta = b.params.get_metadata().derived_units

        try:
            b.pressure_bubble = Var(
                    b.params._pe_pairs,
                    doc="Bubble point pressure",
                    bounds=(b.pressure.lb, b.pressure.ub),
                    units=units_meta["pressure"])

            b._mole_frac_pbub = Var(
                    b.params._pe_pairs,
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Vapor mole fractions at bubble pressure",
                    units=None)

            b.params.config.bubble_dew_method.pressure_bubble(b)
        except AttributeError:
            b.del_component(b.pressure_bubble)
            b.del_component(b._mole_frac_pbub)
            raise

    def _pressure_dew(b):
        if b.params.config.bubble_dew_method is None:
            raise GenericPropertyPackageError(b, "pressure_dew")

        units_meta = b.params.get_metadata().derived_units

        try:
            b.pressure_dew = Var(
                    b.params._pe_pairs,
                    doc="Dew point pressure",
                    bounds=(b.pressure.lb, b.pressure.ub),
                    units=units_meta["pressure"])

            b._mole_frac_pdew = Var(
                    b.params._pe_pairs,
                    b.params.component_list,
                    initialize=1/len(b.params.component_list),
                    bounds=(0, None),
                    doc="Liquid mole fractions at dew pressure",
                    units=None)

            b.params.config.bubble_dew_method.pressure_dew(b)
        except AttributeError:
            b.del_component(b.pressure_dew)
            b.del_component(b._mole_frac_pdew)
            raise

    # -------------------------------------------------------------------------
    # Property Methods
    def _compress_fact_phase(self):
        try:
            def rule_Z_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.compress_fact_phase(b, p)
            self.compress_fact_phase = Expression(
                    self.params.phase_list,
                    doc="Compressibility of each phase",
                    rule=rule_Z_phase)
        except AttributeError:
            self.del_component(self.compress_fact_phass)
            raise

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

    def _cp_mol(self):
        try:
            def rule_cp_mol(b):
                return sum(b.cp_mol_phase[p]*b.phase_frac[p]
                           for p in b.params.phase_list)
            self.cp_mol = Expression(rule=rule_cp_mol,
                                     doc="Mixture molar heat capacity")
        except AttributeError:
            self.del_component(self.cp_mol)
            raise

    def _cp_mol_phase(self):
        try:
            def rule_cp_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.cp_mol_phase(b, p)
            self.cp_mol_phase = Expression(self.params.phase_list,
                                           rule=rule_cp_mol_phase)
        except AttributeError:
            self.del_component(self.cp_mol_phase)
            raise

    def _cp_mol_phase_comp(self):
        try:
            def rule_cp_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.cp_mol_phase_comp(b, p, j)
            self.cp_mol_phase_comp = Expression(
                self.params._phase_component_set,
                rule=rule_cp_mol_phase_comp)
        except AttributeError:
            self.del_component(self.cp_mol_phase_comp)
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
                self.params._phase_component_set,
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
                self.params._phase_component_set,
                rule=rule_entr_mol_phase_comp)
        except AttributeError:
            self.del_component(self.entr_mol_phase_comp)
            raise

    def _fug_phase_comp(self):
        try:
            def rule_fug_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.fug_phase_comp(b, p, j)
            self.fug_phase_comp = Expression(self.params._phase_component_set,
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
                    self.params._phase_component_set,
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
                self.params._phase_component_set,
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


def _valid_VL_component_list(blk, pp):
    valid_comps = []
    # Only need to do this for V-L pairs, so check
    pparams = blk.params
    if ((pparams.get_phase(pp[0]).is_liquid_phase() and
         pparams.get_phase(pp[1]).is_vapor_phase()) or
        (pparams.get_phase(pp[0]).is_vapor_phase() and
         pparams.get_phase(pp[1]).is_liquid_phase())):

        for j in blk.params.component_list:
            if ((pp[0], j) in pparams._phase_component_set and
                    (pp[1], j) in pparams._phase_component_set):
                valid_comps.append(j)

    return valid_comps
