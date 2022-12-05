#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Framework for generic property packages
"""
# Import Pyomo libraries
from pyomo.environ import (
    Block,
    check_optimal_termination,
    Constraint,
    exp,
    Expression,
    log,
    Set,
    Param,
    value,
    Var,
    units as pyunits,
    Reference,
)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    MaterialFlowBasis,
)
from idaes.core.base.components import Component, __all_components__
from idaes.core.base.phases import (
    Phase,
    AqueousPhase,
    LiquidPhase,
    VaporPhase,
    __all_phases__,
)
from idaes.core.util.initialization import (
    fix_state_vars,
    revert_state_vars,
    solve_indexed_blocks,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_activated_constraints,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyPackageError,
    PropertyNotSupportedError,
    InitializationError,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

from idaes.models.properties.modular_properties.base.generic_reaction import (
    equil_rxn_config,
)
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_phase_method,
    GenericPropertyPackageError,
    StateIndex,
)
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
    _valid_VL_component_list as bub_dew_VL_comp_list,
)
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType

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
            config[0], from_units=config[1], to_units=units
        )
    else:
        _log.debug(
            "{} no units provided for parameter {} - assuming default "
            "units".format(b.name, param)
        )
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
    CONFIG.declare(
        "components",
        ConfigValue(
            domain=dict,
            description="Dictionary of components in material",
            doc="""A dict of the components of interest in the mixture.
        Keys are component names and values are configuration arguments to
        be passed to Component on construction.
        """,
        ),
    )
    CONFIG.declare(
        "phases",
        ConfigValue(
            description="Dictionary of phases of interest",
            doc="""A dict of the phases of interest in the mixture.
        Keys are phases names and values are configuration arguments to
        be passed to Phase on construction.
        """,
        ),
    )

    # TODO : Should we allow different state variables in each phase?
    CONFIG.declare(
        "state_definition",
        ConfigValue(
            # default=FPhx,
            description="Choice of State Variables",
            doc="""Flag indicating the set of state variables to use for property
        package. Values should be a valid Python method which creates the
        required state variables.""",
        ),
    )
    CONFIG.declare(
        "state_bounds",
        ConfigValue(
            domain=dict,
            description="Bounds for state variables",
            doc="""A dict containing bounds to use for state variables.""",
        ),
    )
    CONFIG.declare(
        "state_components",
        ConfigValue(
            default=StateIndex.true,
            domain=In(StateIndex),
            doc="Index state variables by true or apparent components",
            description="Argument indicating whether the true or apparent species "
            "set should be used for indexing state variables. Must be "
            "StateIndex.true or StateIndex.apparent.",
        ),
    )

    # Reference State
    CONFIG.declare(
        "pressure_ref", ConfigValue(description="Pressure at reference state")
    )
    CONFIG.declare(
        "temperature_ref", ConfigValue(description="Temperature at reference state")
    )

    # Phase equilibrium config arguments
    CONFIG.declare(
        "phases_in_equilibrium",
        ConfigValue(
            default=None,
            domain=list,
            description="List of phase pairs which are in equilibrium",
            doc="""List of phase pairs for which equilibrium constraints should be
        constructed. Values should be a 2-tuples containing valid phase
        names. Default = None.""",
        ),
    )
    CONFIG.declare(
        "phase_equilibrium_state",
        ConfigValue(
            default=None,
            domain=dict,
            description="Formulation to use when calculating equilibrium state",
            doc="""Method to use for calculating phase equilibrium state and
        how to handle disappearing phases. Value should be a valid Python
        method or None. Default = None, indicating no phase equilibrium will
        occur.""",
        ),
    )

    # Bubble and dew point methods
    CONFIG.declare(
        "bubble_dew_method",
        ConfigValue(
            default=LogBubbleDew,
            description="Method to use to calculate bubble and dew points",
            doc="""Flag indicating what formulation to use for calculating bubble
        and dew points. Value should be a valid Python class.""",
        ),
    )

    # General parameter data dict
    CONFIG.declare(
        "parameter_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict containing initialization data for parameters",
        ),
    )

    # Base units of measurement
    CONFIG.declare(
        "base_units",
        ConfigValue(
            default={},
            domain=dict,
            description="Base units for property package",
            doc="Dict containing definition of base units of measurement to use "
            "with property package.",
        ),
    )

    # Property package options
    CONFIG.declare(
        "include_enthalpy_of_formation",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Include enthalpy of formation in property calculations",
            doc="Flag indiciating whether enthalpy of formation should be included"
            " when calculating specific enthalpies.",
        ),
    )

    # Config arguments for inherent reactions
    CONFIG.declare(
        "reaction_basis",
        ConfigValue(
            default=MaterialFlowBasis.molar,
            domain=In(MaterialFlowBasis),
            doc="Basis of reactions",
            description="Argument indicating basis of reaction terms. Should be "
            "an instance of a MaterialFlowBasis Enum",
        ),
    )

    CONFIG.declare(
        "inherent_reactions",
        ConfigBlock(implicit=True, implicit_domain=equil_rxn_config),
    )

    # User-defined default scaling factors
    CONFIG.declare(
        "default_scaling_factors",
        ConfigValue(
            domain=dict,
            description="User-defined default scaling factors",
            doc="Dict of user-defined properties and associated default "
            "scaling factors",
        ),
    )

    def build(self):
        """
        Callable method for Block construction.
        """
        # Call super.build() to initialize Block
        super(GenericParameterData, self).build()

        # Set base units of measurement
        self.get_metadata().add_default_units(self.config.base_units)

        # Call configure method to set construction arguments
        self.configure()

        # Build core components
        self._state_block_class = GenericStateBlock

        # Add Phase objects
        if self.config.phases is None:
            raise ConfigurationError(
                "{} was not provided with a phases argument.".format(self.name)
            )

        # Add a flag indicating whether this is an electrolyte system or not
        self._electrolyte = False

        for p, d in self.config.phases.items():
            # Create a copy of the phase config dict
            d = dict(d)

            ptype = d.pop("type", None)

            if ptype is None:
                _log.warning(
                    "{} phase {} was not assigned a type. "
                    "Using generic Phase object.".format(self.name, p)
                )
                ptype = Phase
            elif ptype not in __all_phases__:
                raise TypeError(
                    "{} phase {} was assigned unrecognised type {}.".format(
                        self.name, p, str(ptype)
                    )
                )
            elif ptype is AqueousPhase:
                # If there is an aqueous phase, set _electrolyte = True
                self._electrolyte = True
                # Check that specified property package supports electrolytes
                eos = d["equation_of_state"]
                if (
                    not hasattr(eos, "electrolyte_support")
                    or not eos.electrolyte_support
                ):
                    raise ConfigurationError(
                        "{} aqueous phase {} was set to use an equation of "
                        "state which does not support electrolytes: {}".format(
                            self.name, p, eos
                        )
                    )

            self.add_component(str(p), ptype(**d))

        # Check if we need to create electrolyte component lists
        if self._electrolyte:
            self.add_component(
                "anion_set",
                Set(ordered=True, doc="Set of anions present in aqueous phase"),
            )
            self.add_component(
                "cation_set",
                Set(ordered=True, doc="Set of cations present in aqueous phase"),
            )
            self.add_component(
                "solvent_set",
                Set(ordered=True, doc="Set of solvent species in aqueous phase"),
            )
            self.add_component(
                "solute_set",
                Set(ordered=True, doc="Set of molecular solutes in aqueous phase"),
            )
            self.add_component(
                "_apparent_set",
                Set(ordered=True, doc="Set of apparent-only species in aqueous phase"),
            )
            self.add_component(
                "_non_aqueous_set",
                Set(ordered=True, doc="Set of components not present in aqueous phase"),
            )

        # Add Component objects
        if self.config.components is None:
            raise ConfigurationError(
                "{} was not provided with a components argument.".format(self.name)
            )

        for c, d in self.config.components.items():
            # Create a copy of the component config dict
            d = dict(d)

            ctype = d.pop("type", None)
            d["_electrolyte"] = self._electrolyte

            if ctype is None:
                _log.warning(
                    "{} component {} was not assigned a type. "
                    "Using generic Component object.".format(self.name, c)
                )
                ctype = Component
            elif ctype not in __all_components__:
                raise TypeError(
                    "{} component {} was assigned unrecognised type {}.".format(
                        self.name, c, str(ctype)
                    )
                )

            self.add_component(c, ctype(**d))

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
                Set(
                    initialize=all_species,
                    ordered=True,
                    doc="Master set of all components in mixture",
                ),
            )
            self.add_component(
                "true_species_set",
                Set(
                    initialize=true_species,
                    ordered=True,
                    doc="Set of true components in mixture",
                ),
            )
            self.add_component(
                "apparent_species_set",
                Set(
                    initialize=apparent_species,
                    ordered=True,
                    doc="Set of apparent components in mixture",
                ),
            )
            self.add_component(
                "ion_set",
                Set(
                    initialize=self.anion_set | self.cation_set,
                    ordered=True,
                    doc="Master set of all ions in mixture",
                ),
            )

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
                            # If component says phase is valid, add to set
                            pc_set.append((p, j))
                else:
                    # Validate that component names are valid and add to pc_set
                    for j in pc_list:
                        if j not in self.component_list:
                            # Unrecognised component
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {} which is not in the "
                                "master component list".format(self.name, p, j)
                            )
                        # Check that phase is valid for component
                        if not self.get_component(j)._is_phase_valid(pobj):
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {}, however this "
                                "component is not valid for the given "
                                "PhaseType".format(self.name, p, j)
                            )
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
                            # If component says phase is valid, add to set
                            pc_set_true.append((p, j))
                    for j in self.apparent_species_set:
                        if self.get_component(j)._is_phase_valid(pobj):
                            # If component says phase is valid, add to set
                            pc_set_appr.append((p, j))

                            if not isinstance(pobj, AqueousPhase):
                                # Also need to add apparent species
                                if (p, j) not in pc_set_true:
                                    pc_set_true.append((p, j))

                else:
                    # Validate that component names are valid and add to pc_set
                    for j in pc_list:
                        if (
                            j not in self.true_species_set
                            and j not in self.apparent_species_set
                        ):
                            # Unrecognised component
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {} which is not in the "
                                "master component list".format(self.name, p, j)
                            )
                        # Check that phase is valid for component
                        if not self.get_component(j)._is_phase_valid(pobj):
                            raise ConfigurationError(
                                "{} phase-component list for phase {} "
                                "contained component {}, however this "
                                "component is not valid for the given "
                                "PhaseType".format(self.name, p, j)
                            )
                        if j in self.true_species_set:
                            pc_set_true.append((p, j))
                        if j in self.apparent_species_set:
                            pc_set_appr.append((p, j))
            self.true_phase_component_set = Set(initialize=pc_set_true, ordered=True)
            self.apparent_phase_component_set = Set(
                initialize=pc_set_appr, ordered=True
            )
            add_object_reference(
                self, "_phase_component_set", self.true_phase_component_set
            )

        # Check that each component appears phase-component set
        for j in self.component_list:
            count = 0
            for p in self.phase_list:
                if self._electrolyte:
                    if (p, j) in self.true_phase_component_set or (
                        p,
                        j,
                    ) in self.apparent_phase_component_set:
                        count += 1
                elif (p, j) in self._phase_component_set:
                    count += 1
            if count == 0:
                raise ConfigurationError(
                    "{} Component {} does not appear to be valid in any "
                    "phase. Please check the component lists defined for each "
                    "phase, and be sure you do not have generic Components "
                    "in single-phase aqueous systems.".format(self.name, j)
                )

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
                            "integers (not floats): {}: {}.".format(
                                self.name, k, str(v)
                            )
                        )
                element_comp[c] = e_comp

        if len(element_comp) == 0:
            # No elemental compositions defined, don't define components
            pass
        elif len(element_comp) != len(self.component_list):
            # Not all components defined elemental compositions
            raise ConfigurationError(
                "{} not all Components declared an elemental_composition "
                "argument. Either all Components must declare this, or none.".format(
                    self.name
                )
            )
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
                        self.element_comp[c][e] = cobj.config.elemental_composition[e]

        # Validate state definition
        if self.config.state_definition is None:
            raise ConfigurationError(
                "{} Generic Property Package was not provided with a "
                "state_definition configuration argument. Please fix "
                "your property parameter definition to include this.".format(self.name)
            )

        units = self.get_metadata().derived_units

        # Validate reference state and create Params
        if self.config.pressure_ref is None:
            raise ConfigurationError(
                "{} Generic Property Package was not provided with a "
                "pressure_ref configuration argument. Please fix "
                "your property parameter definition to include this.".format(self.name)
            )
        else:
            self.pressure_ref = Param(mutable=True, units=units.PRESSURE)
            set_param_value(self, "pressure_ref", units.PRESSURE)

        if self.config.temperature_ref is None:
            raise ConfigurationError(
                "{} Generic Property Package was not provided with a "
                "temperature_ref configuration argument. Please fix "
                "your property parameter definition to include this.".format(self.name)
            )
        else:
            self.temperature_ref = Param(mutable=True, units=units.TEMPERATURE)
            set_param_value(self, "temperature_ref", units.TEMPERATURE)

        # Validate equations of state
        for p in self.phase_list:
            if self.get_phase(p).config.equation_of_state is None:
                raise ConfigurationError(
                    "{} phase {} was not provided with an "
                    "equation_of_state configuration argument. Please fix "
                    "your property parameter definition to include this.".format(
                        self.name, p
                    )
                )

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
                    "specified for phase_equilibrium_state.".format(self.name)
                )
            pie_config = self.config.phase_equilibrium_state

            for pp in self.config.phases_in_equilibrium:
                if (
                    pp not in pie_config.keys()
                    and (pp[1], pp[0]) not in pie_config.keys()
                ):
                    raise ConfigurationError(
                        "{} Generic Property Package provided with a "
                        "phases_in_equilibrium argument but "
                        "phase_equilibrium_state was not specified "
                        "for all phase pairs.".format(self.name)
                    )

                for j in self.component_list:
                    if (pp[0], j) in self._phase_component_set and (
                        pp[1],
                        j,
                    ) in self._phase_component_set:
                        # Component j is in both phases, in equilibrium
                        pe_dict["PE" + str(counter)] = {j: (pp[0], pp[1])}
                        pe_set.append("PE" + str(counter))
                        counter += 1

                        # Validate that component has an equilibrium form
                        a = self.get_component(j).config.phase_equilibrium_form
                        if a is None:
                            raise ConfigurationError(
                                "{} Generic Property Package component {} is "
                                "in equilibrium but phase_equilibrium_form "
                                "was not specified.".format(self.name, j)
                            )
                        elif pp not in a.keys() and (pp[1], pp[0]) not in a.keys():
                            raise ConfigurationError(
                                "{} Generic Property Package component {} is "
                                "in equilibrium but phase_equilibrium_form "
                                "was not specified for all appropriate phase "
                                "pairs.".format(self.name, j)
                            )

            # Construct phase_equilibrium_list and phase_equilibrium_idx
            self._pe_pairs = Set(
                initialize=self.config.phases_in_equilibrium, ordered=True
            )
            self.phase_equilibrium_list = pe_dict
            self.phase_equilibrium_idx = Set(initialize=pe_set, ordered=True)

        # Construct parameters
        for c in self.component_list:
            cobj = self.get_component(c)
            for a, v in cobj.config.items():
                # Check to see if v has an attribute build_parameters
                if hasattr(v, "build_parameters"):
                    build_parameters = v.build_parameters
                else:
                    # If not, guess v is a class holding property subclasses
                    try:
                        build_parameters = getattr(v, a).build_parameters
                    except AttributeError:
                        # If all else fails, assume no build_parameters method
                        build_parameters = None

                # Call build_parameters if it exists
                if build_parameters is not None:
                    try:
                        build_parameters(cobj)
                    except KeyError:
                        raise ConfigurationError(
                            "{} values were not defined for parameter {} in "
                            "component {}. Please check the parameter_data "
                            "argument to ensure values are provided.".format(
                                self.name, a, c
                            )
                        )

            # Validate and construct Henry parameters (indexed by phase)
            if cobj.config.henry_component is not None:
                for p, d in cobj.config.henry_component.items():
                    # First validate that p is a phase
                    if p not in self.phase_list:
                        raise ConfigurationError(
                            "{} component {} was marked as a Henry's Law "
                            "component in phase {}, but this is not a valid "
                            "phase name.".format(self.name, c, p)
                        )
                    elif not self.get_phase(p).is_liquid_phase():
                        raise ConfigurationError(
                            "{} component {} was marked as a Henry's Law "
                            "component in phase {}, but this is not a Liquid "
                            "phase.".format(self.name, c, p)
                        )
                    else:
                        # Check that dict has necessary information
                        if "method" not in d.keys():
                            raise ConfigurationError(
                                f"{self.name} component {c} was marked as a "
                                f"Henry's Law component in phase {p}, but no "
                                f"method argument was provided."
                            )
                        elif "type" not in d.keys():
                            raise ConfigurationError(
                                f"{self.name} component {c} was marked as a "
                                f"Henry's Law component in phase {p}, but no "
                                f"type argument was provided."
                            )
                        elif not isinstance(d["type"], HenryType):
                            raise ConfigurationError(
                                f"{self.name} component {c} was marked as a "
                                f"Henry's Law component in phase {p}, but "
                                f"type argument was not an instance of "
                                f"HenryType."
                            )
                        elif (
                            self.config.phases_in_equilibrium is not None
                            and d["type"] != HenryType.Kpx
                        ):
                            raise PropertyNotSupportedError(
                                f"{self.name} currently only Kpx type Henry's "
                                f"constants are supported with full phase "
                                f"equilibrium. Support for other forms is a "
                                f"work-in-progress."
                            )
                        elif self._electrolyte and "basis" not in d.keys():
                            raise ConfigurationError(
                                f"{self.name} component {c} was marked as a "
                                f"Henry's Law component in phase {p}, but no "
                                f"basis argument was provided. Property "
                                f"packages using electrolytes must provide a "
                                f"basis in addition to a method and type."
                            )
                        try:
                            d["method"].build_parameters(cobj, p, d["type"])
                        except AttributeError:
                            # Method provided has no build_parameters method
                            # Assume it is not needed and continue
                            pass

        # Validate and other phase indexed props
        phase_indexed_props = ["diffus_phase_comp"]
        for prop in phase_indexed_props:
            for j in self.component_list:
                cobj = self.get_component(j)
                if cobj.config[prop] is not None:
                    for p, meth in cobj.config[prop].items():
                        # First validate that p is a phase
                        if p not in self.phase_list:
                            raise ConfigurationError(
                                f"{self.name} property {prop} definition "
                                f"contained unrecognised phase {p}."
                            )
                        else:
                            try:
                                meth.build_parameters(cobj, p)
                            except AttributeError:
                                # No build_parameters method
                                # Assume it is not needed and continue
                                continue

        for p in self.phase_list:
            pobj = self.get_phase(p)

            for a, v in pobj.config.items():
                # Check to see if v has an attribute build_parameters
                if hasattr(v, "build_parameters"):
                    build_parameters = v.build_parameters
                else:
                    # If not, guess v is a class holding property subclasses
                    try:
                        build_parameters = getattr(v, a).build_parameters
                    except AttributeError:
                        # If all else fails, assume no build_parameters method
                        build_parameters = None

                # Call build_parameters if it exists
                if build_parameters is not None:
                    try:
                        build_parameters(pobj)
                    except KeyError:
                        raise ConfigurationError(
                            "{} values were not defined for parameter {} in "
                            "phase {}. Please check the parameter_data "
                            "argument to ensure values are provided.".format(
                                self.name, a, p
                            )
                        )

        # Next, add inherent reactions if they exist
        if len(self.config.inherent_reactions) > 0:
            # Set has_inherent_reactions flag
            self._has_inherent_reactions = True

            # Construct inherent reaction index
            self.inherent_reaction_idx = Set(
                initialize=self.config.inherent_reactions.keys()
            )

            # Construct inherent reaction stoichiometry dict
            if self._electrolyte:
                pcset = self.true_phase_component_set
            else:
                pcset = self._phase_component_set

            self.inherent_reaction_stoichiometry = {}
            for r, rxn in self.config.inherent_reactions.items():
                for p, j in pcset:
                    self.inherent_reaction_stoichiometry[(r, p, j)] = 0

                if rxn.stoichiometry is None:
                    raise ConfigurationError(
                        "{} inherent reaction {} was not provided with a "
                        "stoichiometry configuration argument.".format(self.name, r)
                    )
                else:
                    for k, v in rxn.stoichiometry.items():
                        if k[0] not in self.phase_list:
                            raise ConfigurationError(
                                "{} stoichiometry for inherent reaction {} "
                                "included unrecognised phase {}.".format(
                                    self.name, r, k[0]
                                )
                            )
                        if k[1] not in self.component_list:
                            raise ConfigurationError(
                                "{} stoichiometry for inherent reaction {} "
                                "included unrecognised component {}.".format(
                                    self.name, r, k[1]
                                )
                            )
                        self.inherent_reaction_stoichiometry[(r, k[0], k[1])] = v

                # Check that a method was provided for the equilibrium form
                if rxn.equilibrium_form is None:
                    raise ConfigurationError(
                        "{} inherent reaction {} was not provided with a "
                        "equilibrium_form configuration argument.".format(self.name, r)
                    )

                # Construct blocks to contain parameters for each reaction
                self.add_component("reaction_" + str(r), Block())

                rblock = getattr(self, "reaction_" + r)
                r_config = self.config.inherent_reactions[r]

                order_init = {}
                for p, j in pcset:
                    if "reaction_order" in r_config.parameter_data:
                        try:
                            order_init[p, j] = r_config.parameter_data[
                                "reaction_order"
                            ][p, j]
                        except KeyError:
                            order_init[p, j] = 0
                    else:
                        # Assume elementary reaction and use stoichiometry
                        try:
                            # Here we use the stoic. coeff. directly
                            # However, solids should be excluded as they
                            # normally do not appear in the equilibrium
                            # relationship
                            pobj = self.get_phase(p)
                            if not pobj.is_solid_phase():
                                order_init[p, j] = r_config.stoichiometry[p, j]
                            else:
                                order_init[p, j] = 0
                        except KeyError:
                            order_init[p, j] = 0

                rblock.reaction_order = Var(
                    pcset, initialize=order_init, doc="Reaction order", units=None
                )

                for val in self.config.inherent_reactions[r].values():
                    try:
                        val.build_parameters(rblock, self.config.inherent_reactions[r])
                    except AttributeError:
                        pass

        # Call custom user parameter method
        self.parameters()

        # For safety, fix all Vars in Component objects
        for v in self.component_objects(Var, descend_into=True):
            for i in v:
                if v[i].value is None:
                    if i is None:  # Scalar Var
                        raise ConfigurationError(
                            "{} parameter {} was not assigned"
                            " a value. Please check your configuration "
                            "arguments.".format(self.name, v.local_name)
                        )
                    else:  # Indexed Var
                        raise ConfigurationError(
                            "{} parameter {}[{}] was not assigned"
                            " a value. Please check your configuration "
                            "arguments.".format(self.name, v.local_name, i)
                        )
                v[i].fix()

        self.config.state_definition.set_metadata(self)

        # Set default scaling factors
        # First, call set_default_scaling_factors method from state definiton
        try:
            self.config.state_definition.define_default_scaling_factors(self)
        except AttributeError:
            pass
        # Next, apply any user-defined scaling factors
        if self.config.default_scaling_factors is not None:
            self.default_scaling_factor.update(self.config.default_scaling_factors)
        # Finally, call populate_default_scaling_factors method to fill blanks
        iscale.populate_default_scaling_factors(self)

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
            {
                "flow_mol": {"method": "_flow_mol"},
                "flow_vol": {"method": "_flow_vol"},
                "flow_mass": {"method": "_flow_mass"},
                "flow_mass_phase": {"method": "_flow_mass_phase"},
                "flow_vol_phase": {"method": "_flow_vol_phase"},
                "flow_mol_phase": {"method": "_flow_mol_phase"},
                "flow_mass_comp": {"method": "_flow_mass_comp"},
                "flow_mol_comp": {"method": "_flow_mol_comp"},
                "flow_mass_phase_comp": {"method": "_flow_mass_phase_comp"},
                "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
                "mole_frac_comp": {"method": "_mole_frac_comp"},
                "mole_frac_phase_comp": {"method": None},
                "phase_frac": {"method": None},
                "temperature": {"method": None},
                "pressure": {"method": None},
                "act_phase_comp": {"method": "_act_phase_comp"},
                "act_phase_comp_true": {"method": "_act_phase_comp_true"},
                "act_phase_comp_apparent": {"method": "_act_phase_comp_apparent"},
                "act_coeff_phase_comp": {"method": "_act_coeff_phase_comp"},
                "act_coeff_phase_comp_true": {"method": "_act_coeff_phase_comp_true"},
                "act_coeff_phase_comp_apparent": {
                    "method": "_act_coeff_phase_comp_apparent"
                },
                "compress_fact_phase": {"method": "_compress_fact_phase"},
                "conc_mol_comp": {"method": "_conc_mol_comp"},
                "conc_mol_phase_comp": {"method": "_conc_mol_phase_comp"},
                "conc_mol_phase_comp_apparent": {
                    "method": "_conc_mol_phase_comp_apparent"
                },
                "conc_mol_phase_comp_true": {"method": "_conc_mol_phase_comp_true"},
                "cp_mol": {"method": "_cp_mol"},
                "cp_mol_phase": {"method": "_cp_mol_phase"},
                "cp_mol_phase_comp": {"method": "_cp_mol_phase_comp"},
                "cv_mol": {"method": "_cv_mol"},
                "cv_mol_phase": {"method": "_cv_mol_phase"},
                "cv_mol_phase_comp": {"method": "_cv_mol_phase_comp"},
                "diffus_phase_comp": {"method": "_diffus_phase_comp"},
                "diffus_phase_comp_apparent": {"method": "_diffus_phase_comp_apparent"},
                "diffus_phase_comp_true": {"method": "_diffus_phase_comp_true"},
                "heat_capacity_ratio_phase": {"method": "_heat_capacity_ratio_phase"},
                "dens_mass": {"method": "_dens_mass"},
                "dens_mass_phase": {"method": "_dens_mass_phase"},
                "dens_mol": {"method": "_dens_mol"},
                "dens_mol_phase": {"method": "_dens_mol_phase"},
                "energy_internal_mol": {"method": "_energy_internal_mol"},
                "energy_internal_mol_phase": {"method": "_energy_internal_mol_phase"},
                "energy_internal_mol_phase_comp": {
                    "method": "_energy_internal_mol_phase_comp"
                },
                "enth_mol": {"method": "_enth_mol"},
                "enth_mol_phase": {"method": "_enth_mol_phase"},
                "enth_mol_phase_comp": {"method": "_enth_mol_phase_comp"},
                "entr_mol": {"method": "_entr_mol"},
                "entr_mol_phase": {"method": "_entr_mol_phase"},
                "entr_mol_phase_comp": {"method": "_entr_mol_phase_comp"},
                "fug_phase_comp": {"method": "_fug_phase_comp"},
                "fug_coeff_phase_comp": {"method": "_fug_coeff_phase_comp"},
                "gibbs_mol": {"method": "_gibbs_mol"},
                "gibbs_mol_phase": {"method": "_gibbs_mol_phase"},
                "gibbs_mol_phase_comp": {"method": "_gibbs_mol_phase_comp"},
                "isentropic_speed_sound_phase": {
                    "method": "_isentropic_speed_sound_phase"
                },
                "isothermal_speed_sound_phase": {
                    "method": "_isothermal_speed_sound_phase"
                },
                "henry": {"method": "_henry"},
                "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
                "mass_frac_phase_comp_apparent": {
                    "method": "_mass_frac_phase_comp_apparent"
                },
                "mass_frac_phase_comp_true": {"method": "_mass_frac_phase_comp_true"},
                "molality_phase_comp": {"method": "_molality_phase_comp"},
                "molality_phase_comp_apparent": {
                    "method": "_molality_phase_comp_apparent"
                },
                "molality_phase_comp_true": {"method": "_molality_phase_comp_true"},
                "mw": {"method": "_mw"},
                "mw_comp": {"method": "_mw_comp"},
                "mw_phase": {"method": "_mw_phase"},
                "pressure_phase_comp": {"method": "_pressure_phase_comp"},
                "pressure_phase_comp_true": {"method": "_pressure_phase_comp_true"},
                "pressure_phase_comp_apparent": {
                    "method": "_pressure_phase_comp_apparent"
                },
                "pressure_bubble": {"method": "_pressure_bubble"},
                "pressure_dew": {"method": "_pressure_dew"},
                "pressure_osm_phase": {"method": "_pressure_osm_phase"},
                "pressure_sat_comp": {"method": "_pressure_sat_comp"},
                "surf_tens_phase": {"method": "_surf_tens_phase"},
                "temperature_bubble": {"method": "_temperature_bubble"},
                "temperature_dew": {"method": "_temperature_dew"},
                "therm_cond_phase": {"method": "_therm_cond_phase"},
                "visc_d_phase": {"method": "_visc_d_phase"},
                "vol_mol_phase": {"method": "_vol_mol_phase"},
                "vol_mol_phase_comp": {"method": "_vol_mol_phase_comp"},
                "dh_rxn": {"method": "_dh_rxn"},
                "log_act_phase_comp": {"method": "_log_act_phase_comp"},
                "log_act_phase_solvents": {"method": "_log_act_phase_solvents"},
                "log_act_phase_comp_true": {"method": "_log_act_phase_comp_true"},
                "log_act_phase_comp_apparent": {
                    "method": "_log_act_phase_comp_apparent"
                },
                "log_conc_mol_phase_comp": {"method": "_log_conc_mol_phase_comp"},
                "log_conc_mol_phase_comp_true": {
                    "method": "_log_conc_mol_phase_comp_true"
                },
                "log_mass_frac_phase_comp": {"method": "_log_mass_frac_phase_comp"},
                "log_mass_frac_phase_comp_apparent": {
                    "method": "_log_mass_frac_phase_comp_apparent"
                },
                "log_mass_frac_phase_comp_true": {
                    "method": "_log_mass_frac_phase_comp_true"
                },
                "log_molality_phase_comp": {"method": "_log_molality_phase_comp"},
                "log_molality_phase_comp_apparent": {
                    "method": "_log_molality_phase_comp_apparent"
                },
                "log_molality_phase_comp_true": {
                    "method": "_log_molality_phase_comp_true"
                },
                "log_mole_frac_comp": {"method": "_log_mole_frac_comp"},
                "log_mole_frac_tbub": {"method": "_log_mole_frac_tbub"},
                "log_mole_frac_tdew": {"method": "_log_mole_frac_tdew"},
                "log_mole_frac_pbub": {"method": "_log_mole_frac_pbub"},
                "log_mole_frac_pdew": {"method": "_log_mole_frac_pdew"},
                "log_mole_frac_phase_comp": {"method": "_log_mole_frac_phase_comp"},
                "log_mole_frac_phase_comp_apparent": {
                    "method": "_log_mole_frac_phase_comp_apparent"
                },
                "log_mole_frac_phase_comp_true": {
                    "method": "_log_mole_frac_phase_comp_true"
                },
                "log_pressure_phase_comp": {"method": "_log_pressure_phase_comp"},
                "log_pressure_phase_comp_apparent": {
                    "method": "_log_pressure_phase_comp_apparent"
                },
                "log_pressure_phase_comp_true": {
                    "method": "_log_pressure_phase_comp_true"
                },
                "log_k_eq": {"method": "_log_k_eq"},
            }
        )


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

        if params.config["state_components"] == StateIndex.true:
            return params.true_species_set
        elif params.config["state_components"] == StateIndex.apparent:
            return params.apparent_species_set
        else:
            raise BurntToast(
                "{} unrecognized value for configuration argument "
                "'state_components'; this should never happen. Please contact "
                "the IDAES developers with this bug.".format(self.name)
            )

    def _return_phase_component_set(self):
        # Overload the _return_phase_component_set method to handle electrolyte
        # systems where we have two component lists to choose from
        params = self._get_parameter_block()

        if not params._electrolyte:
            return params._phase_component_set

        if params.config["state_components"] == StateIndex.true:
            return params.true_phase_component_set
        elif params.config["state_components"] == StateIndex.apparent:
            return params.apparent_phase_component_set
        else:
            raise BurntToast(
                "{} unrecognized value for configuration argument "
                "'state_components'; this should never happen. Please contact "
                "the IDAES developers with this bug.".format(self.name)
            )

    def _include_inherent_reactions(self):
        params = self._get_parameter_block()

        if params.config["state_components"] == StateIndex.true:
            return params.has_inherent_reactions
        elif params.config["state_components"] == StateIndex.apparent:
            # If using apparent species basis, ignore inherent reactions
            return False
        else:
            raise BurntToast(
                "{} unrecognized value for configuration argument "
                "'state_components'; this should never happen. Please contact "
                "the IDAES developers with this bug.".format(self.name)
            )

    def initialize(
        blk,
        state_args=None,
        state_vars_fixed=False,
        hold_state=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Initialization routine for property package.
        Keyword Arguments:
            state_args : a dict of initial values for the state variables
                    defined by the property package.
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
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
                     initialization (default = None, use default solver)
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

        init_log.info("Starting initialization")

        res = None

        for k in blk.keys():
            # Deactivate the constraints specific for outlet block i.e.
            # when defined state is False
            if blk[k].config.defined_state is False:
                try:
                    blk[k].sum_mole_frac_out.deactivate()
                except AttributeError:
                    pass

                if hasattr(blk[k], "inherent_equilibrium_constraint") and (
                    not blk[k].params._electrolyte
                    or blk[k].params.config.state_components == StateIndex.true
                ):
                    blk[k].inherent_equilibrium_constraint.deactivate()

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
                            "developers.".format(dof)
                        )
        else:
            # When state vars are fixed, check that DoF is 0
            for k in blk.keys():
                if degrees_of_freedom(blk[k]) != 0:
                    raise Exception(
                        "State vars fixed but degrees of "
                        "freedom for state block is not zero "
                        "during initialization."
                    )

        # Create solver
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations
        for k in blk.keys():
            T_units = blk[k].params.get_metadata().default_units.TEMPERATURE
            # Bubble temperature initialization
            if hasattr(blk[k], "_mole_frac_tbub"):
                blk._init_Tbub(blk[k], T_units)

            # Dew temperature initialization
            if hasattr(blk[k], "_mole_frac_tdew"):
                blk._init_Tdew(blk[k], T_units)

            # Bubble pressure initialization
            if hasattr(blk[k], "_mole_frac_pbub"):
                blk._init_Pbub(blk[k], T_units)

            # Dew pressure initialization
            if hasattr(blk[k], "_mole_frac_pdew"):
                blk._init_Pdew(blk[k], T_units)

            # Solve bubble and dew point constraints
            for c in blk[k].component_objects(Constraint):
                # Deactivate all constraints not associated wtih bubble and dew
                # points
                if c.local_name not in (
                    "eq_pressure_dew",
                    "eq_pressure_bubble",
                    "eq_temperature_dew",
                    "eq_temperature_bubble",
                    "eq_mole_frac_tbub",
                    "eq_mole_frac_tdew",
                    "eq_mole_frac_pbub",
                    "eq_mole_frac_pdew",
                    "log_mole_frac_tbub_eqn",
                    "log_mole_frac_tdew_eqn",
                    "log_mole_frac_pbub_eqn",
                    "log_mole_frac_pdew_eqn",
                    "mole_frac_comp_eq",
                    "log_mole_frac_comp_eqn",
                ):
                    c.deactivate()

        # If StateBlock has active constraints (i.e. has bubble and/or dew
        # point calculations), solve the block to converge these
        n_cons = 0
        dof = 0
        for k in blk:
            n_cons += number_activated_constraints(blk[k])
            dof += degrees_of_freedom(blk[k])
        if n_cons > 0:
            if dof > 0:
                raise InitializationError(
                    f"{blk.name} Unexpected degrees of freedom during "
                    f"initialization at bubble and dew point step: {dof}."
                )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Dew and bubble point initialization: {}.".format(
                    idaeslog.condition(res)
                )
            )
        # ---------------------------------------------------------------------
        # Calculate _teq if required
        if blk[k].params.config.phases_in_equilibrium is not None and (
            not blk[k].config.defined_state or blk[k].always_flash
        ):
            for k in blk.keys():
                for pp in blk[k].params._pe_pairs:
                    blk[k].params.config.phase_equilibrium_state[pp].calculate_teq(
                        blk[k], pp
                    )

            init_log.info("Equilibrium temperature initialization completed.")

        # ---------------------------------------------------------------------
        # Initialize flow rates and compositions
        for k in blk.keys():

            blk[k].params.config.state_definition.state_initialization(blk[k])

            if blk[k].params._electrolyte:
                if blk[k].params.config.state_components == StateIndex.true:
                    # First calculate initial values for apparent species flows
                    for p, j in blk[k].params.apparent_phase_component_set:
                        calculate_variable_from_constraint(
                            blk[k].flow_mol_phase_comp_apparent[p, j],
                            blk[k].true_to_appr_species[p, j],
                        )
                    # Need to calculate all flows before doing mole fractions
                    for p, j in blk[k].params.apparent_phase_component_set:
                        sum_flow = sum(
                            blk[k].flow_mol_phase_comp_apparent[p, jj]
                            for jj in blk[k].params.apparent_species_set
                            if (p, jj) in blk[k].params.apparent_phase_component_set
                        )
                        if value(sum_flow) == 0:
                            x = 1
                        else:
                            x = value(
                                blk[k].flow_mol_phase_comp_apparent[p, j] / sum_flow
                            )
                        lb = blk[k].mole_frac_phase_comp_apparent[p, j].lb
                        if lb is not None and x <= lb:
                            blk[k].mole_frac_phase_comp_apparent[p, j].set_value(lb)
                        else:
                            blk[k].mole_frac_phase_comp_apparent[p, j].set_value(x)
                elif blk[k].params.config.state_components == StateIndex.apparent:
                    # First calculate initial values for true species flows
                    for p, j in blk[k].params.true_phase_component_set:
                        calculate_variable_from_constraint(
                            blk[k].flow_mol_phase_comp_true[p, j],
                            blk[k].appr_to_true_species[p, j],
                        )
                    # Need to calculate all flows before doing mole fractions
                    for p, j in blk[k].params.true_phase_component_set:
                        sum_flow = sum(
                            blk[k].flow_mol_phase_comp_true[p, jj]
                            for jj in blk[k].params.true_species_set
                            if (p, jj) in blk[k].params.true_phase_component_set
                        )
                        if value(sum_flow) == 0:
                            x = 1
                        else:
                            x = value(blk[k].flow_mol_phase_comp_true[p, j] / sum_flow)
                        lb = blk[k].mole_frac_phase_comp_true[p, j].lb
                        if lb is not None and x <= lb:
                            blk[k].mole_frac_phase_comp_true[p, j].set_value(lb)
                        else:
                            blk[k].mole_frac_phase_comp_true[p, j].set_value(x)

            # If state block has phase equilibrium, use the average of all
            # _teq's as an initial guess for T
            if (
                blk[k].params.config.phases_in_equilibrium is not None
                and isinstance(blk[k].temperature, Var)
                and not blk[k].temperature.fixed
            ):
                blk[k].temperature.value = value(
                    sum(blk[k]._teq[i] for i in blk[k].params._pe_pairs)
                    / len(blk[k].params._pe_pairs)
                )

        if outlvl > 0:
            init_log.info("State variable initialization completed.")

        # ---------------------------------------------------------------------
        n_cons = 0
        dof = 0
        skip = False
        Tfix = {}  # In enth based state defs, need to also fix T until later
        for k in blk.keys():
            if blk[k].params.config.phase_equilibrium_state is not None and (
                not blk[k].config.defined_state or blk[k].always_flash
            ):
                if not blk[k].temperature.fixed:
                    blk[k].temperature.fix()
                    Tfix[k] = True
                for c in blk[k].component_objects(Constraint):
                    # Activate common constraints
                    if c.local_name in (
                        "total_flow_balance",
                        "component_flow_balances",
                        "sum_mole_frac",
                        "phase_fraction_constraint",
                        "mole_frac_phase_comp_eq",
                        "mole_frac_comp_eq",
                    ):
                        c.activate()
                    if c.local_name == "log_mole_frac_phase_comp_eqn":
                        c.activate()
                        for p, j in blk[k].params._phase_component_set:
                            calculate_variable_from_constraint(
                                blk[k].log_mole_frac_phase_comp[p, j],
                                blk[k].log_mole_frac_phase_comp_eqn[p, j],
                            )
                    elif c.local_name == "equilibrium_constraint":
                        # For systems where the state variables fully define the
                        # phase equilibrium, we cannot activate the equilibrium
                        # constraint at this stage.
                        if "flow_mol_phase_comp" not in blk[k].define_state_vars():
                            c.activate()

                for pp in blk[k].params._pe_pairs:
                    # Activate formulation specific constraints
                    blk[k].params.config.phase_equilibrium_state[
                        pp
                    ].phase_equil_initialization(blk[k], pp)

            n_cons += number_activated_constraints(blk[k])
            dof += degrees_of_freedom(blk[k])
            if degrees_of_freedom(blk[k]) < 0:
                # Skip solve if DoF < 0 - this is probably due to a
                # phase-component flow state with flash
                skip = True

        if n_cons > 0 and not skip:
            if dof > 0:
                raise InitializationError(
                    f"{blk.name} Unexpected degrees of freedom during "
                    f"initialization at phase equilibrium step: {dof}."
                )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Phase equilibrium initialization: {}.".format(idaeslog.condition(res))
            )

        # ---------------------------------------------------------------------
        # Initialize other properties
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                # Activate all constraints except flagged do_not_initialize
                if c.local_name not in (
                    blk[k].params.config.state_definition.do_not_initialize
                ):
                    c.activate()
            if k in Tfix:
                blk[k].temperature.unfix()

            # Initialize log-form variables
            log_form_vars = [
                "act_phase_comp",
                "act_phase_comp_apparent",
                "act_phase_comp_true",
                "conc_mol_phase_comp",
                "conc_mol_phase_comp_apparent",
                "conc_mol_phase_comp_true",
                "mass_frac_phase_comp",
                "mass_frac_phase_comp_apparent",
                "mass_frac_phase_comp_true",
                "molality_phase_comp",
                "molality_phase_comp_apparent",
                "molality_phase_comp_true",
                "mole_frac_comp",  # Might have already been initialized
                "mole_frac_phase_comp",  # Might have already been initialized
                "mole_frac_phase_comp_apparent",
                "mole_frac_phase_comp_true",
                "pressure_phase_comp",
                "pressure_phase_comp_apparent",
                "pressure_phase_comp_true",
            ]

            for prop in log_form_vars:
                if blk[k].is_property_constructed("log_" + prop):
                    comp = getattr(blk[k], prop)
                    lcomp = getattr(blk[k], "log_" + prop)
                    for k2, v in lcomp.items():
                        c = value(comp[k2])
                        if c <= 0:
                            c = 1e-8
                        lc = log(c)
                        v.set_value(value(lc))

        n_cons = 0
        dof = 0
        skip = False
        for k in blk:
            if degrees_of_freedom(blk[k]) < 0:
                # Skip solve if DoF < 0 - this is probably due to a
                # phase-component flow state with flash
                skip = True
            n_cons += number_activated_constraints(blk[k])
            dof += degrees_of_freedom(blk[k])
        if n_cons > 0 and not skip:
            if dof > 0:
                raise InitializationError(
                    f"{blk.name} Unexpected degrees of freedom during "
                    f"initialization at property initialization step: {dof}."
                )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = solve_indexed_blocks(opt, [blk], tee=slc.tee)
            init_log.info(
                "Property initialization: {}.".format(idaeslog.condition(res))
            )

        # ---------------------------------------------------------------------
        # Return constraints to initial state
        for k in blk.keys():
            for c in blk[k].component_objects(Constraint):
                if c.local_name in (
                    blk[k].params.config.state_definition.do_not_initialize
                ):
                    c.activate()

        if res is not None and not check_optimal_termination(res):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )

        if state_vars_fixed is False:
            if hold_state is True:
                return flag_dict
            else:
                blk.release_state(flag_dict)

        init_log.info(
            "Property package initialization: {}.".format(idaeslog.condition(res))
        )

    def release_state(blk, flags, outlvl=idaeslog.NOTSET):
        """
        Method to relase state variables fixed during initialization.
        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state=True.
            outlvl : sets output level of initialization routine
        """
        revert_state_vars(blk, flags)
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="properties")
        init_log.info_high("State released.")

    def _init_Tbub(self, blk, T_units):
        for pp in blk.params._pe_pairs:
            raoult_comps, henry_comps = _valid_VL_component_list(blk, pp)

            if raoult_comps == []:
                continue
            if henry_comps != []:
                # Need to get liquid phase name
                # TODO This logic probably needs to be updated to break
                # when LLE is added
                if blk.params.get_phase(pp[0]).is_liquid_phase():
                    l_phase = pp[0]
                else:
                    l_phase = pp[1]

            # Use lowest component temperature_crit as starting point
            # Starting high and moving down generally works better,
            # as it under-predicts next step due to exponential form of
            # Psat.
            # Subtract 1 to avoid potential singularities at Tcrit
            Tbub0 = (
                min(
                    blk.params.get_component(j).temperature_crit.value
                    for j in raoult_comps
                )
                - 1
            )

            err = 1
            counter = 0

            # Newton solver with step limiter to prevent overshoot
            # Tolerance only needs to be ~1e-1
            # Iteration limit of 30
            while err > 1e-1 and counter < 30:
                f = value(
                    sum(
                        get_method(blk, "pressure_sat_comp", j)(
                            blk, blk.params.get_component(j), Tbub0 * T_units
                        )
                        * blk.mole_frac_comp[j]
                        for j in raoult_comps
                    )
                    + sum(
                        blk.mole_frac_comp[j]
                        * blk.params.get_component(j)
                        .config.henry_component[l_phase]["method"]
                        .return_expression(blk, l_phase, j, Tbub0 * T_units)
                        for j in henry_comps
                    )
                    - blk.pressure
                )
                df = value(
                    sum(
                        get_method(blk, "pressure_sat_comp", j)(
                            blk, blk.params.get_component(j), Tbub0 * T_units, dT=True
                        )
                        * blk.mole_frac_comp[j]
                        for j in raoult_comps
                    )
                    + sum(
                        blk.mole_frac_comp[j]
                        * blk.params.get_component(j)
                        .config.henry_component[l_phase]["method"]
                        .dT_expression(blk, l_phase, j, Tbub0 * T_units)
                        for j in henry_comps
                    )
                )

                # Limit temperature step to avoid excessive overshoot
                if f / df < -50:
                    Tbub1 = Tbub0 + 50
                elif f / df > 50:
                    Tbub1 = Tbub0 - 50
                else:
                    Tbub1 = Tbub0 - f / df

                err = abs(Tbub1 - Tbub0)
                Tbub0 = Tbub1
                counter += 1

            blk.temperature_bubble[pp].value = Tbub0

            for j in raoult_comps:
                blk._mole_frac_tbub[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * get_method(blk, "pressure_sat_comp", j)(
                        blk, blk.params.get_component(j), Tbub0 * T_units
                    )
                    / blk.pressure
                )
                if blk.is_property_constructed("log_mole_frac_tbub"):
                    blk.log_mole_frac_tbub[pp, j].value = value(
                        log(blk._mole_frac_tbub[pp, j])
                    )

            for j in henry_comps:
                blk._mole_frac_tbub[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * blk.params.get_component(j)
                    .config.henry_component[l_phase]["method"]
                    .return_expression(blk, l_phase, j, Tbub0 * T_units)
                    / blk.pressure
                )
                if blk.is_property_constructed("log_mole_frac_tbub"):
                    blk.log_mole_frac_tbub[pp, j].value = value(
                        log(blk._mole_frac_tbub[pp, j])
                    )

    def _init_Tdew(self, blk, T_units):
        for pp in blk.params._pe_pairs:
            raoult_comps, henry_comps = _valid_VL_component_list(blk, pp)

            if raoult_comps == []:
                continue
            if henry_comps != []:
                # Need to get liquid phase name
                if blk.params.get_phase(pp[0]).is_liquid_phase():
                    l_phase = pp[0]
                else:
                    l_phase = pp[1]

            if (
                hasattr(blk, "_mole_frac_tbub")
                and blk.temperature_bubble[pp].value is not None
            ):
                # If Tbub has been calculated above, use this as the
                # starting point
                Tdew0 = blk.temperature_bubble[pp].value
            else:
                # Otherwise, use lowest component critical temperature
                # as starting point
                # Subtract 1 to avoid potential singularities at Tcrit
                Tdew0 = (
                    min(
                        blk.params.get_component(j).temperature_crit.value
                        for j in raoult_comps
                    )
                    - 1
                )

            err = 1
            counter = 0

            # Newton solver with step limiter to prevent overshoot
            # Tolerance only needs to be ~1e-1
            # Iteration limit of 30
            while err > 1e-1 and counter < 30:
                f = value(
                    blk.pressure
                    * (
                        sum(
                            blk.mole_frac_comp[j]
                            / get_method(blk, "pressure_sat_comp", j)(
                                blk, blk.params.get_component(j), Tdew0 * T_units
                            )
                            for j in raoult_comps
                        )
                        + sum(
                            blk.mole_frac_comp[j]
                            / blk.params.get_component(j)
                            .config.henry_component[l_phase]["method"]
                            .return_expression(blk, l_phase, j, Tdew0 * T_units)
                            for j in henry_comps
                        )
                    )
                    - 1
                )
                df = -value(
                    blk.pressure
                    * (
                        sum(
                            blk.mole_frac_comp[j]
                            / get_method(blk, "pressure_sat_comp", j)(
                                blk, blk.params.get_component(j), Tdew0 * T_units
                            )
                            ** 2
                            * get_method(blk, "pressure_sat_comp", j)(
                                blk,
                                blk.params.get_component(j),
                                Tdew0 * T_units,
                                dT=True,
                            )
                            for j in raoult_comps
                        )
                        + sum(
                            blk.mole_frac_comp[j]
                            / blk.params.get_component(j)
                            .config.henry_component[l_phase]["method"]
                            .return_expression(blk, l_phase, j, Tdew0 * T_units)
                            ** 2
                            * blk.params.get_component(j)
                            .config.henry_component[l_phase]["method"]
                            .dT_expression(blk, l_phase, j, Tdew0 * T_units)
                            for j in henry_comps
                        )
                    )
                )

                # Limit temperature step to avoid excessive overshoot
                if f / df < -50:
                    Tdew1 = Tdew0 + 50
                elif f / df > 50:
                    Tdew1 = Tdew0 - 50
                else:
                    Tdew1 = Tdew0 - f / df

                err = abs(Tdew1 - Tdew0)
                Tdew0 = Tdew1
                counter += 1

            blk.temperature_dew[pp].value = Tdew0

            for j in raoult_comps:
                blk._mole_frac_tdew[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * blk.pressure
                    / get_method(blk, "pressure_sat_comp", j)(
                        blk, blk.params.get_component(j), Tdew0 * T_units
                    )
                )
                if blk.is_property_constructed("log_mole_frac_tdew"):
                    blk.log_mole_frac_tdew[pp, j].value = value(
                        log(blk._mole_frac_tdew[pp, j])
                    )
            for j in henry_comps:
                blk._mole_frac_tdew[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * blk.pressure
                    / blk.params.get_component(j)
                    .config.henry_component[l_phase]["method"]
                    .return_expression(blk, l_phase, j, Tdew0 * T_units)
                )
                if blk.is_property_constructed("log_mole_frac_tdew"):
                    blk.log_mole_frac_tdew[pp, j].value = value(
                        log(blk._mole_frac_tdew[pp, j])
                    )

    def _init_Pbub(self, blk, T_units):
        for pp in blk.params._pe_pairs:
            raoult_comps, henry_comps = _valid_VL_component_list(blk, pp)

            if raoult_comps == []:
                continue
            if henry_comps != []:
                # Need to get liquid phase name
                if blk.params.get_phase(pp[0]).is_liquid_phase():
                    l_phase = pp[0]
                else:
                    l_phase = pp[1]

            blk.pressure_bubble[pp].value = value(
                sum(
                    blk.mole_frac_comp[j] * blk.pressure_sat_comp[j]
                    for j in raoult_comps
                )
                + sum(
                    blk.mole_frac_comp[j] * blk.henry[l_phase, j] for j in henry_comps
                )
            )

            for j in raoult_comps:
                blk._mole_frac_pbub[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * blk.pressure_sat_comp[j]
                    / blk.pressure_bubble[pp]
                )
                if blk.is_property_constructed("log_mole_frac_pbub"):
                    blk.log_mole_frac_pbub[pp, j].value = value(
                        log(blk._mole_frac_pbub[pp, j])
                    )
            for j in henry_comps:
                blk._mole_frac_pbub[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * blk.henry[l_phase, j]
                    / blk.pressure_bubble[pp]
                )
                if blk.is_property_constructed("log_mole_frac_pbub"):
                    blk.log_mole_frac_pbub[pp, j].value = value(
                        log(blk._mole_frac_pbub[pp, j])
                    )

    def _init_Pdew(self, blk, T_units):
        for pp in blk.params._pe_pairs:
            raoult_comps, henry_comps = _valid_VL_component_list(blk, pp)

            if raoult_comps == []:
                continue
            if henry_comps != []:
                # Need to get liquid phase name
                if blk.params.get_phase(pp[0]).is_liquid_phase():
                    l_phase = pp[0]
                else:
                    l_phase = pp[1]

            blk.pressure_dew[pp].value = value(
                1
                / (
                    sum(
                        blk.mole_frac_comp[j] / blk.pressure_sat_comp[j]
                        for j in raoult_comps
                    )
                    + sum(
                        blk.mole_frac_comp[j] / blk.henry[l_phase, j]
                        for j in henry_comps
                    )
                )
            )

            for j in raoult_comps:
                blk._mole_frac_pdew[pp, j].value = value(
                    blk.mole_frac_comp[j]
                    * blk.pressure_dew[pp]
                    / blk.pressure_sat_comp[j]
                )
                if blk.is_property_constructed("log_mole_frac_pdew"):
                    blk.log_mole_frac_pdew[pp, j].value = value(
                        log(blk._mole_frac_pdew[pp, j])
                    )
            for j in henry_comps:
                blk._mole_frac_pdew[pp, j].value = value(
                    blk.mole_frac_comp[j] * blk.pressure_dew[pp] / blk.henry[l_phase, j]
                )
                if blk.is_property_constructed("log_mole_frac_pdew"):
                    blk.log_mole_frac_pdew[pp, j].value = value(
                        log(blk._mole_frac_pdew[pp, j])
                    )


@declare_process_block_class("GenericStateBlock", block_class=_GenericStateBlock)
class GenericStateBlockData(StateBlockData):
    CONFIG = StateBlockData.CONFIG()

    def build(self):
        super(GenericStateBlockData, self).build()

        # Add state variables and associated methods
        self.params.config.state_definition.define_state(self)

        # Add equilibrium temperature variable if required
        if self.params.config.phases_in_equilibrium is not None and (
            not self.config.defined_state or self.always_flash
        ):

            t_units = self.params.get_metadata().default_units.TEMPERATURE
            self._teq = Var(
                self.params._pe_pairs,
                initialize=value(self.temperature),
                doc="Temperature for calculating phase equilibrium",
                units=t_units,
            )

        # Create common components for each property package
        for p in self.phase_list:
            pobj = self.params.get_phase(p)
            pobj.config.equation_of_state.common(self, pobj)

        # Add phase equilibrium constraints if necessary
        if self.params.config.phases_in_equilibrium is not None and (
            not self.config.defined_state or self.always_flash
        ):

            pe_form_config = self.params.config.phase_equilibrium_state
            for pp in self.params._pe_pairs:
                pe_form_config[pp].phase_equil(self, pp)

            def rule_equilibrium(b, phase1, phase2, j):
                if (phase1, j) not in b.phase_component_set or (
                    phase2,
                    j,
                ) not in b.phase_component_set:
                    return Constraint.Skip
                config = b.params.get_component(j).config
                try:
                    e_mthd = config.phase_equilibrium_form[
                        (phase1, phase2)
                    ].return_expression
                except KeyError:
                    e_mthd = config.phase_equilibrium_form[
                        (phase2, phase1)
                    ].return_expression
                if e_mthd is None:
                    raise GenericPropertyPackageError(b, "phase_equilibrium_form")
                return e_mthd(self, phase1, phase2, j)

            self.equilibrium_constraint = Constraint(
                self.params._pe_pairs, self.component_list, rule=rule_equilibrium
            )

        # Add inherent reaction constraints if necessary
        if self.params.has_inherent_reactions and (
            not self.config.defined_state
            or (
                self.params._electrolyte
                and self.params.config.state_components == StateIndex.apparent
            )
        ):

            def equil_rule(b, r):
                rblock = getattr(b.params, "reaction_" + r)

                carg = b.params.config.inherent_reactions[r]

                return carg["equilibrium_form"].return_expression(
                    b, rblock, r, b.temperature
                )

            def keq_rule(b, r):
                rblock = getattr(b.params, "reaction_" + r)

                carg = b.params.config.inherent_reactions[r]

                return carg["equilibrium_constant"].return_expression(
                    b, rblock, r, b.temperature
                )

            self.k_eq = Expression(
                self.params.inherent_reaction_idx,
                doc="Equilibrium constant for inherent reactions",
                rule=keq_rule,
            )

            self.inherent_equilibrium_constraint = Constraint(
                self.params.inherent_reaction_idx,
                doc="Inherent reaction equilibrium constraint",
                rule=equil_rule,
            )

    def calculate_scaling_factors(self):
        # Get default scale factors and do calculations from base classes
        super().calculate_scaling_factors()

        # Scale state variables and associated constraints
        self.params.config.state_definition.calculate_scaling_factors(self)

        sf_T = iscale.get_scaling_factor(self.temperature, default=1, warning=True)
        sf_P = iscale.get_scaling_factor(self.pressure, default=1, warning=True)
        sf_mf = {}
        for i, v in self.mole_frac_phase_comp.items():
            sf_mf[i] = iscale.get_scaling_factor(v, default=1e3, warning=True)

        # Add scaling for components in build method
        # Phase equilibrium temperature
        if hasattr(self, "_teq"):
            for v in self._teq.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, sf_T)

        # Other EoS variables and constraints
        for p in self.phase_list:
            pobj = self.params.get_phase(p)
            pobj.config.equation_of_state.calculate_scaling_factors(self, pobj)

        # Flow and density terms
        if self.is_property_constructed("_enthalpy_flow_term"):
            for k, v in self._enthalpy_flow_term.items():
                if iscale.get_scaling_factor(v) is None:
                    sf_flow_phase = iscale.get_scaling_factor(
                        self.flow_mol_phase[k],
                        default=1,
                        warning=True,
                        hint="for _enthalpy_flow_term",
                    )
                    sf_h = iscale.get_scaling_factor(
                        self.enth_mol_phase[k],
                        default=1,
                        warning=True,
                        hint="for _enthalpy_flow_term",
                    )
                    iscale.set_scaling_factor(v, sf_flow_phase * sf_h)

        if self.is_property_constructed("_material_density_term"):
            for (p, j), v in self._material_density_term.items():
                if iscale.get_scaling_factor(v) is None:
                    sf_rho = iscale.get_scaling_factor(
                        self.dens_mol_phase[p], default=1, warning=True
                    )
                    sf_x = iscale.get_scaling_factor(
                        self.mole_frac_phase_comp[p, j], default=1, warning=True
                    )
                    iscale.set_scaling_factor(v, sf_rho * sf_x)

        if self.is_property_constructed("_energy_density_term"):
            for k, v in self._enthalpy_flow_term.items():
                if iscale.get_scaling_factor(v) is None:
                    sf_rho = iscale.get_scaling_factor(
                        self.dens_mol_phase[k], default=1, warning=True
                    )
                    sf_u = iscale.get_scaling_factor(
                        self.energy_internal_mol_phase[k], default=1, warning=True
                    )
                    iscale.set_scaling_factor(v, sf_rho * sf_u)

        # Phase equilibrium constraint
        if hasattr(self, "equilibrium_constraint"):
            pe_form_config = self.params.config.phase_equilibrium_state
            for pp in self.params._pe_pairs:
                pe_form_config[pp].calculate_scaling_factors(self, pp)

            for k in self.equilibrium_constraint:
                try:
                    sf_fug = (
                        self.params.get_component(k[2])
                        .config.phase_equilibrium_form[(k[0], k[1])]
                        .calculate_scaling_factors(self, k[0], k[1], k[2])
                    )

                    iscale.constraint_scaling_transform(
                        self.equilibrium_constraint[k], sf_fug, overwrite=False
                    )
                except KeyError:  # component not in phase
                    pass

        # Inherent reactions
        if hasattr(self, "k_eq"):
            for r in self.params.inherent_reaction_idx:
                rblock = getattr(self.params, "reaction_" + r)
                carg = self.params.config.inherent_reactions[r]

                sf_keq = iscale.get_scaling_factor(self.k_eq[r])
                if sf_keq is None:
                    sf_keq = carg["equilibrium_constant"].calculate_scaling_factors(
                        self, rblock
                    )
                    iscale.set_scaling_factor(self.k_eq[r], sf_keq)

                sf_const = carg["equilibrium_form"].calculate_scaling_factors(
                    self, sf_keq
                )

                iscale.constraint_scaling_transform(
                    self.inherent_equilibrium_constraint[r], sf_const, overwrite=False
                )

        # Add scaling for additional Vars and Constraints
        # Bubble and dew points
        def bubble_dew_scaling(b, pt_var):
            # Ditch the m.fs.unit.control_volume...
            short_name = pt_var.name.split(".")[-1]

            if short_name.startswith("temperature"):
                abbrv = "t"
                sf_pt = sf_T
            elif short_name.startswith("pressure"):
                abbrv = "p"
                sf_pt = sf_P
            else:
                _raise_dev_burnt_toast()

            if short_name.endswith("bubble"):
                phase = VaporPhase
                abbrv += "bub"
            elif short_name.endswith("dew"):
                phase = LiquidPhase
                abbrv += "dew"

            x_var = getattr(b, "_mole_frac_" + abbrv)

            if b.is_property_constructed("log_mole_frac_" + abbrv):
                log_eq = getattr(b, "log_mole_frac_" + abbrv + "_eqn")
            else:
                log_eq = None

            # Directly scale the bubble/dew temperature/pressure variable
            for v in pt_var.values():
                if iscale.get_scaling_factor(v) is None:
                    iscale.set_scaling_factor(v, sf_pt)

            # Scale mole fractions for bubble/dew calcs
            for i, v in x_var.items():
                if iscale.get_scaling_factor(v) is None:
                    if b.params.config.phases[i[0]]["type"] is phase:
                        p = i[0]
                    elif b.params.config.phases[i[1]]["type"] is phase:
                        p = i[1]
                    else:
                        # We create bubble/dew variables for all phase
                        # equilibrium pairs, regardless of whether it makes
                        # sense. If the pair doesn't make sense, the constraint
                        # is not created and the scaling factor is arbitrary
                        p = i[0]
                    try:
                        iscale.set_scaling_factor(v, sf_mf[p, i[2]])
                        if log_eq is not None and (
                            iscale.get_scaling_factor(log_eq[i]) is None
                        ):
                            iscale.constraint_scaling_transform(
                                log_eq[i], sf_mf[p, i[2]], overwrite=False
                            )
                    except KeyError:
                        # component i[2] is not in the new phase, so this
                        # variable is likely unused and scale doesn't matter
                        iscale.set_scaling_factor(v, 1)

            scaling_method = getattr(
                b.params.config.bubble_dew_method, "scale_" + short_name
            )
            scaling_method(b, overwrite=False)

            return

        if self.is_property_constructed("temperature_bubble"):
            bubble_dew_scaling(self, self.temperature_bubble)

        if self.is_property_constructed("temperature_dew"):
            bubble_dew_scaling(self, self.temperature_dew)

        if self.is_property_constructed("pressure_bubble"):
            bubble_dew_scaling(self, self.pressure_bubble)

        if self.is_property_constructed("pressure_dew"):
            bubble_dew_scaling(self, self.pressure_dew)

        # Scale log form constraints
        if self.is_property_constructed("log_mole_frac_comp"):
            for j, v in self.log_mole_frac_comp_eqn.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_comp[j],
                    default=1e3,
                    warning=True,
                    hint="for log_mole_frac_comp",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        # Activity is generally of similar order to mole fractions
        if self.is_property_constructed("log_act_phase_comp"):
            for (p, j), v in self.log_act_phase_comp_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mole_frac_phase_comp",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_act_phase_comp_apparent"):
            for (p, j), v in self.log_act_phase_comp_apparent_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp_apparent[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mole_frac_phase_comp_apparent",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_act_phase_comp_true"):
            for (p, j), v in self.log_act_phase_comp_true_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp_true[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mole_frac_phase_comp_true",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if (
            self.is_property_constructed("log_act_phase_solvents")
            and len(self.params.solvent_set) > 1
        ):
            for p, v in self.log_act_phase_solvents_eq.items():
                iscale.constraint_scaling_transform(v, 1e-3, overwrite=False)

        if self.is_property_constructed("log_conc_mol_phase_comp"):
            for (p, j), v in self.log_conc_mol_phase_comp_eq.items():
                sf_dens_mol = iscale.get_scaling_factor(
                    self.dens_mol_phase[p],
                    default=55e3,
                    warning=True,
                    hint="for log_conc_mol_phase_comp",
                )
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp[p, j],
                    default=1,
                    warning=True,
                    hint="for log_conc_mol_phase_comp",
                )
                iscale.constraint_scaling_transform(
                    v, sf_dens_mol * sf_x, overwrite=False
                )

        if self.is_property_constructed("log_conc_mol_phase_comp_true"):
            for (p, j), v in self.log_conc_mol_phase_comp_true_eq.items():
                sf_dens_mol = iscale.get_scaling_factor(
                    self.dens_mol_phase[p],
                    default=55e3,
                    warning=True,
                    hint="for log_conc_mol_phase_comp_true",
                )
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp_true[p, j],
                    default=1,
                    warning=True,
                    hint="for log_conc_mol_phase_comp_true",
                )
                iscale.constraint_scaling_transform(
                    v, sf_dens_mol * sf_x, overwrite=False
                )

        if self.is_property_constructed("log_mass_frac_phase_comp"):
            for (p, j), v in self.log_mass_frac_phase_comp_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mass_frac_phase_comp[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mass_frac_phase_comp",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_mass_frac_phase_comp_apparent"):
            for (p, j), v in self.log_mass_frac_phase_comp_apparent_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mass_frac_phase_comp_apparent[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mass_frac_phase_comp_apparent",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_mass_frac_phase_comp_true"):
            for (p, j), v in self.log_mass_frac_phase_comp_true_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mass_frac_phase_comp_true[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mass_frac_phase_comp_true",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_molality_phase_comp"):
            for (p, j), v in self.log_molality_phase_comp_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.molality_phase_comp[p, j],
                    default=1,
                    warning=True,
                    hint="for log_molality_phase_comp",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_molality_phase_comp_apparent"):
            for (p, j), v in self.log_molality_phase_comp_apparent_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.molality_phase_comp_apparent[p, j],
                    default=1,
                    warning=True,
                    hint="for log_molality_phase_comp_apparent",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_molality_phase_comp_true"):
            for (p, j), v in self.log_molality_phase_comp_true_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.molality_phase_comp_true[p, j],
                    default=1,
                    warning=True,
                    hint="for log_molality_phase_comp_true",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_mole_frac_phase_comp"):
            for (p, j), v in self.log_mole_frac_phase_comp_eqn.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp[p, j],
                    default=1e3,
                    warning=True,
                    hint="for log_mole_frac_phase_comp",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_mole_frac_phase_comp_apparent"):
            for (p, j), v in self.log_mole_frac_phase_comp_apparent_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp_apparent[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mole_frac_phase_comp_apparent",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

        if self.is_property_constructed("log_mole_frac_phase_comp_true"):
            for (p, j), v in self.log_mole_frac_phase_comp_true_eq.items():
                sf_x = iscale.get_scaling_factor(
                    self.mole_frac_phase_comp_true[p, j],
                    default=1e-3,
                    warning=True,
                    hint="for log_mole_frac_phase_comp_true",
                )
                iscale.constraint_scaling_transform(v, sf_x, overwrite=False)

    def components_in_phase(self, phase):
        """
        Generator method which yields components present in a given phase.
        As this methid is used only for property calcuations, it should use the
        true species sset if one exists.

        Args:
            phase - phase for which to yield components

        Yields:
            components present in phase.
        """
        if not self.params._electrolyte:
            component_list = self.component_list
            pc_set = self.phase_component_set
        else:
            config = self.params.get_phase(phase).config
            if (
                config.equation_of_state_options is not None
                and "property_basis" in config.equation_of_state_options
                and config.equation_of_state_options["property_basis"] == "apparent"
            ):
                component_list = self.params.apparent_species_set
                pc_set = self.params.apparent_phase_component_set
            else:
                component_list = self.params.true_species_set
                pc_set = self.params.true_phase_component_set

        for j in component_list:
            if (phase, j) in pc_set:
                yield j

    def get_mole_frac(self, phase=None):
        """
        Property calcuations generally depend on phase_component mole fractions
        for mixing rules, but in some cases there are multiple component lists
        to work from. This method is used to return the correct phase-component
        indexed mole fraction component for the given circumstances.

        Returns:
            mole fraction object
        """
        if not self.params._electrolyte:
            return self.mole_frac_phase_comp
        else:
            if phase is not None:
                config = self.params.get_phase(phase).config
                if (
                    config.equation_of_state_options is not None
                    and "property_basis" in config.equation_of_state_options
                    and config.equation_of_state_options["property_basis"] == "apparent"
                ):
                    return self.mole_frac_phase_comp_apparent
            return self.mole_frac_phase_comp_true

    # -------------------------------------------------------------------------
    # Bubble and Dew Points

    def _temperature_bubble(b):
        _temperature_pressure_bubble_dew(b, "temperature_bubble")

    def _log_mole_frac_tbub(b):
        _log_mole_frac_bubble_dew(b, "log_mole_frac_tbub")

    def _temperature_dew(b):
        _temperature_pressure_bubble_dew(b, "temperature_dew")

    def _log_mole_frac_tdew(b):
        _log_mole_frac_bubble_dew(b, "log_mole_frac_tdew")

    def _pressure_bubble(b):
        _temperature_pressure_bubble_dew(b, "pressure_bubble")

    def _log_mole_frac_pbub(b):
        _log_mole_frac_bubble_dew(b, "log_mole_frac_pbub")

    def _pressure_dew(b):
        _temperature_pressure_bubble_dew(b, "pressure_dew")

    def _log_mole_frac_pdew(b):
        _log_mole_frac_bubble_dew(b, "log_mole_frac_pdew")

    # -------------------------------------------------------------------------
    # Property Methods
    def _act_phase_comp(self):
        try:

            def rule_act_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.act_phase_comp(b, p, j)

            self.act_phase_comp = Expression(
                self.phase_component_set,
                doc="Component activity in each phase",
                rule=rule_act_phase_comp,
            )
        except AttributeError:
            self.del_component(self.act_phase_comp)
            raise

    def _act_phase_comp_true(self):
        try:

            def rule_act_phase_comp_true(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.act_phase_comp_true(b, p, j)

            self.act_phase_comp_true = Expression(
                self.params.true_phase_component_set,
                doc="Component activity in each phase",
                rule=rule_act_phase_comp_true,
            )
        except AttributeError:
            self.del_component(self.act_phase_comp_true)
            raise

    def _act_phase_comp_apparent(self):
        try:

            def rule_act_phase_comp_apparent(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.act_phase_comp_apparent(b, p, j)

            self.act_phase_comp_apparent = Expression(
                self.params.apparent_phase_component_set,
                doc="Component activity in each phase",
                rule=rule_act_phase_comp_apparent,
            )
        except AttributeError:
            self.del_component(self.act_phase_comp_apparent)
            raise

    def _act_coeff_phase_comp(self):
        try:

            def rule_act_coeff_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.act_coeff_phase_comp(b, p, j)

            self.act_coeff_phase_comp = Expression(
                self.phase_component_set,
                doc="Component activity coefficient in each phase",
                rule=rule_act_coeff_phase_comp,
            )
        except AttributeError:
            self.del_component(self.act_coeff_phase_comp)
            raise

    def _act_coeff_phase_comp_true(self):
        try:

            def rule_act_coeff_phase_comp_true(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.act_coeff_phase_comp_true(b, p, j)

            self.act_coeff_phase_comp_true = Expression(
                self.params.true_phase_component_set,
                doc="Component activity coefficient in each phase",
                rule=rule_act_coeff_phase_comp_true,
            )
        except AttributeError:
            self.del_component(self.act_coeff_phase_comp_true)
            raise

    def _act_coeff_phase_comp_apparent(self):
        try:

            def rule_act_coeff_phase_comp_apparent(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.act_coeff_phase_comp_apparent(b, p, j)

            self.act_coeff_phase_comp_apparent = Expression(
                self.params.apparent_phase_component_set,
                doc="Component activity coefficient in each phase",
                rule=rule_act_coeff_phase_comp_apparent,
            )
        except AttributeError:
            self.del_component(self.act_coeff_phase_comp_apparent)
            raise

    def _compress_fact_phase(self):
        try:

            def rule_Z_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.compress_fact_phase(b, p)

            self.compress_fact_phase = Expression(
                self.phase_list, doc="Compressibility of each phase", rule=rule_Z_phase
            )
        except AttributeError:
            self.del_component(self.compress_fact_phase)
            raise

    def _conc_mol_comp(self):
        try:

            def rule_conc_mol_comp(b, j):
                return b.dens_mol * b.mole_frac_comp[j]

            self.conc_mol_comp = Expression(
                self.component_list,
                rule=rule_conc_mol_comp,
                doc="Molar concentration of component",
            )
        except AttributeError:
            self.del_component(self.conc_mol_comp)
            raise

    def _conc_mol_phase_comp(self):
        try:

            def rule_conc_mol_phase_comp(b, p, j):
                return b.dens_mol_phase[p] * b.mole_frac_phase_comp[p, j]

            self.conc_mol_phase_comp = Expression(
                self.phase_component_set,
                rule=rule_conc_mol_phase_comp,
                doc="Molar concentration of component by phase",
            )
        except AttributeError:
            self.del_component(self.conc_mol_phase_comp)
            raise

    def _conc_mol_phase_comp_apparent(self):
        try:

            def rule_conc_mol_phase_comp_apparent(b, p, j):
                return b.dens_mol_phase[p] * b.mole_frac_phase_comp_apparent[p, j]

            self.conc_mol_phase_comp_apparent = Expression(
                self.params.apparent_phase_component_set,
                rule=rule_conc_mol_phase_comp_apparent,
                doc="Molar concentration of apparent component by phase",
            )
        except AttributeError:
            self.del_component(self.conc_mol_phase_comp_apparent)
            raise

    def _conc_mol_phase_comp_true(self):
        try:

            def rule_conc_mol_phase_comp_true(b, p, j):
                return b.dens_mol_phase[p] * b.mole_frac_phase_comp_true[p, j]

            self.conc_mol_phase_comp_true = Expression(
                self.params.true_phase_component_set,
                rule=rule_conc_mol_phase_comp_true,
                doc="Molar concentration of true component by phase",
            )
        except AttributeError:
            self.del_component(self.conc_mol_phase_comp_true)
            raise

    def _cp_mol(self):
        try:

            def rule_cp_mol(b):
                return sum(b.cp_mol_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.cp_mol = Expression(
                rule=rule_cp_mol, doc="Mixture molar heat capacity"
            )
        except AttributeError:
            self.del_component(self.cp_mol)
            raise

    def _cp_mol_phase(self):
        try:

            def rule_cp_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.cp_mol_phase(b, p)

            self.cp_mol_phase = Expression(self.phase_list, rule=rule_cp_mol_phase)
        except AttributeError:
            self.del_component(self.cp_mol_phase)
            raise

    def _cp_mol_phase_comp(self):
        try:

            def rule_cp_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.cp_mol_phase_comp(b, p, j)

            self.cp_mol_phase_comp = Expression(
                self.phase_component_set, rule=rule_cp_mol_phase_comp
            )
        except AttributeError:
            self.del_component(self.cp_mol_phase_comp)
            raise

    def _cv_mol(self):
        try:

            def rule_cv_mol(b):
                return sum(b.cv_mol_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.cv_mol = Expression(
                rule=rule_cv_mol, doc="Mixture molar heat capacity"
            )
        except AttributeError:
            self.del_component(self.cv_mol)
            raise

    def _cv_mol_phase(self):
        try:

            def rule_cv_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.cv_mol_phase(b, p)

            self.cv_mol_phase = Expression(self.phase_list, rule=rule_cv_mol_phase)
        except AttributeError:
            self.del_component(self.cv_mol_phase)
            raise

    def _cv_mol_phase_comp(self):
        try:

            def rule_cv_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.cv_mol_phase_comp(b, p, j)

            self.cv_mol_phase_comp = Expression(
                self.phase_component_set, rule=rule_cv_mol_phase_comp
            )
        except AttributeError:
            self.del_component(self.cv_mol_phase_comp)
            raise

    def _heat_capacity_ratio_phase(self):
        try:

            def rule_heat_capacity_ratio_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.heat_capacity_ratio_phase(b, p)

            self.heat_capacity_ratio_phase = Expression(
                self.phase_list,
                rule=rule_heat_capacity_ratio_phase,
                doc="Heat capacity ratio by phase",
            )
        except AttributeError:
            self.del_component(self.heat_capacity_ratio_phase)
            raise

    def _dens_mass(self):
        try:

            def rule_dens_mass(b):
                return sum(b.dens_mass_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.dens_mass = Expression(doc="Mixture mass density", rule=rule_dens_mass)
        except AttributeError:
            self.del_component(self.dens_mass)
            raise

    def _dens_mass_phase(self):
        try:

            def rule_dens_mass_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.dens_mass_phase(b, p)

            self.dens_mass_phase = Expression(
                self.phase_list,
                doc="Mass density of each phase",
                rule=rule_dens_mass_phase,
            )
        except AttributeError:
            self.del_component(self.dens_mass_phase)
            raise

    def _dens_mol(self):
        try:

            def rule_dens_mol(b):
                return sum(b.dens_mol_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.dens_mol = Expression(doc="Mixture molar density", rule=rule_dens_mol)
        except AttributeError:
            self.del_component(self.dens_mol)
            raise

    def _dens_mol_phase(self):
        try:

            def rule_dens_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.dens_mol_phase(b, p)

            self.dens_mol_phase = Expression(
                self.phase_list,
                doc="Molar density of each phase",
                rule=rule_dens_mol_phase,
            )
        except AttributeError:
            self.del_component(self.dens_mol_phase)
            raise

    def _diffus_phase_comp(self):
        try:

            def rule_diffus_phase_comp(b, p, j):
                cobj = b.params.get_component(j)
                if (
                    cobj.config.diffus_phase_comp is not None
                    and p in cobj.config.diffus_phase_comp
                ):
                    return cobj.config.diffus_phase_comp[p].return_expression(
                        b, p, j, b.temperature
                    )
                else:
                    return Expression.Skip

            self.diffus_phase_comp = Expression(
                self.phase_component_set,
                doc="Diffusivity for each phase-component pair",
                rule=rule_diffus_phase_comp,
            )
        except AttributeError:
            self.del_component(self.diffus_phase_comp)
            raise

    def _diffus_phase_comp_apparent(self):
        try:

            def rule_diffus_phase_comp_apparent(b, p, j):
                cobj = b.params.get_component(j)
                if (
                    cobj.config.diffus_phase_comp is not None
                    and p in cobj.config.diffus_phase_comp
                ):
                    return cobj.config.diffus_phase_comp[p].return_expression(
                        b, p, j, b.temperature
                    )
                else:
                    return Expression.Skip

            self.diffus_phase_comp_apparent = Expression(
                self.params.apparent_phase_component_set,
                doc="Diffusivity for apparent components in each phase",
                rule=rule_diffus_phase_comp_apparent,
            )
        except AttributeError:
            self.del_component(self.diffus_phase_comp_apparent)
            raise

    def _diffus_phase_comp_true(self):
        try:

            def rule_diffus_phase_comp_true(b, p, j):
                cobj = b.params.get_component(j)
                if (
                    cobj.config.diffus_phase_comp is not None
                    and p in cobj.config.diffus_phase_comp
                ):
                    return cobj.config.diffus_phase_comp[p].return_expression(
                        b, p, j, b.temperature
                    )
                else:
                    return Expression.Skip

            self.diffus_phase_comp_true = Expression(
                self.params.true_phase_component_set,
                doc="Diffusivity for true components in each phase",
                rule=rule_diffus_phase_comp_true,
            )
        except AttributeError:
            self.del_component(self.diffus_phase_comp_true)
            raise

    def _energy_internal_mol(self):
        try:

            def rule_energy_internal_mol(b):
                return sum(
                    b.energy_internal_mol_phase[p] * b.phase_frac[p]
                    for p in b.phase_list
                )

            self.energy_internal_mol = Expression(
                rule=rule_energy_internal_mol, doc="Mixture molar internal energy"
            )
        except AttributeError:
            self.del_component(self.energy_internal_mol)
            raise

    def _energy_internal_mol_phase(self):
        try:

            def rule_energy_internal_mol_phase(b, p):
                eos = b.params.get_phase(p).config.equation_of_state
                return eos.energy_internal_mol_phase(b, p)

            self.energy_internal_mol_phase = Expression(
                self.phase_list, rule=rule_energy_internal_mol_phase
            )
        except AttributeError:
            self.del_component(self.energy_internal_mol_phase)
            raise

    def _energy_internal_mol_phase_comp(self):
        try:

            def rule_energy_internal_mol_phase_comp(b, p, j):
                eos = b.params.get_phase(p).config.equation_of_state
                return eos.energy_internal_mol_phase_comp(b, p, j)

            self.energy_internal_mol_phase_comp = Expression(
                self.phase_component_set, rule=rule_energy_internal_mol_phase_comp
            )
        except AttributeError:
            self.del_component(self.energy_internal_mol_phase_comp)
            raise

    def _enth_mol(self):
        try:

            def rule_enth_mol(b):
                return sum(b.enth_mol_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.enth_mol = Expression(rule=rule_enth_mol, doc="Mixture molar enthalpy")
        except AttributeError:
            self.del_component(self.enth_mol)
            raise

    def _enth_mol_phase(self):
        try:

            def rule_enth_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.enth_mol_phase(b, p)

            self.enth_mol_phase = Expression(self.phase_list, rule=rule_enth_mol_phase)
        except AttributeError:
            self.del_component(self.enth_mol_phase)
            raise

    def _enth_mol_phase_comp(self):
        try:

            def rule_enth_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.enth_mol_phase_comp(b, p, j)

            self.enth_mol_phase_comp = Expression(
                self.phase_component_set, rule=rule_enth_mol_phase_comp
            )
        except AttributeError:
            self.del_component(self.enth_mol_phase_comp)
            raise

    def _entr_mol(self):
        try:

            def rule_entr_mol(b):
                return sum(b.entr_mol_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.entr_mol = Expression(rule=rule_entr_mol, doc="Mixture molar entropy")
        except AttributeError:
            self.del_component(self.entr_mol)
            raise

    def _entr_mol_phase(self):
        try:

            def rule_entr_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.entr_mol_phase(b, p)

            self.entr_mol_phase = Expression(self.phase_list, rule=rule_entr_mol_phase)
        except AttributeError:
            self.del_component(self.entr_mol_phase)
            raise

    def _entr_mol_phase_comp(self):
        try:

            def rule_entr_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.entr_mol_phase_comp(b, p, j)

            self.entr_mol_phase_comp = Expression(
                self.phase_component_set, rule=rule_entr_mol_phase_comp
            )
        except AttributeError:
            self.del_component(self.entr_mol_phase_comp)
            raise

    def _flow_mass(self):
        try:
            if self.get_material_flow_basis() == MaterialFlowBasis.mass:
                self.flow_mass = Expression(
                    expr=sum(self.flow_mass_comp[j] for j in self.component_list),
                    doc="Mass flow rate",
                )
            elif self.get_material_flow_basis() == MaterialFlowBasis.molar:
                self.flow_mass = Expression(
                    expr=self.mw * self.flow_mol, doc="Mass flow rate"
                )
            else:
                raise PropertyPackageError(
                    "{} Generic Property Package set to use unsupported "
                    "material flow basis: {}".format(
                        self.name, self.get_material_flow_basis()
                    )
                )
        except AttributeError:
            self.del_component(self.flow_mass)
            raise

    def _flow_mass_phase(self):
        try:

            def rule_flow_mass_phase(b, p):
                if b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return b.flow_mass * b.phase_frac[p]
                elif self.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return b.mw_phase[p] * b.flow_mol_phase[p]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_mass_phase = Expression(
                self.phase_list,
                doc="Mass flow rate of each phase",
                rule=rule_flow_mass_phase,
            )
        except AttributeError:
            self.del_component(self.flow_mass_phase)
            raise

    def _flow_mass_comp(self):
        try:

            def rule_flow_mass_comp(b, i):
                if b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return self.mass_frac_comp[i] * self.flow_mass
                elif b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return b.flow_mol_comp[i] * b.mw_comp[i]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_mass_comp = Expression(
                self.component_list,
                doc="Component mass flow rate",
                rule=rule_flow_mass_comp,
            )
        except AttributeError:
            self.del_component(self.flow_mass_comp)
            raise

    def _flow_mass_phase_comp(self):
        try:

            def rule_flow_mass_phase_comp(b, p, i):
                if b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    raise PropertyPackageError(
                        "{} Generic property Package set to use material flow "
                        "basis {}, but flow_mass_phase_comp was not created "
                        "by state definition."
                    )
                elif b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return b.flow_mol_phase_comp[p, i] * b.mw_comp[i]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_mass_phase_comp = Expression(
                self.phase_component_set,
                doc="Phase-component mass flow rate",
                rule=rule_flow_mass_phase_comp,
            )
        except AttributeError as e:
            print(e)
            self.del_component(self.flow_mass_phase_comp)
            raise

    def _flow_mol(self):
        try:
            if self.get_material_flow_basis() == MaterialFlowBasis.molar:
                self.flow_mol = Expression(
                    expr=sum(self.flow_mol_comp[j] for j in self.component_list),
                    doc="Total molar flow rate",
                )
            elif self.get_material_flow_basis() == MaterialFlowBasis.mass:
                self.flow_mol = Expression(
                    expr=self.flow_mass / self.mw, doc="Total molar flow rate"
                )
            else:
                raise PropertyPackageError(
                    "{} Generic Property Package set to use unsupported "
                    "material flow basis: {}".format(
                        self.name, self.get_material_flow_basis()
                    )
                )
        except AttributeError:
            self.del_component(self.flow_mol)
            raise

    def _flow_mol_phase(self):
        try:

            def rule_flow_mol_phase(b, p):
                if b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return b.flow_mol * b.phase_frac[p]
                elif b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return b.flow_mass_phase[p] / b.mw_phase[p]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_mol_phase = Expression(
                self.phase_list,
                doc="Molar flow rate of each phase",
                rule=rule_flow_mol_phase,
            )
        except AttributeError:
            self.del_component(self.flow_mol_phase)
            raise

    def _flow_mol_comp(self):
        try:

            def rule_flow_mol_comp(b, i):
                if b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return self.mole_frac_comp[i] * self.flow_mol
                elif b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return b.flow_mass_comp[i] / b.mw_comp[i]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_mol_comp = Expression(
                self.component_list,
                doc="Component molar flow rate",
                rule=rule_flow_mol_comp,
            )
        except AttributeError:
            self.del_component(self.flow_mol_comp)
            raise

    def _flow_mol_phase_comp(self):
        try:

            def rule_flow_mol_phase_comp(b, p, i):
                if b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    raise PropertyPackageError(
                        "{} Generic proeprty Package set to use material flow "
                        "basis {}, but flow_mol_phase_comp was not created "
                        "by state definition."
                    )
                elif b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return b.flow_mass_phase_comp[p, i] / b.mw_comp[i]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_mol_phase_comp = Expression(
                self.phase_component_set,
                doc="Phase-component molar flow rate",
                rule=rule_flow_mol_phase_comp,
            )
        except AttributeError:
            self.del_component(self.flow_mol_phase_comp)
            raise

    def _flow_vol(self):
        try:

            def _flow_vol_rule(b):
                if b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return self.flow_mol / self.dens_mol
                elif b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return self.flow_mass / self.dens_mass
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_vol = Expression(rule=_flow_vol_rule, doc="Volumetric flow rate")
        except AttributeError:
            self.del_component(self.flow_vol)
            raise

    def _flow_vol_phase(self):
        try:

            def rule_flow_vol_phase(b, p):
                if b.get_material_flow_basis() == MaterialFlowBasis.molar:
                    return b.flow_mol_phase[p] / b.dens_mol_phase[p]
                elif b.get_material_flow_basis() == MaterialFlowBasis.mass:
                    return b.flow_mass_phase[p] / b.dens_mass_phase[p]
                else:
                    raise PropertyPackageError(
                        "{} Generic Property Package set to use unsupported "
                        "material flow basis: {}".format(
                            self.name, self.get_material_flow_basis()
                        )
                    )

            self.flow_vol_phase = Expression(
                self.phase_list,
                doc="Volumetric flow rate of each phase",
                rule=rule_flow_vol_phase,
            )
        except AttributeError:
            self.del_component(self.flow_vol_phase)
            raise

    def _fug_phase_comp(self):
        try:

            def rule_fug_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.fug_phase_comp(b, p, j)

            self.fug_phase_comp = Expression(
                self.phase_component_set, rule=rule_fug_phase_comp
            )
        except AttributeError:
            self.del_component(self.fug_phase_comp)
            raise

    def _fug_coeff_phase_comp(self):
        try:

            def rule_fug_coeff_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.fug_coeff_phase_comp(b, p, j)

            self.fug_coeff_phase_comp = Expression(
                self.phase_component_set, rule=rule_fug_coeff_phase_comp
            )
        except AttributeError:
            self.del_component(self.fug_coeff_phase_comp)
            raise

    def _gibbs_mol(self):
        try:

            def rule_gibbs_mol(b):
                return sum(b.gibbs_mol_phase[p] * b.phase_frac[p] for p in b.phase_list)

            self.gibbs_mol = Expression(
                rule=rule_gibbs_mol, doc="Mixture molar Gibbs energy"
            )
        except AttributeError:
            self.del_component(self.gibbs_mol)
            raise

    def _gibbs_mol_phase(self):
        try:

            def rule_gibbs_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.gibbs_mol_phase(b, p)

            self.gibbs_mol_phase = Expression(
                self.phase_list, rule=rule_gibbs_mol_phase
            )
        except AttributeError:
            self.del_component(self.gibbs_mol_phase)
            raise

    def _gibbs_mol_phase_comp(self):
        try:

            def rule_gibbs_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.gibbs_mol_phase_comp(b, p, j)

            self.gibbs_mol_phase_comp = Expression(
                self.phase_component_set, rule=rule_gibbs_mol_phase_comp
            )
        except AttributeError:
            self.del_component(self.gibbs_mol_phase_comp)
            raise

    def _isentropic_speed_sound_phase(self):
        try:

            def rule_isentropic_speed_sound_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.isentropic_speed_sound_phase(b, p)

            self.isentropic_speed_sound_phase = Expression(
                self.phase_list,
                doc="Isentropic speed of sound in each phase",
                rule=rule_isentropic_speed_sound_phase,
            )
        except AttributeError:
            self.del_component(self.isentropic_speed_sound_phase)
            raise

    def _isothermal_speed_sound_phase(self):
        try:

            def rule_isothermal_speed_sound_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.isothermal_speed_sound_phase(b, p)

            self.isothermal_speed_sound_phase = Expression(
                self.phase_list,
                doc="Isothermal speed of sound in each phase",
                rule=rule_isothermal_speed_sound_phase,
            )
        except AttributeError:
            self.del_component(self.isothermal_speed_sound_phase)
            raise

    def _henry(self):
        try:

            def henry_rule(b, p, j):
                cobj = b.params.get_component(j)
                if (
                    cobj.config.henry_component is not None
                    and p in cobj.config.henry_component
                ):
                    return cobj.config.henry_component[p]["method"].return_expression(
                        b, p, j, b.temperature
                    )
                else:
                    return Expression.Skip

            self.henry = Expression(self.phase_component_set, rule=henry_rule)
        except AttributeError:
            self.del_component(self.henry)
            raise

    def _mass_frac_phase_comp(self):
        # If this doesn't exist yet, assume it needs to be built from mole_frac
        try:

            def rule_mass_frac(b, p, j):
                return (
                    b.mw_comp[j]
                    * b.mole_frac_phase_comp[p, j]
                    / sum(
                        b.mw_comp[k] * b.mole_frac_phase_comp[p, k]
                        for k in b.component_list
                        if (p, k) in b.phase_component_set
                    )
                )

            self.mass_frac_phase_comp = Expression(
                self.phase_component_set,
                rule=rule_mass_frac,
                doc="Phase-component mass fractions",
            )
        except AttributeError:
            self.del_component(self.mass_frac_phase_comp)
            raise

    def _mass_frac_phase_comp_apparent(self):
        # If this doesn't exist yet, assume it needs to be built from mole_frac
        try:

            def rule_mass_frac_apparent(b, p, j):
                return (
                    b.mw_comp[j]
                    * b.mole_frac_phase_comp_apparent[p, j]
                    / sum(
                        b.mw_comp[k] * b.mole_frac_phase_comp_apparent[p, k]
                        for k in b.params.apparent_species_set
                        if (p, k) in b.params.apparent_phase_component_set
                    )
                )

            self.mass_frac_phase_comp_apparent = Expression(
                self.params.apparent_phase_component_set,
                rule=rule_mass_frac_apparent,
                doc="Phase-component mass fractions",
            )
        except AttributeError:
            self.del_component(self.mass_frac_phase_comp_apparent)
            raise

    def _mass_frac_phase_comp_true(self):
        # If this doesn't exist yet, assume it needs to be built from mole_frac
        try:

            def rule_mass_frac_true(b, p, j):
                return (
                    b.mw_comp[j]
                    * b.mole_frac_phase_comp_true[p, j]
                    / sum(
                        b.mw_comp[k] * b.mole_frac_phase_comp_true[p, k]
                        for k in b.params.true_species_set
                        if (p, k) in b.paramas.ture_phase_component_set
                    )
                )

            self.mass_frac_phase_comp_true = Expression(
                self.params.true_phase_component_set,
                rule=rule_mass_frac_true,
                doc="Phase-component mass fractions",
            )
        except AttributeError:
            self.del_component(self.mass_frac_phase_comp_true)
            raise

    def _mw(self):
        try:
            self.mw = Expression(
                doc="Average molecular weight",
                expr=sum(
                    self.phase_frac[p]
                    * sum(
                        self.mole_frac_phase_comp[p, j]
                        * self.params.get_component(j).mw
                        if (p, j) in self.phase_component_set
                        else 0
                        for j in self.component_list
                    )
                    for p in self.phase_list
                ),
            )
        except AttributeError:
            self.del_component(self.mw)
            raise

    def _molality_phase_comp(self):
        try:

            def rule_molality(b, p, j):
                pobj = b.params.get_phase(p)
                if not isinstance(pobj, LiquidPhase):
                    # Molality is defined per mass of solvent so only makes
                    # sense in liquid phases
                    return Expression.Skip
                else:
                    return b.mole_frac_phase_comp[p, j] / sum(
                        b.mw_comp[k] * b.mole_frac_phase_comp[p, k]
                        for k in b.params.solvent_set
                    )

            self.molality_phase_comp = Expression(
                self.phase_component_set,
                rule=rule_molality,
                doc="Component molalitites",
            )
        except AttributeError:
            self.del_component(self.molality_phase_comp)
            raise

    def _molality_phase_comp_apparent(self):
        try:

            def rule_molality_apparent(b, p, j):
                pobj = b.params.get_phase(p)
                if not isinstance(pobj, LiquidPhase):
                    # Molality is defined per mass of solvent so only makes
                    # sense in liquid phases
                    return Expression.Skip
                else:
                    return b.mole_frac_phase_comp_apparent[p, j] / sum(
                        b.mw_comp[k] * b.mole_frac_phase_comp_apparent[p, k]
                        for k in b.params.solvent_set
                    )

            self.molality_phase_comp_apparent = Expression(
                self.phase_component_set,
                rule=rule_molality_apparent,
                doc="Component molalitites",
            )
        except AttributeError:
            self.del_component(self.molality_phase_comp_apparent)
            raise

    def _molality_phase_comp_true(self):
        try:

            def rule_molality_true(b, p, j):
                pobj = b.params.get_phase(p)
                if not isinstance(pobj, LiquidPhase):
                    # Molality is defined per mass of solvent so only makes
                    # sense in liquid phases
                    return Expression.Skip
                else:
                    return b.mole_frac_phase_comp_true[p, j] / sum(
                        b.mw_comp[k] * b.mole_frac_phase_comp_true[p, k]
                        for k in b.params.solvent_set
                    )

            self.molality_phase_comp_true = Expression(
                self.phase_component_set,
                rule=rule_molality_true,
                doc="Component molalitites",
            )
        except AttributeError:
            self.del_component(self.molality_phase_comp_true)
            raise

    def _mole_frac_comp(self):
        """If mole_frac_comp not state var assume mole_frac_phase_comp is"""
        try:

            def rule_mole_frac_comp(b, i):
                return sum(
                    b.phase_frac[p] * b.mole_frac_phase_comp[p, i] for p in b.phase_list
                )

            self.mole_frac_comp = Expression(
                self.component_list,
                doc="Mole fraction of each component",
                rule=rule_mole_frac_comp,
            )
        except AttributeError:
            self.del_component(self.mole_frac_comp)
            raise

    def _mw_comp(self):
        try:

            def rule_mw_comp(b, j):
                return b.params.get_component(j).mw

            self.mw_comp = Expression(
                self.component_list,
                doc="Molecular weight of each component",
                rule=rule_mw_comp,
            )
        except AttributeError:
            self.del_component(self.mw_comp)
            raise

    def _mw_phase(self):
        try:

            def rule_mw_phase(b, p):
                return sum(
                    b.mole_frac_phase_comp[p, j] * b.params.get_component(j).mw
                    if (p, j) in b.phase_component_set
                    else 0
                    for j in b.component_list
                )

            self.mw_phase = Expression(
                self.phase_list,
                doc="Average molecular weight of each phase",
                rule=rule_mw_phase,
            )
        except AttributeError:
            self.del_component(self.mw_phase)
            raise

    def _pressure_phase_comp(self):
        try:

            def rule_partial_pressure(b, p, i):
                return b.mole_frac_phase_comp[p, i] * b.pressure

            self.pressure_phase_comp = Expression(
                self.phase_component_set,
                doc="Component partial pressure in phase",
                rule=rule_partial_pressure,
            )
        except AttributeError:
            self.del_component(self.pressure_phase_comp)
            raise

    def _pressure_phase_comp_true(self):
        try:

            def rule_partial_pressure(b, p, i):
                return b.mole_frac_phase_comp_true[p, i] * b.pressure

            self.pressure_phase_comp_true = Expression(
                self.phase_component_set,
                doc="Component partial pressure in phase",
                rule=rule_partial_pressure,
            )
        except AttributeError:
            self.del_component(self.pressure_phase_comp_true)
            raise

    def _pressure_phase_comp_apparent(self):
        try:

            def rule_partial_pressure(b, p, i):
                return b.mole_frac_phase_comp_apparent[p, i] * b.pressure

            self.pressure_phase_comp_apparent = Expression(
                self.phase_component_set,
                doc="Component partial pressure in phase",
                rule=rule_partial_pressure,
            )
        except AttributeError:
            self.del_component(self.pressure_phase_comp_apparent)
            raise

    def _pressure_osm_phase(self):
        try:

            def rule_posm_phase(b, p):
                pobj = b.params.get_phase(p)
                if isinstance(pobj, LiquidPhase):
                    p_config = pobj.config
                    return p_config.equation_of_state.pressure_osm_phase(b, p)
                else:
                    return Expression.Skip

            self.pressure_osm_phase = Expression(
                self.phase_list,
                doc="Osmotic pressure in each phase",
                rule=rule_posm_phase,
            )
        except AttributeError:
            self.del_component(self.pressure_osm_phase)
            raise

    def _pressure_sat_comp(self):
        try:

            def rule_pressure_sat_comp(b, j):
                cobj = b.params.get_component(j)
                try:
                    return get_method(b, "pressure_sat_comp", j)(b, cobj, b.temperature)
                except GenericPropertyPackageError:
                    # There is the possibility this is a Henry component
                    if cobj.config.henry_component is not None:
                        # Assume it is a Henry component and skip
                        _log.debug(
                            "{} Component {} does not have a method for"
                            " pressure_sat_comp, but is marked as being"
                            " Henry component in at least one phase. "
                            "It will be assumed that saturation "
                            "pressure is not required for this "
                            "component.".format(b.name, j)
                        )
                        return Expression.Skip
                    else:
                        raise

            self.pressure_sat_comp = Expression(
                self.component_list, rule=rule_pressure_sat_comp
            )
        except AttributeError:
            self.del_component(self.pressure_sat_comp)
            raise

    def _surf_tens_phase(self):
        try:

            def rule_surf_tens_phase(b, p):
                pobj = b.params.get_phase(p)
                if isinstance(pobj, LiquidPhase):
                    return get_phase_method(b, "surf_tens_phase", p)(b, p)
                else:
                    return Expression.Skip

            self.surf_tens_phase = Expression(
                self.phase_list,
                doc="Surface tension for each phase",
                rule=rule_surf_tens_phase,
            )
        except AttributeError:
            self.del_component(self.surf_tens_phase)
            raise

    def _therm_cond_phase(self):
        try:

            def rule_therm_cond_phase(b, p):
                return get_phase_method(b, "therm_cond_phase", p)(b, p)

            self.therm_cond_phase = Expression(
                self.phase_list,
                doc="Thermal conductivity for each phase",
                rule=rule_therm_cond_phase,
            )
        except AttributeError:
            self.del_component(self.therm_cond_phase)
            raise

    def _visc_d_phase(self):
        try:

            def rule_visc_d_phase(b, p):
                return get_phase_method(b, "visc_d_phase", p)(b, p)

            self.visc_d_phase = Expression(
                self.phase_list,
                doc="Dynamic viscosity for each phase",
                rule=rule_visc_d_phase,
            )
        except AttributeError:
            self.del_component(self.visc_d_phase)
            raise

    def _vol_mol_phase(self):
        try:

            def rule_vol_mol_phase(b, p):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.vol_mol_phase(b, p)

            self.vol_mol_phase = Expression(
                self.phase_list,
                doc="Molar volume of each phase",
                rule=rule_vol_mol_phase,
            )
        except AttributeError:
            self.del_component(self.vol_mol_phase)
            raise

    def _vol_mol_phase_comp(self):
        try:

            def rule_vol_mol_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return p_config.equation_of_state.vol_mol_phase_comp(b, p, j)

            self.vol_mol_phase_comp = Expression(
                self.phase_list,
                self.component_list,
                doc="Partial molar volume of component j",
                rule=rule_vol_mol_phase_comp,
            )
        except AttributeError:
            self.del_component(self.vol_mol_phase_comp)
            raise

    def _dh_rxn(self):
        def dh_rule(b, r):
            rblock = getattr(b.params, "reaction_" + r)

            carg = b.params.config.inherent_reactions[r]

            return carg["heat_of_reaction"].return_expression(
                b, rblock, r, b.temperature
            )

        self.dh_rxn = Expression(
            self.params.inherent_reaction_idx,
            doc="Specific heat of reaction for inherent reactions",
            rule=dh_rule,
        )

    def _log_act_phase_comp(self):
        try:
            self.log_act_phase_comp = Var(
                self.phase_component_set,
                initialize=1,
                bounds=(-50, 1),
                units=pyunits.dimensionless,
                doc="Log of activity of component by phase",
            )

            def rule_log_act_phase_comp(b, p, j):
                p_config = b.params.get_phase(p).config
                return exp(
                    b.log_act_phase_comp[p, j]
                ) == p_config.equation_of_state.act_phase_comp(b, p, j)

            self.log_act_phase_comp_eq = Constraint(
                self.phase_component_set,
                doc="Natural log of component activity in each phase",
                rule=rule_log_act_phase_comp,
            )
        except AttributeError:
            self.del_component(self.log_act_phase_comp)
            self.del_component(self.log_act_phase_comp_eq)
            raise

    def _log_act_phase_solvents(self):
        if len(self.params.solvent_set) == 1:
            self.log_act_phase_solvents = Reference(
                self.log_act_phase_comp[:, self.params.solvent_set.first()]
            )
        else:
            try:
                self.log_act_phase_solvents = Var(
                    self.phase_list,
                    initialize=1,
                    bounds=(-50, 1),
                    units=pyunits.dimensionless,
                    doc="Log of activities summed across solvents by phase",
                )

                def rule_log_act_phase_solvents(b, p):
                    p_obj = b.params.get_phase(p)
                    p_config = b.params.get_phase(p).config
                    if not isinstance(p_obj, LiquidPhase):
                        return Expression.Skip
                    else:
                        return exp(b.log_act_phase_solvents[p]) == sum(
                            p_config.equation_of_state.act_phase_comp(b, p, j)
                            for j in b.params.solvent_set
                        )

                self.log_act_phase_solvents_eq = Constraint(
                    self.phase_list,
                    doc="Natural log of solvent activity in each phase",
                    rule=rule_log_act_phase_solvents,
                )
            except AttributeError:
                self.del_component(self.log_act_phase_solvents)
                self.del_component(self.log_act_phase_solvents_eq)
                raise

    def _log_act_phase_comp_true(self):
        try:
            self.log_act_phase_comp_true = Var(
                self.params.true_phase_component_set,
                initialize=1,
                bounds=(-50, 1),
                units=pyunits.dimensionless,
                doc="Log of activity of component by phase",
            )

            def rule_log_act_phase_comp_true(b, p, j):
                p_config = b.params.get_phase(p).config
                return exp(
                    b.log_act_phase_comp_true[p, j]
                ) == p_config.equation_of_state.act_phase_comp_true(b, p, j)

            self.log_act_phase_comp_true_eq = Constraint(
                self.params.true_phase_component_set,
                doc="Natural log of component activity in each phase",
                rule=rule_log_act_phase_comp_true,
            )
        except AttributeError:
            self.del_component(self.log_act_phase_comp_true)
            self.del_component(self.log_act_phase_comp_true_eq)
            raise

    def _log_act_phase_comp_apparent(self):
        try:
            self.log_act_phase_comp_apparent = Var(
                self.params.apparent_phase_component_set,
                initialize=1,
                bounds=(-50, 1),
                units=pyunits.dimensionless,
                doc="Log of activity of component by phase",
            )

            def rule_log_act_phase_comp_apparent(b, p, j):
                p_config = b.params.get_phase(p).config
                return exp(
                    b.log_act_phase_comp_apparent[p, j]
                ) == p_config.equation_of_state.act_phase_comp_apparent(b, p, j)

            self.log_act_phase_comp_apparent_eq = Constraint(
                self.params.apparent_phase_component_set,
                doc="Natural log of component activity in each phase",
                rule=rule_log_act_phase_comp_apparent,
            )
        except AttributeError:
            self.del_component(self.log_act_phase_comp_apparent)
            self.del_component(self.log_act_phase_comp_apparent_eq)
            raise

    def _log_conc_mol_phase_comp(self):
        try:
            self.log_conc_mol_phase_comp = Var(
                self.phase_component_set,
                initialize=1,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of molar concentration of component by phase",
            )

            def rule_log_conc_mol_phase_comp(b, p, j):
                return exp(b.log_conc_mol_phase_comp[p, j]) == (
                    b.conc_mol_phase_comp[p, j]
                    / pyunits.get_units(b.conc_mol_phase_comp[p, j])
                )

            self.log_conc_mol_phase_comp_eq = Constraint(
                self.phase_component_set,
                rule=rule_log_conc_mol_phase_comp,
                doc="Constraint for log of molar concentration",
            )
        except AttributeError:
            self.del_component(self.log_conc_mol_phase_comp)
            self.del_component(self.log_conc_mol_phase_comp_eq)
            raise

    def _log_conc_mol_phase_comp_true(self):
        try:
            self.log_conc_mol_phase_comp_true = Var(
                self.params.true_phase_component_set,
                initialize=1,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of molar concentration of component by phase",
            )

            def rule_log_conc_mol_phase_comp_true(b, p, j):
                return exp(b.log_conc_mol_phase_comp_true[p, j]) == (
                    b.conc_mol_phase_comp_true[p, j]
                    / pyunits.get_units(b.conc_mol_phase_comp_true[p, j])
                )

            self.log_conc_mol_phase_comp_true_eq = Constraint(
                self.params.true_phase_component_set,
                rule=rule_log_conc_mol_phase_comp_true,
                doc="Constraint for log of molar concentration",
            )
        except AttributeError:
            self.del_component(self.log_conc_mol_phase_comp_true)
            self.del_component(self.log_conc_mol_phase_comp_true_eq)
            raise

    def _log_mass_frac_phase_comp(self):
        try:
            self.log_mass_frac_phase_comp = Var(
                self.phase_component_set,
                initialize=0,
                bounds=(-50, 0),
                units=pyunits.dimensionless,
                doc="Log of mass fractions of component by phase",
            )

            def rule_log_mass_frac_phase_comp(b, p, j):
                return exp(b.log_mass_frac_phase_comp[p, j]) == (
                    b.mass_frac_phase_comp[p, j]
                )

            self.log_mass_frac_phase_comp_eq = Constraint(
                self.phase_component_set,
                rule=rule_log_mass_frac_phase_comp,
                doc="Constraint for log of mass fractions",
            )
        except AttributeError:
            self.del_component(self.log_mass_frac_phase_comp)
            self.del_component(self.log_mass_frac_phase_comp_eq)
            raise

    def _log_mass_frac_phase_comp_apparent(self):
        try:
            self.log_mass_frac_phase_comp_apparent = Var(
                self.params.apparent_phase_component_set,
                initialize=0,
                bounds=(-50, 0),
                units=pyunits.dimensionless,
                doc="Log of mass fractions of component by phase",
            )

            def rule_log_mass_frac_phase_comp_appr(b, p, j):
                return exp(b.log_mass_frac_phase_comp_apparent[p, j]) == (
                    b.mass_frac_phase_comp_apparent[p, j]
                )

            self.log_mass_frac_phase_comp_apparent_eq = Constraint(
                self.params.apparent_phase_component_set,
                rule=rule_log_mass_frac_phase_comp_appr,
                doc="Constraint for log of mass fractions",
            )
        except AttributeError:
            self.del_component(self.log_mass_frac_phase_comp_apparent)
            self.del_component(self.log_mass_frac_phase_comp_apparent_eq)
            raise

    def _log_mass_frac_phase_comp_true(self):
        try:
            self.log_mass_frac_phase_comp_true = Var(
                self.params.true_phase_component_set,
                initialize=0,
                bounds=(-50, 0),
                units=pyunits.dimensionless,
                doc="Log of mass fractions of component by phase",
            )

            def rule_log_mass_frac_phase_comp_true(b, p, j):
                return exp(b.log_mass_frac_phase_comp_true[p, j]) == (
                    b.mass_frac_phase_comp_true[p, j]
                )

            self.log_mass_frac_phase_comp_true_eq = Constraint(
                self.params.true_phase_component_set,
                rule=rule_log_mass_frac_phase_comp_true,
                doc="Constraint for log of mass fractions",
            )
        except AttributeError:
            self.del_component(self.log_mass_frac_phase_comp_true)
            self.del_component(self.log_mass_frac_phase_comp_true_eq)
            raise

    def _log_molality_phase_comp(self):
        try:
            self.log_molality_phase_comp = Var(
                self.phase_component_set,
                initialize=1,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of molality of component by phase",
            )

            def rule_log_molality_phase_comp(b, p, j):
                pobj = b.params.get_phase(p)
                if not isinstance(pobj, LiquidPhase):
                    # Molality is defined per mass of solvent so only makes
                    # sense in liquid phases
                    return Expression.Skip
                return exp(b.log_molality_phase_comp[p, j]) == (
                    b.molality_phase_comp[p, j]
                )

            self.log_molality_phase_comp_eq = Constraint(
                self.phase_component_set,
                rule=rule_log_molality_phase_comp,
                doc="Constraint for log of molality",
            )
        except AttributeError:
            self.del_component(self.log_molality_phase_comp)
            self.del_component(self.log_molality_phase_comp_eq)
            raise

    def _log_molality_phase_comp_apparent(self):
        try:
            self.log_molality_phase_comp_apparent = Var(
                self.params.apparent_phase_component_set,
                initialize=1,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of molality of component by phase",
            )

            def rule_log_molality_phase_comp_appr(b, p, j):
                pobj = b.params.get_phase(p)
                if not isinstance(pobj, LiquidPhase):
                    # Molality is defined per mass of solvent so only makes
                    # sense in liquid phases
                    return Expression.Skip
                return exp(b.log_molality_phase_comp_apparent[p, j]) == (
                    b.molality_phase_comp_apparent[p, j]
                )

            self.log_molality_phase_comp_apparent_eq = Constraint(
                self.params.apparent_phase_component_set,
                rule=rule_log_molality_phase_comp_appr,
                doc="Constraint for log of molality",
            )
        except AttributeError:
            self.del_component(self.log_molality_phase_comp_apparent)
            self.del_component(self.log_molality_phase_comp_apparent_eq)
            raise

    def _log_molality_phase_comp_true(self):
        try:
            self.log_molality_phase_comp_true = Var(
                self.params.true_phase_component_set,
                initialize=1,
                units=pyunits.dimensionless,
                doc="Log of molality of component by phase",
            )

            def rule_log_molality_phase_comp_true(b, p, j):
                pobj = b.params.get_phase(p)
                if not isinstance(pobj, LiquidPhase):
                    # Molality is defined per mass of solvent so only makes
                    # sense in liquid phases
                    return Expression.Skip
                return exp(b.log_molality_phase_comp_true[p, j]) == (
                    b.molality_phase_comp_true[p, j]
                )

            self.log_molality_phase_comp_true_eq = Constraint(
                self.params.true_phase_component_set,
                rule=rule_log_molality_phase_comp_true,
                doc="Constraint for log of molality",
            )
        except AttributeError:
            self.del_component(self.log_molality_phase_comp_true)
            self.del_component(self.log_molality_phase_comp_true_eq)
            raise

    def _log_mole_frac_comp(self):
        try:
            self.log_mole_frac_comp = Var(
                self.component_list,
                initialize=-log(len(self.component_list)),
                bounds=(-100, log(1.001)),
                units=pyunits.dimensionless,
                doc="Log of mole fractions of component",
            )

            def rule_log_mole_frac_comp(b, j):
                return exp(b.log_mole_frac_comp[j]) == (b.mole_frac_comp[j])

            self.log_mole_frac_comp_eqn = Constraint(
                self.component_list,
                rule=rule_log_mole_frac_comp,
                doc="Constraint for log of mole fractions",
            )
        except AttributeError:
            self.del_component(self.log_mole_frac_comp)
            self.del_component(self.log_mole_frac_comp_eqn)
            raise

    def _log_mole_frac_phase_comp(self):
        try:
            self.log_mole_frac_phase_comp = Var(
                self.phase_component_set,
                initialize=0,
                bounds=(-100, log(1.001)),
                units=pyunits.dimensionless,
                doc="Log of mole fractions of component by phase",
            )

            def rule_log_mole_frac_phase_comp(b, p, j):
                return exp(b.log_mole_frac_phase_comp[p, j]) == (
                    b.mole_frac_phase_comp[p, j]
                )

            self.log_mole_frac_phase_comp_eqn = Constraint(
                self.phase_component_set,
                rule=rule_log_mole_frac_phase_comp,
                doc="Constraint for log of mole fractions",
            )
        except AttributeError:
            self.del_component(self.log_mole_frac_phase_comp)
            self.del_component(self.log_mole_frac_phase_comp_eqn)
            raise

    def _log_mole_frac_phase_comp_apparent(self):
        try:
            self.log_mole_frac_phase_comp_apparent = Var(
                self.params.apparent_phase_component_set,
                initialize=0,
                bounds=(-50, 0),
                units=pyunits.dimensionless,
                doc="Log of mole fractions of component by phase",
            )

            def rule_log_mole_frac_phase_comp_appr(b, p, j):
                return exp(b.log_mole_frac_phase_comp_apparent[p, j]) == (
                    b.mole_frac_phase_comp_apparent[p, j]
                )

            self.log_mole_frac_phase_comp_apparent_eq = Constraint(
                self.params.apparent_phase_component_set,
                rule=rule_log_mole_frac_phase_comp_appr,
                doc="Constraint for log of mole fractions",
            )
        except AttributeError:
            self.del_component(self.log_mole_frac_phase_comp_apparent)
            self.del_component(self.log_mole_frac_phase_comp_apparent_eq)
            raise

    def _log_mole_frac_phase_comp_true(self):
        try:
            self.log_mole_frac_phase_comp_true = Var(
                self.params.true_phase_component_set,
                initialize=0,
                bounds=(-50, 0),
                units=pyunits.dimensionless,
                doc="Log of mole fractions of component by phase",
            )

            def rule_log_mole_frac_phase_comp_true(b, p, j):
                return exp(b.log_mole_frac_phase_comp_true[p, j]) == (
                    b.mole_frac_phase_comp_true[p, j]
                )

            self.log_mole_frac_phase_comp_true_eq = Constraint(
                self.params.true_phase_component_set,
                rule=rule_log_mole_frac_phase_comp_true,
                doc="Constraint for log of mole fractions",
            )
        except AttributeError:
            self.del_component(self.log_mole_frac_phase_comp_true)
            self.del_component(self.log_mole_frac_phase_comp_true_eq)
            raise

    def _log_pressure_phase_comp(self):
        try:
            self.log_pressure_phase_comp = Var(
                self.phase_component_set,
                initialize=12,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of partial pressures of component by phase",
            )

            def rule_log_pressure_phase_comp(b, p, j):
                return exp(b.log_pressure_phase_comp[p, j]) == (
                    b.pressure_phase_comp[p, j]
                )

            self.log_pressure_phase_comp_eq = Constraint(
                self.phase_component_set,
                rule=rule_log_pressure_phase_comp,
                doc="Constraint for log of partial pressures",
            )
        except AttributeError:
            self.del_component(self.log_pressure_phase_comp)
            self.del_component(self.log_pressure_phase_comp_eq)
            raise

    def _log_pressure_phase_comp_apparent(self):
        try:
            self.log_pressure_phase_comp_apparent = Var(
                self.params.apparent_phase_component_set,
                initialize=12,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of partial pressures of component by phase",
            )

            def rule_log_pressure_phase_comp_appr(b, p, j):
                return exp(b.log_pressure_phase_comp_apparent[p, j]) == (
                    b.pressure_phase_comp_apparent[p, j]
                )

            self.log_pressure_phase_comp_apparent_eq = Constraint(
                self.params.apparent_phase_component_set,
                rule=rule_log_pressure_phase_comp_appr,
                doc="Constraint for log of partial pressures",
            )
        except AttributeError:
            self.del_component(self.log_pressure_phase_comp_apparent)
            self.del_component(self.log_pressure_phase_comp_apparent_eq)
            raise

    def _log_pressure_phase_comp_true(self):
        try:
            self.log_pressure_phase_comp_true = Var(
                self.params.true_phase_component_set,
                initialize=12,
                bounds=(-50, None),
                units=pyunits.dimensionless,
                doc="Log of partial pressures of component by phase",
            )

            def rule_log_pressure_phase_comp_true(b, p, j):
                return exp(b.log_pressure_phase_comp_true[p, j]) == (
                    b.pressure_phase_comp_true[p, j]
                )

            self.log_pressure_phase_comp_true_eq = Constraint(
                self.params.true_phase_component_set,
                rule=rule_log_pressure_phase_comp_true,
                doc="Constraint for log of partial pressures",
            )
        except AttributeError:
            self.del_component(self.log_pressure_phase_comp_true)
            self.del_component(self.log_pressure_phase_comp_true_eq)
            raise

    def _log_k_eq(self):
        try:
            self.log_k_eq = Var(
                self.params.inherent_reaction_idx,
                initialize=1,
                doc="Log of equilibrium constant for inherent reactions",
            )

            def rule_log_k_eq(b, r):
                rblock = getattr(b.params, "reaction_" + r)

                carg = b.params.config.inherent_reactions[r]
                return carg["equilibrium_constant"].return_log_expression(
                    b, rblock, r, b.temperature
                )

            self.log_k_eq_constraint = Constraint(
                self.params.inherent_reaction_idx,
                rule=rule_log_k_eq,
                doc="Constraint for log of equilibrium constant",
            )
        except AttributeError:
            self.del_component(self.log_k_eq)
            self.del_component(self.log_k_eq_constraint)
            raise


def _raise_dev_burnt_toast():
    raise BurntToast(
        "Users shouldn't be calling this function. "
        "If you're a dev, you know what you did."
    )


def _valid_VL_component_list(blk, pp):
    raoult_comps = []
    henry_comps = []
    # Only need to do this for V-L pairs, so check
    pparams = blk.params
    if (
        pparams.get_phase(pp[0]).is_liquid_phase()
        and pparams.get_phase(pp[1]).is_vapor_phase()
    ) or (
        pparams.get_phase(pp[0]).is_vapor_phase()
        and pparams.get_phase(pp[1]).is_liquid_phase()
    ):

        for j in blk.component_list:
            if (pp[0], j) in blk.phase_component_set and (
                pp[1],
                j,
            ) in blk.phase_component_set:
                cobj = blk.params.get_component(j)
                if cobj.config.henry_component is not None and (
                    pp[0] in cobj.config.henry_component
                    or pp[1] in cobj.config.henry_component
                ):
                    henry_comps.append(j)
                else:
                    raoult_comps.append(j)

    return raoult_comps, henry_comps


def _temperature_pressure_bubble_dew(b, name):
    #  temperature/pressure bubble/dew
    splt = name.split("_")
    docstring_var = " point "
    docstring_mole_frac = " mole fractions at "
    if splt[1] == "dew":
        abbrv = "dew"
        docstring_var = "Dew" + docstring_var
        docstring_mole_frac = "Liquid" + docstring_mole_frac + "at dew "
    elif splt[1] == "bubble":
        abbrv = "bub"
        docstring_var = "Bubble" + docstring_var
        docstring_mole_frac = "Vapor" + docstring_mole_frac + "at bubble "
    else:
        _raise_dev_burnt_toast()

    if splt[0] == "temperature":
        abbrv = "t" + abbrv
        bounds = (b.temperature.lb, b.temperature.ub)
        units = b.params.get_metadata().default_units.TEMPERATURE
    elif splt[0] == "pressure":
        abbrv = "p" + abbrv
        bounds = (b.pressure.lb, b.pressure.ub)
        units_meta = b.params.get_metadata().derived_units
        units = units_meta.PRESSURE
    else:
        _raise_dev_burnt_toast()

    docstring_var += splt[0]
    docstring_mole_frac += splt[0]

    if b.params.config.bubble_dew_method is None:
        raise GenericPropertyPackageError(b, name)

    try:
        setattr(
            b,
            name,
            Var(b.params._pe_pairs, doc=docstring_var, bounds=bounds, units=units),
        )
        setattr(
            b,
            "_mole_frac_" + abbrv,
            Var(
                b.params._pe_pairs,
                b.component_list,
                initialize=1 / len(b.component_list),
                bounds=(0, None),
                doc=docstring_mole_frac,
                units=None,
            ),
        )

        tmp = getattr(b.params.config.bubble_dew_method, name)
        tmp(b)

    except AttributeError:
        if hasattr(b, name):
            var = getattr(b, name)
            b.del_component(var)
        if hasattr(b, "_mole_frac_" + abbrv):
            mole_frac = getattr(b, "_mole_frac_" + abbrv)
            b.del_component(mole_frac)
        raise


def _log_mole_frac_bubble_dew(b, name):
    abbrv = name.split("_")[-1]
    docstring_var = " log mole fractions at "
    docstring_eqn = "Constraint for log of "

    if abbrv.endswith("dew"):
        docstring_var = "Liquid" + docstring_var + "dew "
        docstring_eqn += "dew mole fractions"
    elif abbrv.endswith("bub"):
        docstring_var = "Vapor" + docstring_var + "bubble "
        docstring_eqn += "bubble mole fractions"
    else:
        _raise_dev_burnt_toast()

    if abbrv.startswith("t"):
        docstring_var += "temperature"
    elif abbrv.startswith("p"):
        docstring_var += "pressure"
    else:
        _raise_dev_burnt_toast()

    try:
        setattr(
            b,
            "log_mole_frac_" + abbrv,
            Var(
                b.params._pe_pairs,
                b.component_list,
                initialize=log(1 / len(b.component_list)),
                bounds=(None, 0),
                doc=docstring_var,
                units=None,
            ),
        )

        def rule_log_mole_frac(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = bub_dew_VL_comp_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif j not in set(vl_comps + henry_comps):
                return Constraint.Skip

            if abbrv.endswith("dew") and (l_only_comps != []):
                # Non-vaporisables present, no dew point
                return Constraint.Skip
            elif abbrv.endswith("bub") and (v_only_comps != []):
                # Non-condensables present, no bubble point
                return Constraint.Skip
            log_mole_frac = getattr(b, "log_mole_frac_" + abbrv)
            mole_frac = getattr(b, "_mole_frac_" + abbrv)

            return exp(log_mole_frac[p1, p2, j]) == mole_frac[p1, p2, j]

        setattr(
            b,
            "log_mole_frac_" + abbrv + "_eqn",
            Constraint(
                b.params._pe_pairs,
                b.component_list,
                rule=rule_log_mole_frac,
                doc=docstring_eqn,
            ),
        )
    except AttributeError:
        if hasattr(b, "log_mole_frac_" + abbrv):
            var = getattr(b, "log_mole_frac_" + abbrv)
            b.del_component(var)
        if hasattr(b, "log_mole_frac_" + abbrv + "_eqn"):
            eqn = getattr(b, "log_mole_frac_" + abbrv + "_eqn")
            b.del_component(eqn)
        raise
