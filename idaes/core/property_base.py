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
This module contains classes for property blocks and property parameter blocks.
"""

import sys

# Import Pyomo libraries
from pyomo.environ import Set, value, Var, Expression, Constraint
from pyomo.core.base.var import _VarData
from pyomo.core.base.expression import _ExpressionData
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.core.base.misc import tabular_writer

# Import IDAES cores
from idaes.core.process_block import ProcessBlock
from idaes.core import ProcessBlockData
from idaes.core import property_meta
from idaes.core import MaterialFlowBasis
from idaes.core.phases import Phase, PhaseData
from idaes.core.components import Component, ComponentData
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import (BurntToast,
                                        ConfigurationError,
                                        PropertyNotSupportedError,
                                        PropertyPackageError)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              number_variables,
                                              number_activated_constraints,
                                              number_activated_blocks)
from idaes.core.util import scaling as iscale
import idaes.logger as idaeslog

# Some more information about this module
__author__ = "Andrew Lee, John Eslick"

__all__ = ['StateBlockData',
           'StateBlock',
           'PhysicalParameterBlock']

# Set up logger
_log = idaeslog.getLogger(__name__)


class _lock_attribute_creation_context(object):
    """Context manager to lock creation of new attributes on a state block"""
    def __init__(self, block):
        self.block = block
    def __enter__(self):
        self.block._lock_attribute_creation = True
    def __exit__(self, exc_type, exc_value, traceback):
        self.block._lock_attribute_creation = False


class PhysicalParameterBlock(ProcessBlockData,
                             property_meta.HasPropertyClassMetadata):
    """
        This is the base class for thermophysical parameter blocks. These are
        blocks that contain a set of parameters associated with a specific
        thermophysical property package, and are linked to by all instances of
        that property package.
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("default_arguments", ConfigBlock(
            implicit=True,
            description="Default arguments to use with Property Package"))

    def build(self):
        """
        General build method for PropertyParameterBlocks. Inheriting models
        should call super().build.

        Args:
            None

        Returns:
            None
        """
        super(PhysicalParameterBlock, self).build()

        # Need this to work with the Helmholtz EoS package
        if not hasattr(self, "_state_block_class"):
            self._state_block_class = None

        # By default, property packages do not include inherent reactions
        self._has_inherent_reactions = False

        # This is a dict to store default property scaling factors. They are
        # defined in the parameter block to provide a universal default for
        # quantities in a particular kind of state block.  For example, you can
        # set flow scaling once instead of for every state block. Some of these
        # may be left for the user to set and some may be defined in a property
        # module where reasonable defaults can be defined a priori. See
        # set_default_scaling, get_default_scaling, and unset_default_scaling
        self.default_scaling_factor = {}

    def set_default_scaling(self, attrbute, value, index=None):
        """Set a default scaling factor for a property.

        Args:
            attribute: property attribute name
            value: default scaling factor
            index: for indexed properties, if this is not provied the scaling
                factor default applies to all indexed elements where specific
                indexes are no specifcally specified.

        Returns:
            None
        """
        self.default_scaling_factor[(attrbute, index)] = value

    def unset_default_scaling(self, attrbute, index=None):
        """Remove a previously set default value

        Args:
            attribute: property attribute name
            index: optional index for indexed properties

        Returns:
            None
        """
        try:
            del self.default_scaling_factor[(attrbute, index)]
        except KeyError:
            pass

    def get_default_scaling(self, attrbute, index=None):
        """ Returns a default scale factor for a property

        Args:
            attribute: property attribute name
            index: optional index for indexed properties

        Returns:
            None
        """
        try:
            # If a specific component data index exists
            return self.default_scaling_factor[(attrbute, index)]
        except KeyError:
            try:
                # indexed, but no specifc index?
                return self.default_scaling_factor[(attrbute, None)]
            except KeyError:
                # Can't find a default scale factor for what you asked for
                return None

    @property
    def state_block_class(self):
        if self._state_block_class is not None:
            return self._state_block_class
        else:
            raise AttributeError(
                "{} has not assigned a StateBlock class to be associated "
                "with this property package. Please contact the developer of "
                "the property package.".format(self.name))

    @state_block_class.setter
    def state_block_class(self, val):
        _log.warning("DEPRECATED: state_block_class should not be set "
                     "directly. Property package developers should set the "
                     "_state_block_class attribute instead.")
        self._state_block_class = val

    @property
    def has_inherent_reactions(self):
        return self._has_inherent_reactions

    def build_state_block(self, *args, **kwargs):
        """
        Methods to construct a StateBlock assoicated with this
        PhysicalParameterBlock. This will automatically set the parameters
        construction argument for the StateBlock.

        Returns:
            StateBlock

        """
        default = kwargs.pop("default", {})
        initialize = kwargs.pop("initialize", {})

        if initialize == {}:
            default["parameters"] = self
        else:
            for i in initialize.keys():
                initialize[i]["parameters"] = self

        return self.state_block_class(*args,
                                      **kwargs,
                                      default=default,
                                      initialize=initialize)

    def get_phase_component_set(self):
        """
        Method to get phase-component set for property package. If a phase-
        component set has not been constructed yet, this method will construct
        one.

        Args:
            None

        Returns:
            Phase-Component Set object
        """
        try:
            return self._phase_component_set
        except AttributeError:
            # Phase-component set does not exist, so create one.
            pc_set = []
            for p in self.phase_list:
                p_obj = self.get_phase(p)
                if p_obj.config.component_list is not None:
                    c_list = p_obj.config.component_list
                else:
                    c_list = self.component_list
                for j in c_list:
                    pc_set.append((p, j))

            self._phase_component_set = Set(initialize=pc_set, ordered=True)

            return self._phase_component_set

    def get_component(self, comp):
        """
        Method to retrieve a Component object based on a name from the
        component_list.

        Args:
            comp: name of Component object to retrieve

        Returns:
            Component object
        """
        obj = getattr(self, comp)
        if not isinstance(obj, ComponentData):
            raise PropertyPackageError(
                    "{} get_component found an attribute {}, but it does not "
                    "appear to be an instance of a Component object."
                    .format(self.name, comp))
        return obj

    def get_phase(self, phase):
        """
        Method to retrieve a Phase object based on a name from the phase_list.

        Args:
            phase: name of Phase object to retrieve

        Returns:
            Phase object
        """
        obj = getattr(self, phase)
        if not isinstance(obj, PhaseData):
            raise PropertyPackageError(
                    "{} get_phase found an attribute {}, but it does not "
                    "appear to be an instance of a Phase object."
                    .format(self.name, phase))
        return obj

    # TODO : Deprecate this code at some point
    def _validate_parameter_block(self):
        """
        Backwards compatability checks.

        This is code to check for old-style property packages and create
        the necessary Phase and Component objects.

        It also tries to catch some possible mistakes and provide the user with
        useful error messages.
        """
        try:
            # Check names in component list have matching Component objects
            for c in self.component_list:
                try:
                    obj = getattr(self, str(c))
                    if not isinstance(obj, ComponentData):
                        raise TypeError(
                                "Property package {} has an object {} whose "
                                "name appears in component_list but is not an "
                                "instance of Component".format(self.name, c))
                except AttributeError:
                    # No object with name c, must be old-style package
                    self._make_component_objects()
                    break
        except AttributeError:
            # No component list
            raise PropertyPackageError("Property package {} has not defined a "
                                       "component list.".format(self.name))

        try:
            # Valdiate that names in phase list have matching Phase objects
            for p in self.phase_list:
                try:
                    obj = getattr(self, str(p))
                    if not isinstance(obj, PhaseData):
                        raise TypeError(
                                "Property package {} has an object {} whose "
                                "name appears in phase_list but is not an "
                                "instance of Phase".format(self.name, p))
                except AttributeError:
                    # No object with name p, must be old-style package
                    self._make_phase_objects()
                    break
        except AttributeError:
            # No phase list
            raise PropertyPackageError("Property package {} has not defined a "
                                       "phase list.".format(self.name))

        # Also check that the phase-component set has been created.
        self.get_phase_component_set()

    def _make_component_objects(self):
        _log.warning("DEPRECATED: {} appears to be an old-style property "
                     "package. It will be automatically converted to a "
                     "new-style package, however users are strongly encouraged"
                     " to convert their property packages to use phase and "
                     "component objects."
                     .format(self.name))
        for c in self.component_list:
            if hasattr(self, c):
                # An object with this name already exists, raise exception
                raise PropertyPackageError(
                    "{} could not add Component object {} - an object with "
                    "that name already exists.".format(self.name, c))

            self.add_component(str(c), Component(
                default={"_component_list_exists": True}))

    def _make_phase_objects(self):
        _log.warning("DEPRECATED: {} appears to be an old-style property "
                     "package. It will be automatically converted to a "
                     "new-style package, however users are strongly encouraged"
                     " to convert their property packages to use phase and "
                     "component objects."
                     .format(self.name))
        for p in self.phase_list:
            if hasattr(self, p):
                # An object with this name already exists, raise exception
                raise PropertyPackageError(
                    "{} could not add Phase object {} - an object with "
                    "that name already exists.".format(self.name, p))

            try:
                pc_list = self.phase_comp[p]
            except AttributeError:
                pc_list = None
            self.add_component(str(p), Phase(
                default={"component_list": pc_list,
                         "_phase_list_exists": True}))


class StateBlock(ProcessBlock):
    """
        This is the base class for state block objects. These are used when
        constructing the SimpleBlock or IndexedBlock which will contain the
        PropertyData objects, and contains methods that can be applied to
        multiple StateBlockData objects simultaneously.
    """

    @property
    def component_list(self):
        return self._return_component_list()

    def _return_component_list(self):
        return self._get_parameter_block().component_list

    @property
    def phase_list(self):
        return self._return_phase_list()

    def _return_phase_list(self):
        return self._get_parameter_block().phase_list

    @property
    def phase_component_set(self):
        return self._return_phase_component_set()

    def _return_phase_component_set(self):
        return self._get_parameter_block().get_phase_component_set()

    @property
    def has_inherent_reactions(self):
        return self._has_inherent_reactions()

    def _has_inherent_reactions(self):
        return self._get_parameter_block().has_inherent_reactions

    # Need to separate the existence of inherent reactions from whether they
    # should be included in material balances
    # For some cases, Using an apparent species basis means they can be ignored
    @property
    def include_inherent_reactions(self):
        return self._include_inherent_reactions()

    def _include_inherent_reactions(self):
        return self._get_parameter_block().has_inherent_reactions

    def _get_parameter_block(self):
        try:
            return self._block_data_config_default["parameters"]
        except (KeyError, TypeError):
            # Need to get parameters from initialize dict
            # We will also confirm these are all the same whilst we are at it
            param = None
            if self._block_data_config_default is None:
                self._block_data_config_default = {}
            for v in self._block_data_config_initialize.values():
                if param is None:
                    param = v["parameters"]
                elif param is not v["parameters"]:
                    raise ConfigurationError(
                        "{} StateBlock must use the same parameter block for "
                        "elements. When using the initialize argument, please "
                        "ensure that the same value is used for all parameter "
                        "keys.".format(self.name))
            self._block_data_config_default["parameters"] = param
        return param

    def initialize(self, *args, **kwargs):
        """
        This is a default initialization routine for StateBlocks to ensure
        that a routine is present. All StateBlockData classes should
        overload this method with one suited to the particular property package

        Args:
            None

        Returns:
            None
        """
        raise NotImplementedError('{} property package has not implemented an'
                                  ' initialize method. Please contact '
                                  'the property package developer'
                                  .format(self.name))

    def report(self, index=(0), true_state=False,
               dof=False, ostream=None, prefix=""):
        """
        Default report method for StateBlocks. Returns a Block report populated
        with either the display or state variables defined in the
        StateBlockData class.

        Args:
            index : tuple of Block indices indicating which point in time (and
                    space if applicable) to report state at.
            true_state : whether to report the display variables (False
                    default) or the actual state variables (True)
            dof : whether to show local degrees of freedom in the report
                    (default=False)
            ostream : output stream to write report to
            prefix : string to append to the beginning of all output lines

        Returns:
            Printed output to ostream
        """

        if ostream is None:
            ostream = sys.stdout

        # Get DoF and model stats
        if dof:
            dof_stat = degrees_of_freedom(self[index])
            nv = number_variables(self[index])
            nc = number_activated_constraints(self[index])
            nb = number_activated_blocks(self[index])

        # Create stream table
        if true_state:
            disp_dict = self[index].define_state_vars()
        else:
            disp_dict = self[index].define_display_vars()

        stream_attributes = {}

        for k in disp_dict:
            for i in disp_dict[k]:
                if i is None:
                    stream_attributes[k] = disp_dict[k][i]
                else:
                    stream_attributes[k+" "+i] = disp_dict[k][i]

        # Write output
        max_str_length = 84
        tab = " "*4
        ostream.write("\n"+"="*max_str_length+"\n")

        lead_str = f"{prefix}State : {self.name}"
        trail_str = f"Index: {index}"
        mid_str = " "*(max_str_length-len(lead_str)-len(trail_str))
        ostream.write(lead_str+mid_str+trail_str)

        if dof:
            ostream.write("\n"+"="*max_str_length+"\n")
            ostream.write(f"{prefix}{tab}Local Degrees of Freedom: {dof_stat}")
            ostream.write('\n')
            ostream.write(f"{prefix}{tab}Total Variables: {nv}{tab}"
                          f"Activated Constraints: {nc}{tab}"
                          f"Activated Blocks: {nb}")

        ostream.write("\n"+"-"*max_str_length+"\n")
        ostream.write(f"{prefix}{tab}State Report")

        if any(isinstance(v, _VarData) for k, v in stream_attributes.items()):
            ostream.write("\n"*2)
            ostream.write(f"{prefix}{tab}Variables: \n\n")
            tabular_writer(
                    ostream,
                    prefix+tab,
                    ((k, v) for k, v in stream_attributes.items()
                        if isinstance(v, _VarData)),
                    ("Value", "Fixed", "Bounds"),
                    lambda k, v: ["{:#.5g}".format(value(v)),
                                  v.fixed,
                                  v.bounds])

        if any(isinstance(v, _ExpressionData) for
               k, v in stream_attributes.items()):
            ostream.write("\n"*2)
            ostream.write(f"{prefix}{tab}Expressions: \n\n")
            tabular_writer(
                    ostream,
                    prefix+tab,
                    ((k, v) for k, v in stream_attributes.items()
                        if isinstance(v, _ExpressionData)),
                    ("Value",),
                    lambda k, v: ["{:#.5g}".format(value(v))])

        ostream.write("\n"+"="*max_str_length+"\n")


class StateBlockData(ProcessBlockData):
    """
        This is the base class for state block data objects. These are
        blocks that contain the Pyomo components associated with calculating a
        set of thermophysical and transport properties for a given material.
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("parameters", ConfigValue(
            domain=is_physical_parameter_block,
            description="""A reference to an instance of the Property Parameter
Block associated with this property package."""))
    CONFIG.declare("defined_state", ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Flag indicating if incoming state is fully defined",
            doc="""Flag indicating whether the state should be considered fully
defined, and thus whether constraints such as sum of mass/mole fractions should
be included,
**default** - False.
**Valid values:** {
**True** - state variables will be fully defined,
**False** - state variables will not be fully defined.}"""))
    CONFIG.declare("has_phase_equilibrium", ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Phase equilibrium constraint flag",
            doc="""Flag indicating whether phase equilibrium constraints
should be constructed in this state block,
**default** - True.
**Valid values:** {
**True** - StateBlock should calculate phase equilibrium,
**False** - StateBlock should not calculate phase equilibrium.}"""))

    def __init__(self, *args, **kwargs):
        self._lock_attribute_creation = False
        super().__init__(*args, **kwargs)

    def lock_attribute_creation_context(self):
        """Returns a context manager that does not allow attributes to be created
        while in the context and allows attributes to be created normally outside
        the context.
        """
        return _lock_attribute_creation_context(self)

    def is_property_constructed(self, attr):
        """Returns True if the attribute ``attr`` already exists, or false if it
        would be added in ``__getattr__``, or does not exist.

        Args:
            attr (str): Attribute name to check

        Return:
            True if the attribute is already constructed, False otherwise
        """
        with self.lock_attribute_creation_context():
            return hasattr(self, attr)

    @property
    def component_list(self):
        return self.parent_component()._return_component_list()

    @property
    def phase_list(self):
        return self.parent_component()._return_phase_list()

    @property
    def phase_component_set(self):
        return self.parent_component()._return_phase_component_set()

    @property
    def has_inherent_reactions(self):
        return self.parent_component()._has_inherent_reactions()

    @property
    def include_inherent_reactions(self):
        return self.parent_component()._include_inherent_reactions()

    def build(self):
        """
        General build method for StateBlockDatas.

        Args:
            None

        Returns:
            None
        """
        super(StateBlockData, self).build()
        add_object_reference(self, "_params", self.config.parameters)

        # TODO: Deprecate this at some point
        # Backwards compatability check for old-style property packages
        self._params._validate_parameter_block()

    @property
    def params(self):
        return self._params

    def define_state_vars(self):
        """
        Method that returns a dictionary of state variables used in property
        package. Implement a placeholder method which returns an Exception to
        force users to overload this.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' define_state_vars method. Please contact '
                                  'the property package developer.')

    def define_port_members(self):
        """
        Method used to specify components to populate Ports with. Defaults to
        define_state_vars, and developers should overload as required.
        """
        return self.define_state_vars()

    def define_display_vars(self):
        """
        Method used to specify components to use to generate stream tables and
        other outputs. Defaults to define_state_vars, and developers should
        overload as required.
        """
        return self.define_state_vars()

    def get_material_flow_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for material flow to use in
        the material balances.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' get_material_flow_terms method. Please '
                                  'contact the property package developer.')

    def get_material_density_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for material density to use in
        the material balances .
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' get_material_density_terms method. Please '
                                  'contact the property package developer.')

    def get_material_diffusion_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for material diffusion to use
        in the material balances.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' get_material_diffusion_terms method. '
                                  'Please contact the property package '
                                  'developer.')

    def get_enthalpy_flow_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for enthalpy flow to use in
        the energy balances.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' get_enthalpy_flow_terms method. Please '
                                  'contact the property package developer.')

    def get_energy_density_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for enthalpy density to use in
        the energy balances.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' get_energy_density_terms method. Please '
                                  'contact the property package developer.')

    def get_energy_diffusion_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for energy diffusion to use in
        the energy balances.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' get_energy_diffusion_terms method. '
                                  'Please contact the property package '
                                  'developer.')

    def get_material_flow_basis(self, *args, **kwargs):
        """
        Method which returns an Enum indicating the basis of the material flow
        term.
        """
        return MaterialFlowBasis.other

    def calculate_bubble_point_temperature(self, *args, **kwargs):
        """
        Method which computes the bubble point temperature for a multi-
        component mixture given a pressure and mole fraction.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' calculate_bubble_point_temperature method.'
                                  ' Please contact the property package '
                                  'developer.')

    def calculate_dew_point_temperature(self, *args, **kwargs):
        """
        Method which computes the dew point temperature for a multi-
        component mixture given a pressure and mole fraction.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' calculate_dew_point_temperature method.'
                                  ' Please contact the property package '
                                  'developer.')

    def calculate_bubble_point_pressure(self, *args, **kwargs):
        """
        Method which computes the bubble point pressure for a multi-
        component mixture given a temperature and mole fraction.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' calculate_bubble_point_pressure method.'
                                  ' Please contact the property package '
                                  'developer.')

    def calculate_dew_point_pressure(self, *args, **kwargs):
        """
        Method which computes the dew point pressure for a multi-
        component mixture given a temperature and mole fraction.
        """
        raise NotImplementedError('{} property package has not implemented the'
                                  ' calculate_dew_point_pressure method.'
                                  ' Please contact the property package '
                                  'developer.')

    def __getattr__(self, attr):
        """
        This method is used to avoid generating unnecessary property
        calculations in state blocks. __getattr__ is called whenever a
        property is called for, and if a propery does not exist, it looks for
        a method to create the required property, and any associated
        components.

        Create a property calculation if needed. Return an attrbute error if
        attr == 'domain' or starts with a _ . The error for _ prevents a
        recursion error if trying to get a function to create a property and
        that function doesn't exist.  Pyomo also ocasionally looks for things
        that start with _ and may not exist.  Pyomo also looks for the domain
        attribute, and it may not exist.
        This works by creating a property calculation by calling the "_"+attr
        function.

        A list of __getattr__ calls is maintained in self.__getattrcalls to
        check for recursive loops which maybe useful for debugging. This list
        is cleared after __getattr__ completes successfully.

        Args:
            attr: an attribute to create and return. Should be a property
                  component.
        """
        if self._lock_attribute_creation:
            raise AttributeError(
                f"{attr} does not exist, and attribute creation is locked.")

        def clear_call_list(self, attr):
            """Local method for cleaning up call list when a call is handled.

                Args:
                    attr: attribute currently being handled
            """
            if self.__getattrcalls[-1] == attr:
                if len(self.__getattrcalls) <= 1:
                    del self.__getattrcalls
                else:
                    del self.__getattrcalls[-1]
            else:
                raise PropertyPackageError(
                        "{} Trying to remove call {} from __getattr__"
                        " call list, however this is not the most "
                        "recent call in the list ({}). This indicates"
                        " a bug in the __getattr__ calls. Please "
                        "contact the IDAES developers with this bug."
                        .format(self.name, attr, self.__getattrcalls[-1]))

        # Check that attr is not something we shouldn't touch
        if attr == "domain" or attr.startswith("_"):
            # Don't interfere with anything by getting attributes that are
            # none of my business
            raise PropertyPackageError(
                    '{} {} does not exist, but is a protected '
                    'attribute. Check the naming of your '
                    'components to avoid any reserved names'
                    .format(self.name, attr))

        if attr == "config":
            try:
                self._get_config_args()
                return self.config
            except:
                raise BurntToast("{} getattr method was triggered by a call "
                                 "to the config block, but _get_config_args "
                                 "failed. This should never happen.")

        # Check for recursive calls
        try:
            # Check if __getattrcalls is initialized
            self.__getattrcalls
        except AttributeError:
            # Initialize it
            self.__getattrcalls = [attr]
        else:
            # Check to see if attr already appears in call list
            if attr in self.__getattrcalls:
                # If it does, indicates a recursive loop.
                if attr == self.__getattrcalls[-1]:
                    # attr method is calling itself
                    self.__getattrcalls.append(attr)
                    raise PropertyPackageError(
                                    '{} _{} made a recursive call to '
                                    'itself, indicating a potential '
                                    'recursive loop. This is generally '
                                    'caused by the {} method failing to '
                                    'create the {} component.'
                                    .format(self.name, attr, attr, attr))
                else:
                    self.__getattrcalls.append(attr)
                    raise PropertyPackageError(
                                    '{} a potential recursive loop has been '
                                    'detected whilst trying to construct {}. '
                                    'A method was called, but resulted in a '
                                    'subsequent call to itself, indicating a '
                                    'recursive loop. This may be caused by a '
                                    'method trying to access a component out '
                                    'of order for some reason (e.g. it is '
                                    'declared later in the same method). See '
                                    'the __getattrcalls object for a list of '
                                    'components called in the __getattr__ '
                                    'sequence.'
                                    .format(self.name, attr))
            # If not, add call to list
            self.__getattrcalls.append(attr)

        # Get property information from properties metadata
        try:
            m = self.config.parameters.get_metadata().properties

            if m is None:
                raise PropertyPackageError(
                        '{} property package get_metadata()'
                        ' method returned None when trying to create '
                        '{}. Please contact the developer of the '
                        'property package'.format(self.name, attr))
        except KeyError:
            # If attr not in metadata, assume package does not
            # support property
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                    '{} {} is not supported by property package (property is '
                    'not listed in package metadata properties).'
                    .format(self.name, attr, attr))

        # Get method name from resulting properties
        try:
            if m[attr]['method'] is None:
                # If method is none, property should be constructed
                # by property package, so raise PropertyPackageError
                clear_call_list(self, attr)
                raise PropertyPackageError(
                        '{} {} should be constructed automatically '
                        'by property package, but is not present. '
                        'This can be caused by methods being called '
                        'out of order.'.format(self.name, attr))
            elif m[attr]['method'] is False:
                # If method is False, package does not support property
                # Raise NotImplementedError
                clear_call_list(self, attr)
                raise PropertyNotSupportedError(
                        '{} {} is not supported by property package '
                        '(property method is listed as False in '
                        'package property metadata).'
                        .format(self.name, attr))
            elif isinstance(m[attr]['method'], str):
                # Try to get method name in from PropertyBlock object
                try:
                    f = getattr(self, m[attr]['method'])
                except AttributeError:
                    # If fails, method does not exist
                    clear_call_list(self, attr)
                    raise PropertyPackageError(
                            '{} {} package property metadata method '
                            'returned a name that does not correspond'
                            ' to any method in the property package. '
                            'Please contact the developer of the '
                            'property package.'.format(self.name, attr))
            else:
                # Otherwise method name is invalid
                clear_call_list(self, attr)
                raise PropertyPackageError(
                             '{} {} package property metadata method '
                             'returned invalid value for method name. '
                             'Please contact the developer of the '
                             'property package.'
                             .format(self.name, attr))
        except KeyError:
            # No method key - raise Exception
            # Need to use an AttributeError so Pyomo.DAE will handle this
            clear_call_list(self, attr)
            raise PropertyNotSupportedError(
                    '{} package property metadata method '
                    'does not contain a method for {}. '
                    'Please select a package which supports '
                    'the necessary properties for your process.'
                    .format(self.name, attr))

        # Call attribute if it is callable
        # If this fails, it should return a meaningful error.
        if callable(f):
            try:
                f()
            except Exception:
                # Clear call list and reraise error
                clear_call_list(self, attr)
                raise
        else:
            # If f is not callable, inform the user and clear call list
            clear_call_list(self, attr)
            raise PropertyPackageError(
                    '{} tried calling attribute {} in order to create '
                    'component {}. However the method is not callable.'
                    .format(self.name, f, attr))

        # Clear call list, and return
        comp = getattr(self, attr)
        clear_call_list(self, attr)
        return comp

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # Get scaling factor defaults, if no scaling factor set
        for v in self.component_data_objects(
            (Constraint, Var, Expression),
            descend_into=False):
            if iscale.get_scaling_factor(v) is None: # don't replace if set
                name = v.getname().split("[")[0]
                index = v.index()
                sf = self.config.parameters.get_default_scaling(name, index)
                if sf is not None:
                    iscale.set_scaling_factor(v, sf)
