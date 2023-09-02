#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains classes for property blocks and property parameter blocks.
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

import sys

# Import Pyomo libraries
from pyomo.environ import Set, value, Var, Expression, Constraint, Reference
from pyomo.core.base.var import _VarData
from pyomo.core.base.expression import _ExpressionData
from pyomo.common.config import ConfigBlock, ConfigValue, Bool
from pyomo.common.formatting import tabular_writer
from pyomo.network import Port

# Import IDAES cores
from idaes.core.base.process_block import ProcessBlock
from idaes.core import ProcessBlockData, MaterialFlowBasis
from idaes.core.base import property_meta
from idaes.core.base.phases import PhaseData
from idaes.core.base.components import ComponentData
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import (
    ConfigurationError,
    PropertyPackageError,
)
from idaes.core.util.misc import add_object_reference
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_activated_constraints,
    number_activated_blocks,
)
from idaes.core.util import scaling as iscale
import idaes.logger as idaeslog
from idaes.core.base.util import build_on_demand
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
)


# Some more information about this module
__author__ = "Andrew Lee, John Eslick"

__all__ = ["StateBlockData", "StateBlock", "PhysicalParameterBlock"]

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


class PhysicalParameterBlock(ProcessBlockData, property_meta.HasPropertyClassMetadata):
    """
    This is the base class for thermophysical parameter blocks. These are
    blocks that contain a set of parameters associated with a specific
    thermophysical property package, and are linked to by all instances of
    that property package.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "default_arguments",
        ConfigBlock(
            implicit=True, description="Default arguments to use with Property Package"
        ),
    )

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

    @property
    def state_block_class(self):
        if self._state_block_class is not None:
            return self._state_block_class
        else:
            raise AttributeError(
                "{} has not assigned a StateBlock class to be associated "
                "with this property package. Please contact the developer of "
                "the property package.".format(self.name)
            )

    @property
    def has_inherent_reactions(self):
        return self._has_inherent_reactions

    def build_state_block(self, *args, **kwargs):
        """
        Methods to construct a StateBlock associated with this
        PhysicalParameterBlock. This will automatically set the parameters
        construction argument for the StateBlock.

        Returns:
            StateBlock

        """
        # default = kwargs.pop("default", {})
        initialize = kwargs.pop("initialize", {})

        # TODO: Should find a better way to do this
        # For now, trigger get_phase_component_set to make sure
        # it has been constructed before building the state block
        self.get_phase_component_set()

        if initialize == {}:
            kwargs["parameters"] = self
        else:
            for i in initialize.keys():
                initialize[i]["parameters"] = self

        return self.state_block_class(  # pylint: disable=not-callable
            *args, **kwargs, initialize=initialize
        )

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
                "appear to be an instance of a Component object.".format(
                    self.name, comp
                )
            )
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
                "appear to be an instance of a Phase object.".format(self.name, phase)
            )
        return obj


class StateBlock(ProcessBlock):
    """
    This is the base class for state block objects. These are used when
    constructing the SimpleBlock or IndexedBlock which will contain the
    PropertyData objects, and contains methods that can be applied to
    multiple StateBlockData objects simultaneously.
    """

    # Set default initializer
    default_initializer = BlockTriangularizationInitializer

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
        # PYLINT-WHY: self._block_data_config_default is set elsewhere,
        # and is supposed to already exist as an attribute by the time
        # (or, it will raise AttributeError at L410)
        # this method is called, so the pylint errors here are false positives
        # pylint: disable=no-member,access-member-before-definition
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
                        "keys.".format(self.name)
                    )
            self._block_data_config_default["parameters"] = param
        return param

    @property
    def params(self):
        return self._get_parameter_block()

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        raise NotImplementedError

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
        raise NotImplementedError(
            "{} property package has not implemented an"
            " initialize method. Please contact "
            "the property package developer".format(self.name)
        )

    def report(self, index=(0), true_state=False, dof=False, ostream=None, prefix=""):
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
                    stream_attributes[k + " " + i] = disp_dict[k][i]

        # Write output
        max_str_length = 84
        tab = " " * 4
        ostream.write("\n" + "=" * max_str_length + "\n")

        lead_str = f"{prefix}State : {self.name}"
        trail_str = f"Index: {index}"
        mid_str = " " * (max_str_length - len(lead_str) - len(trail_str))
        ostream.write(lead_str + mid_str + trail_str)

        if dof:
            ostream.write("\n" + "=" * max_str_length + "\n")
            ostream.write(f"{prefix}{tab}Local Degrees of Freedom: {dof_stat}")
            ostream.write("\n")
            ostream.write(
                f"{prefix}{tab}Total Variables: {nv}{tab}"
                f"Activated Constraints: {nc}{tab}"
                f"Activated Blocks: {nb}"
            )

        ostream.write("\n" + "-" * max_str_length + "\n")
        ostream.write(f"{prefix}{tab}State Report")

        if any(isinstance(v, _VarData) for k, v in stream_attributes.items()):
            ostream.write("\n" * 2)
            ostream.write(f"{prefix}{tab}Variables: \n\n")
            tabular_writer(
                ostream,
                prefix + tab,
                (
                    (k, v)
                    for k, v in stream_attributes.items()
                    if isinstance(v, _VarData)
                ),
                ("Value", "Fixed", "Bounds"),
                lambda k, v: ["{:#.5g}".format(value(v)), v.fixed, v.bounds],
            )

        if any(isinstance(v, _ExpressionData) for k, v in stream_attributes.items()):
            ostream.write("\n" * 2)
            ostream.write(f"{prefix}{tab}Expressions: \n\n")
            tabular_writer(
                ostream,
                prefix + tab,
                (
                    (k, v)
                    for k, v in stream_attributes.items()
                    if isinstance(v, _ExpressionData)
                ),
                ("Value",),
                lambda k, v: ["{:#.5g}".format(value(v))],
            )

        ostream.write("\n" + "=" * max_str_length + "\n")

    def get_port_reference_name(self, component_name, port_name):
        """
        Get the standard name of a "port reference", the component
        accessed on the port by "port_name.component_name".

        Args:
            component_name - name of port member of interest (str)
            port_name - name of Port object (str)

        Returns:
            str with name for Reference used for Port member
        """
        return f"_{component_name}_{port_name}_ref"

    def build_port(
        self,
        doc=None,
        slice_index=None,
        index=None,
    ):
        """
        Constructs a Port based on this StateBlock attached to the target block.

        Args:
            doc - doc string or Prot object
            slice_index - Slice index (e.g. (slice(None), 0.0) that will be
                used to index self when constructing port references. Default = None.
            index - time index to use when calling define_port_members. Default = None.

        Returns:
            Port object and list of tuples with form (Reference, member name)
        """
        if slice_index is None:
            slice_index = Ellipsis
        if index is None:
            index = self.index_set().first()

        # Create empty Port
        p = Port(doc=doc)

        # Get dict of Port members and names
        # Need to get a representative member of StateBlockDatas
        port_member_dict = self[index].define_port_members()

        # Create References for port members
        ref_name_list = []
        for name, component in port_member_dict.items():
            if not component.is_indexed():
                slicer = self[slice_index].component(component.local_name)
            else:
                slicer = self[slice_index].component(component.local_name)[...]

            ref = Reference(slicer)
            ref_name_list.append((ref, name))

            # Add Reference to Port
            p.add(ref, name)

        return p, ref_name_list


class StateBlockData(ProcessBlockData):
    """
    This is the base class for state block data objects. These are
    blocks that contain the Pyomo components associated with calculating a
    set of thermophysical and transport properties for a given material.
    """

    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "parameters",
        ConfigValue(
            domain=is_physical_parameter_block,
            description="""A reference to an instance of the Property Parameter
Block associated with this property package.""",
        ),
    )
    CONFIG.declare(
        "defined_state",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Flag indicating if incoming state is fully defined",
            doc="""Flag indicating whether the state should be considered fully
defined, and thus whether constraints such as sum of mass/mole fractions should
be included,
**default** - False.
**Valid values:** {
**True** - state variables will be fully defined,
**False** - state variables will not be fully defined.}""",
        ),
    )
    CONFIG.declare(
        "has_phase_equilibrium",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Phase equilibrium constraint flag",
            doc="""Flag indicating whether phase equilibrium constraints
should be constructed in this state block,
**default** - True.
**Valid values:** {
**True** - StateBlock should calculate phase equilibrium,
**False** - StateBlock should not calculate phase equilibrium.}""",
        ),
    )

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
        # TODO: Should refactor parent so this is not private
        # pylint: disable-next=protected-access
        return self.parent_component()._return_component_list()

    @property
    def phase_list(self):
        # TODO: Should refactor parent so this is not private
        # pylint: disable-next=protected-access
        return self.parent_component()._return_phase_list()

    @property
    def phase_component_set(self):
        # TODO: Should refactor parent so this is not private
        # pylint: disable-next=protected-access
        return self.parent_component()._return_phase_component_set()

    @property
    def has_inherent_reactions(self):
        # TODO: Should refactor parent so this is not private
        # pylint: disable-next=protected-access
        return self.parent_component()._has_inherent_reactions()

    @property
    def include_inherent_reactions(self):
        # TODO: Should refactor parent so this is not private
        # pylint: disable-next=protected-access
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

    @property
    def params(self):
        return self._params

    def define_state_vars(self):
        """
        Method that returns a dictionary of state variables used in property
        package. Implement a placeholder method which returns an Exception to
        force users to overload this.
        """
        raise NotImplementedError(
            f"{self.name} property package has not implemented the"
            f" define_state_vars method. Please contact "
            f"the property package developer."
        )

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
        raise NotImplementedError(
            "{} property package has not implemented the"
            " get_material_flow_terms method. Please "
            "contact the property package developer."
        )

    def get_material_density_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for material density to use in
        the material balances .
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " get_material_density_terms method. Please "
            "contact the property package developer."
        )

    def get_material_diffusion_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for material diffusion to use
        in the material balances.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " get_material_diffusion_terms method. "
            "Please contact the property package "
            "developer."
        )

    def get_enthalpy_flow_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for enthalpy flow to use in
        the energy balances.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " get_enthalpy_flow_terms method. Please "
            "contact the property package developer."
        )

    def get_energy_density_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for enthalpy density to use in
        the energy balances.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " get_energy_density_terms method. Please "
            "contact the property package developer."
        )

    def get_energy_diffusion_terms(self, *args, **kwargs):
        """
        Method which returns a valid expression for energy diffusion to use in
        the energy balances.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " get_energy_diffusion_terms method. "
            "Please contact the property package "
            "developer."
        )

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
        raise NotImplementedError(
            "{} property package has not implemented the"
            " calculate_bubble_point_temperature method."
            " Please contact the property package "
            "developer."
        )

    def calculate_dew_point_temperature(self, *args, **kwargs):
        """
        Method which computes the dew point temperature for a multi-
        component mixture given a pressure and mole fraction.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " calculate_dew_point_temperature method."
            " Please contact the property package "
            "developer."
        )

    def calculate_bubble_point_pressure(self, *args, **kwargs):
        """
        Method which computes the bubble point pressure for a multi-
        component mixture given a temperature and mole fraction.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " calculate_bubble_point_pressure method."
            " Please contact the property package "
            "developer."
        )

    def calculate_dew_point_pressure(self, *args, **kwargs):
        """
        Method which computes the dew point pressure for a multi-
        component mixture given a temperature and mole fraction.
        """
        raise NotImplementedError(
            "{} property package has not implemented the"
            " calculate_dew_point_pressure method."
            " Please contact the property package "
            "developer."
        )

    def __getattr__(self, attr):
        """
        This method is used to avoid generating unnecessary property
        calculations in state blocks. __getattr__ is called whenever a
        property is called for, and if a property does not exist, it looks for
        a method to create the required property, and any associated
        components.

        Create a property calculation if needed. Return an attribute error if
        attr == 'domain' or starts with a _ . The error for _ prevents a
        recursion error if trying to get a function to create a property and
        that function doesn't exist.  Pyomo also occasionally looks for things
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
        try:
            # try the Pyomo Block's __getattr__ method first which will return
            # decorators for creating components on the block (e.g. Expression,
            # Constraint, ...).
            return super().__getattr__(attr)
        except AttributeError:
            return build_on_demand(self, attr)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        # Get scaling factor defaults, if no scaling factor set
        for v in self.component_data_objects(
            (Constraint, Var, Expression), descend_into=False
        ):
            if iscale.get_scaling_factor(v) is None:  # don't replace if set
                name = v.getname().split("[")[0]
                index = v.index()
                sf = self.config.parameters.get_default_scaling(name, index)
                if sf is not None:
                    iscale.set_scaling_factor(v, sf)
