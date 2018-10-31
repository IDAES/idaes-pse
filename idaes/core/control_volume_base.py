##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Base class for control volumes
"""

from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyutilib.enum import Enum

# Import IDAES cores
from idaes.core import ProcessBlockData, useDefault
from idaes.core.util.config import (is_property_parameter_block,
                                    is_reaction_parameter_block)
from idaes.core.util.exceptions import (ConfigurationError,
                                        DynamicError,
                                        IDAESError,
                                        PropertyPackageError)

__author__ = "Andrew Lee"

# Set up logger
logger = logging.getLogger(__name__)


# Enumerate options for material balances
MaterialBalanceType = Enum(
    'none',
    'componentPhase',
    'componentTotal',
    'elementTotal',
    'total')

# Enumerate options for energy balances
EnergyBalanceType = Enum(
    'none',
    'enthalpyPhase',
    'enthalpyTotal',
    'energyPhase',
    'energyTotal')

# Enumerate options for momentum balances
MomentumBalanceType = Enum(
    'none',
    'pressureTotal',
    'pressurePhase',
    'momentumTotal',
    'momentumPhase')

# Enumerate options for flow direction
FlowDirection = Enum(
    'forward',
    'backward')

# Set up example ConfigBlock that will work with ControlVolume autobuild method
CONFIG_Base = ProcessBlockData.CONFIG()
CONFIG_Base.declare("dynamic", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic,
**default** - useDefault.  **Valid values:** {
**useDefault** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}"""))
CONFIG_Base.declare("include_holdup", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True, **default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
CONFIG_Base.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.componentPhase,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase. **Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase - use phase component balances,
**MaterialBalanceType.componentTotal - use total component balances,
**MaterialBalanceType.elementTotal - use total element balances,
**MaterialBalanceType.total - use total material balance.}"""))
CONFIG_Base.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyPhase,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyPhase. **Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
CONFIG_Base.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal. **Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
CONFIG_Base.declare("has_rate_reactions", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Rate reaction construction flag",
        doc="""Indicates whether terms for rate controlled reactions should be
constructed, **default** - False. **Valid values:** {
**True** - include kinetic reaction terms,
**False** - exclude kinetic reaction terms.}"""))
CONFIG_Base.declare("has_equilibrium_reactions", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Equilibrium reaction construction flag",
        doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed, **default** - False. **Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}"""))
CONFIG_Base.declare("has_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Phase equilibrium construction flag",
        doc="""Indicates whether terms for phase equilibrium should be
constructed, **default** = False. **Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}"""))
CONFIG_Base.declare("has_mass_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Mass transfer term construction flag",
        doc="""Indicates whether terms for mass transfer should be constructed,
**default** - False. **Valid values:** {
**True** - include mass transfer terms,
**False** - exclude mass transfer terms.}"""))
CONFIG_Base.declare("has_heat_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False. **Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}"""))
CONFIG_Base.declare("has_work_transfer", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Work transfer term construction flag",
        doc="""Indicates whether terms for work transfer should be constructed,
**default** - False. **Valid values** {,
**True** - include work transfer terms,
**False** - exclude work transfer terms.}"""))
CONFIG_Base.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed, **default** - False. **Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
CONFIG_Base.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_property_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault. **Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
CONFIG_Base.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
 and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}"""))
CONFIG_Base.declare("reaction_package", ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for control volume",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None. **Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}"""))
CONFIG_Base.declare("reaction_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
 and used when constructing these, **default** - None. **Valid values:** {
see reaction package for documentation.}"""))


class ControlVolumeBase(ProcessBlockData):
    """
    The ControlVolumeBase Class forms the base class for all IDAES
    ControlVolume models. The purpose of this class is to automate the tasks
    common to all control volume blockss and ensure that the necessary
    attributes of a control volume block are present.

    The most signfiicant role of the ControlVolumeBase class is to set up the
    bconstruction arguments for the control volume block, automatically link to
    the time domain of the parent block, and to get the information about the
    property and reaction packages.
    """

    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("dynamic", ConfigValue(
        domain=In([useDefault, True, False]),
        default=useDefault,
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic, **default** -
useDefault. **Valid values:** {
**useDefault** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_property_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault. **Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
 and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}"""))
    CONFIG.declare("reaction_package", ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for control volume",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None. **Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}"""))
    CONFIG.declare("reaction_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
 and used when constructing these, **default** - None. **Valid values:** {
see reaction package for documentation.}"""))
    CONFIG.declare("auto_construct", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Argument indicating whether ControlVolume should "
                    "automatically construct balance equations",
        doc="""If set to True, this argument will trigger the auto_construct
method which will attempt to construct a set of material, energy and momentum
balance equations based on the parent unit's config block. THe parent unit must
have a config block which derives from CONFIG_Base, **default** - False.
**Valid values:** {**True** - use automatic construction,
**False** - do not use automatic construciton.}"""))

    def build(self):
        """
        General build method for Control Volumes blocks. This method calls a
        number of sub-methods which automate the construction of expected
        attributes of all ControlVolume blocks.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        # Setup dynamics flag and time domain
        self._setup_dynamics()

        # Get property package details
        self._get_property_package()

        # Get indexing sets
        self._get_indexing_sets()

        # Get reaction package details (as necessary)
        self._get_reaction_package()

        if self.config.auto_construct is True:
            self._auto_construct()

    def _auto_construct(self):
        """
        Placeholder _auto_construct method to ensure a useful exception is
        returned if auto_build is set to True but something breaks in the
        process. Derived ControlVolume classes should overload this.

        Args:
            None

        Returns:
            None
        """
        raise IDAESError("{} auto-construct failed as ControlVolume "
                         "class failed to create _auto_construct method."
                         "Please contact the IDAES developers with this bug."
                         .format(self.name))

    def _setup_dynamics(self):
        """
        This method automates the setting of the dynamic flag and time domain
        for control volume blocks.

        If dynamic flag is 'use_parent_value', method attempts to get the value
        of the dynamic flag from the parent model, otherwise the local value is
        used. The time domain is always collected from the parent model.

        Finally, if dynamic = True, the include_holdup flag is checked to
        ensure it is also True.

        Args:
            None

        Returns:
            None
        """
        # Check the dynamic flag, and retrieve if necessary
        if self.config.dynamic == useDefault:
            # Get dynamic flag from parent
            try:
                self.config.dynamic = self.parent_block().config.dynamic
            except AttributeError:
                # If parent does not have dynamic flag, raise Exception
                raise DynamicError('{} has a parent model '
                                   'with no dynamic attribute.'
                                   .format(self.name))

        # Try to get reference to time object from parent
        try:
            object.__setattr__(self, "time", self.parent_block().time)
        except AttributeError:
            raise DynamicError('{} has a parent model '
                               'with no time domain'.format(self.name))

    def _get_property_package(self):
        """
        This method gathers the necessary information about the property
        package to be used in the control volume block.

        If a property package has not been provided by the user, the method
        searches up the model tree until it finds an object with the
        'default_property_package' attribute and uses this package for the
        control volume block.

        The method also gathers any default construction arguments specified
        for the property package and combines these with any arguments
        specified by the user for the control volume block (user specified
        arguments take priority over defaults).

        Args:
            None

        Returns:
            None
        """
        # Get property_package block if not provided in arguments
        parent = self.parent_block()
        if self.config.property_package == useDefault:
            # Try to get property_package from parent
            try:
                if parent.config.property_package is None:
                    parent.config.property_package = \
                        self._get_default_prop_pack()

                self.config.property_package = parent.config.property_package
            except AttributeError:
                self.config.property_package = self._get_default_prop_pack()

        # Get module of property package
        self.property_module = self.config.property_package.property_module

        # Check for any flowsheet level build arguments
        for k in self.config.property_package.config.default_arguments:
            if k not in self.config.property_package_args:
                self.config.property_package_args[k] = \
                    self.config.property_package.config.default_arguments[k]

    def _get_default_prop_pack(self):
        """
        This method is used to find a default property package defined at the
        flowsheet level if a package is not provided as an argument when
        instantiating the control volume block.

        Args:
            None

        Returns:
            None
        """
        parent = self.parent_block()
        while True:
            if hasattr(parent.config, "default_property_package"):
                break
            else:
                if parent.parent_block() is None:
                    raise ConfigurationError(
                            '{} no property package provided and '
                            'no default defined. Found end of '
                            'parent tree.'.format(self.name))
                elif parent.parent_block() == parent:
                    raise ConfigurationError(
                            '{} no property package provided and '
                            'no default defined. Found recursive '
                            'loop in parent tree.'.format(self.name))
                parent = parent.parent_block()

        logger.info('{} Using default property package'
                    .format(self.name))

        if parent.config.default_property_package is None:
            raise ConfigurationError(
                             '{} no default property package has been '
                             'specified at flowsheet level ('
                             'default_property_package = None)'
                             .format(self.name))

        return parent.config.default_property_package

    def _get_indexing_sets(self):
        """
        This method collects all necessary indexing sets from property
        parameter block and makes references to these for use within the
        control volume block. Collected indexing sets are phase_list and
        component_list.

        Args:
            None

        Returns:
            None
        """
        # Get phase and component list(s)
        try:
            # TODO : Look at ways to use Pyomo references, or create new Set
            object.__setattr__(self, "phase_list",
                               self.config.property_package.phase_list)
        except AttributeError:
            raise PropertyPackageError(
                    '{} property_package provided does not '
                    'contain a phase_list. '
                    'Please contact the developer of the property package.'
                    .format(self.name))
        try:
            # TODO : Look at ways to use Pyomo references, or create new Set
            object.__setattr__(self, "component_list",
                               self.config.property_package.component_list)
        except AttributeError:
            raise PropertyPackageError(
                    '{} property_package provided does not '
                    'contain a component_list. '
                    'Please contact the developer of the property package.'
                    .format(self.name))

    def _get_reaction_package(self):
        """
        This method gathers the necessary information about the reaction
        package to be used in the control volume block (if required).

        If a reaction package has been provided by the user, the method
        gathers any default construction arguments specified
        for the reaction package and combines these with any arguments
        specified by the user for the control volume block (user specified
        arguments take priority over defaults).

        Args:
            None

        Returns:
            None
        """
        if self.config.reaction_package is not None:
            # Get module of reaction package
            self.reaction_module = self.config.reaction_package.property_module

            # Check for any flowsheet level build arguments
            for k in self.config.reaction_package.config.default_arguments:
                if k not in self.config.reaction_package_args:
                    self.config.reaction_package_args[k] = \
                       self.config.reaction_package.config.default_arguments[k]

    def _validate_add_balance_arguments(self, dynamic, has_holdup):
        """
        Method to validate dynamic and has_holdup arguments used by many
        balance equation methods.

        Args:
            dynamic, has_holdup

        Returns:
            Validated values of dynamic and has_holdup
        """
        # If dynamic argument not provided, try to get argument from parent
        if dynamic == useDefault:
            dynamic = self.config.dynamic
        elif dynamic and not self.config.dynamic:
            raise DynamicError("{} cannot have dynamic balance equations "
                               "within a steady-state control volume."
                               .format(self.name))

        # If dynamic = True, has_holdup must also be True
        if dynamic and not has_holdup:
            raise ConfigurationError(
                    "{} invalid arguments for dynamic and has_holdup. "
                    "If dynamic = True, has_holdup must also be True (was "
                    "False)".format(self.name))

        return dynamic, has_holdup

    # Add placeholder methods for all types of material, energy and momentum
    # balance equations which return NotImplementedErrors
    def add_material_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_material_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_phase_component_material_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_phase_component_material_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_component_material_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_component_material_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_element_material_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_element_material_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_material_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_material_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_energy_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_energy_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_phase_enthalpy_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_phase_enthalpy_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_enthalpy_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_enthalpy_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_phase_energy_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_phase_energy_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_energy_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_energy_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_momentum_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_momentum_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_phase_pressure_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_phase_pressure_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_pressure_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_pressure_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_phase_momentum_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_phase_momentum_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))

    def add_total_momentum_balances(self, *args, **kwargs):
        raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "add_total_momentum_balances. Please contact the "
                "developer of the ControlVolume class you are using."
                .format(self.name))
