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
from idaes.core import ProcessBlockData
from idaes.core.util.config import is_property_parameter_block

__author__ = "Andrew Lee"

# Set up logger
logger = logging.getLogger('idaes.unit_model.holdup')


# Enumerate options for material balances
MaterialBalanceType = Enum(
    'None',
    'ComponentPhase',
    'ComponentTotal',
    'ElementTotal',
    'Total')

# Enumerate options for energy balances
EnergyBalanceType = Enum(
    'None',
    'EnthalpyPhase',
    'EnthalpyTotal',
    'EnergyPhase',
    'EnergyTotal')

# Enumerate options for momentum balances
MomentumBalanceType = Enum(
    'None',
    'PressureTotal',
    'PressurePhase',
    'MomentumTotal',
    'MomentumPhase')

# Set up example ConfigBlock that will work with ControlVolume autobuild method
CONFIG_Base = ProcessBlockData.CONFIG()
CONFIG_Base.declare("dynamic", ConfigValue(
        default='use_parent_value',
        domain=In(['use_parent_value', True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic,
**default** - 'use_parent_value'.  **Valid values:** {
**'use_parent_value'** - get flag from parent,
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
        default=MaterialBalanceType.ComponentPhase,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.ComponentPhase. **Valid values:** {
**MaterialBalanceType.None** - exclude material balances,
**MaterialBalanceType.ComponentPhase - use phase component balances,
**MaterialBalanceType.ComponentTotal - use total component balances,
**MaterialBalanceType.ElementTotal - use total element balances,
**MaterialBalanceType.Total - use total material balance.}"""))
CONFIG_Base.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.EnthalpyPhase,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.EnthalpyPhase. **Valid values:** {
**EnergyBalanceType.None** - exclude energy balances,
**EnergyBalanceType.EnthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.EnthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.EnergyTotal** - single energy balance for material,
**EnergyBalanceType.EnergyPhase** - energy balances for each phase.}"""))
CONFIG_Base.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.PressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.PressureTotal. **Valid values:** {
**MomentumBalanceType.None** - exclude momentum balances,
**MomentumBalanceType.PressureTotal** - single pressure balance for material,
**MomentumBalanceType.PressurePhase** - pressure balances for each phase,
**MomentumBalanceType.MomentumTotal** - single momentum balance for material,
**MomentumBalanceType.MomentumPhase** - momentum balances for each phase.}"""))
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
        domain=is_property_parameter_block,
        description="Property package to use for holdup",
        doc="""Property parameter object used to define property calculations,
**default** - None. **Valid values:** {
**a ParameterBlock object**.}"""))
CONFIG_Base.declare("property_package_args", ConfigValue(
        domain=ConfigBlock,
        description="Arguments to use for constructing property packages",
        doc="""A dict of arguments to be passed to a property block and used
when constructing these, **default** - empty ConfigBlock. **Valid values:** {
**ConfigBlock** - see property package for documentation.}"""))


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
        domain=In([True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic, **default** -
'use_parent_value'. **Valid values:** {
**'use_parent_value'** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}"""))
    CONFIG.declare("property_package", ConfigValue(
        domain=is_property_parameter_block,
        description="Property package to use for holdup",
        doc="""Property parameter object used to define property calculations.
**Valid values:** a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigValue(
        domain=ConfigBlock,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
 and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """
        General build method for Holdup blocks. This method calls a number
        of sub-methods which automate the construction of expected attributes
        of all Holdup blocks.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        pass
#        # Check construction flags
#        self._inherit_default_build_arguments()
#
#        # Call UnitBlockData.build to setup dynamics
#        self._setup_dynamics()
#
#        # Get property pacakge details
#        self.get_property_package()

#    def _inherit_default_build_arguments(self):
#        """
#        Method to collect build arguments from parent blocks as required.
#
#        If an argument is set as 'use_parent_value', this method attempts to
#        set the argument to that of the parent model, otherwise a default value
#        is used.
#
#        Args:
#            None
#
#        Returns:
#            None
#        """
#        parent = self.parent_block()
#        build_arguments = {'property_package': None,
#                           'property_package_args': {},
#                           'include_holdup': True,
#                           'material_balance_type': 'component_phase',
#                           'energy_balance_type': 'enthalpy_total',
#                           'momentum_balance_type': 'pressure',
#                           'has_rate_reactions': False,
#                           'has_equilibrium_reactions': False,
#                           'has_phase_equilibrium': False,
#                           'has_mass_transfer': False,
#                           'has_heat_transfer': False,
#                           'has_work_transfer': False,
#                           'has_pressure_change': False
#                           }
#        for arg in build_arguments:
#            # If argument is 'use_parent_value', try to get from parent
#            if getattr(self.config, arg) is 'use_parent_value':
#                try:
#                    setattr(self.config, arg, getattr(parent.config, arg))
#                except AttributeError:
#                    # If parent does not have flag, resort to defaults
#                    self.config[arg] = build_arguments[arg]
#
#    def _setup_dynamics(self):
#        """
#        This method automates the setting of the dynamic flag and time domain
#        for holdup blocks.
#
#        If dynamic flag is 'use_parent_value', method attempts to get the value
#        of the dynamic flag from the parent model, otherwise the local vlaue is
#        used. The time domain is aloways collected from the parent model.
#
#        Finally, if dynamic = True, the include_holdup flag is checked to
#        ensure it is also True.
#
#        Args:
#            None
#
#        Returns:
#            None
#        """
#        # Check the dynamic flag, and retrieve if necessary
#        if self.config.dynamic == 'use_parent_value':
#            # Get dynamic flag from parent
#            try:
#                self.config.dynamic = self.parent_block().config.dynamic
#            except AttributeError:
#                # If parent does not have dynamic flag, raise Exception
#                raise AttributeError('{} has a parent model '
#                                     'with no dynamic attribute.'
#                                     .format(self.name))
#
#        # Try to get reference to time object from parent
#        try:
#            object.__setattr__(self, "time", self.parent_block().time)
#        except AttributeError:
#            raise AttributeError('{} has a parent model '
#                                 'with no time domain'.format(self.name))
#
#        # Check include_holdup, if present
#        if self.config.dynamic:
#            if not self.config.include_holdup:
#                # Dynamic model must have include_holdup = True
#                logger.warning('{} Dynamic holdup blocks must have '
#                               'include_holdup = True. '
#                               'Overwritting argument.'
#                               .format(self.name))
#                self.config.include_holdup = True
#
#    def get_property_package(self):
#        """
#        This method gathers the necessary information about the property
#        package to be used in the holdup block.
#
#        If a property package has not been provided by the user, the method
#        searches up the model tree until it finds an object with the
#        'default_property_package' attribute and uses this package for the
#        holdup block.
#
#        The method also gathers any default construction arguments specified
#        for the property package and combines these with any arguments
#        specified by the user for the holdup block (user specified arguments
#        take priority over defaults).
#
#        Args:
#            None
#
#        Returns:
#            None
#        """
#        # Get property_package block if not provided in arguments
#        parent = self.parent_block()
#        if self.config.property_package in (None, 'use_parent_value'):
#            # Try to get property_package from parent
#            if hasattr(parent.config, "property_package"):
#                if parent.config.property_package is None:
#                    parent.config.property_package = \
#                        self._get_default_prop_pack()
#
#                self.config.property_package = parent.config.property_package
#            else:
#                self.config.property_package = self._get_default_prop_pack()
#
#        # Get module of property package
#        self.property_module = self.config.property_package.property_module
#
#        # Check for any flowsheet level build arguments
#        try:
#            # If flowsheet arguments exist, merge with local arguments
#            # Local arguments take precedence
#            arg_dict = self.config.property_package.config.default_arguments
#            arg_dict.update(self.config.property_package_args)
#            self.config.property_package_args = arg_dict
#        except AttributeError:
#            # Otherwise, just use local arguments
#            pass
#
#    def _get_default_prop_pack(self):
#        """
#        This method is used to find a default property package defined at the
#        flowsheet level if a package is not provided as an argument when
#        instantiating the holdup block.
#
#        Args:
#            None
#
#        Returns:
#            None
#        """
#        parent = self.parent_block()
#        while True:
#            if hasattr(parent.config, "default_property_package"):
#                break
#            else:
#                if parent.parent_block() is None:
#                    raise AttributeError(
#                            '{} no property package provided and '
#                            'no default defined. Found end of '
#                            'parent tree.'.format(self.name))
#                elif parent.parent_block() == parent:
#                    raise ValueError(
#                            '{} no property package provided and '
#                            'no default defined. Found recursive '
#                            'loop in parent tree.'.format(self.name))
#                parent = parent.parent_block()
#
#        logger.info('{} Using default property package'
#                    .format(self.name))
#
#        if parent.config.default_property_package is None:
#            raise ValueError('{} no default property package has been '
#                             'specified at flowsheet level ('
#                             'default_property_package = None)'
#                             .format(self.name))
#
#        return parent.config.default_property_package
#
#    def _get_indexing_sets(self):
#        """
#        This method collects all necessary indexing sets from property
#        parameter block and makes references to these for use within the holdup
#        block. Collected indexing sets are: phase_list, component_list,
#        rate_reaction_idx, equilibrium_reaction_idx, and element_list.
#
#        Args:
#            None
#
#        Returns:
#            None
#        """
#        # Get outlet property block object
#        try:
#            prop_block = self.properties_out[0]
#        except AttributeError:
#            try:
#                prop_block = self.properties[0, self.ldomain.last()]
#            except AttributeError:
#                prop_block = self.properties[0]
#
#        # Get phase and component list(s)
#        try:
#            add_object_ref(self, "phase_list",
#                           prop_block.phase_list)
#            add_object_ref(self, "component_list",
#                           prop_block.component_list)
#        except AttributeError:
#            raise AttributeError('{} property_package provided does not appear'
#                                 ' to be a valid property package (does not '
#                                 'contain a component_list and/or phase_list)'
#                                 .format(self.name))
#
#        # Get reaction indices and stoichiometry
#        if self.config.has_rate_reactions:
#            try:
#                add_object_ref(self, "rate_reaction_idx",
#                               prop_block.rate_reaction_idx)
#                add_object_ref(self,
#                               "rate_reaction_stoichiometry",
#                               prop_block.rate_reaction_stoichiometry)
#            except AttributeError:
#                self.config.has_rate_reactions = False
#                logger.info('{} property package does not support kinetic '
#                            'reactions. has_rate_reactions set to False.'
#                            .format(self.name))
#
#        if self.config.has_equilibrium_reactions:
#            try:
#                add_object_ref(self,
#                               "equilibrium_reaction_idx",
#                               prop_block.equilibrium_reaction_idx)
#                add_object_ref(self,
#                               "equilibrium_reaction_stoichiometry",
#                               prop_block.equilibrium_reaction_stoichiometry)
#            except AttributeError:
#                self.config.has_equilibrium_reactions = False
#                logger.info('{} property package does not support equilibrium '
#                            'reactions. has_equilibrium_reactions set to '
#                            'False.'.format(self.name))
#
#        if self.config.has_phase_equilibrium:
#            try:
#                add_object_ref(self,
#                               "phase_equilibrium_idx",
#                               prop_block.phase_equilibrium_idx)
#                add_object_ref(self,
#                               "phase_equilibrium_list",
#                               prop_block.phase_equilibrium_list)
#            except AttributeError:
#                self.config.has_phase_equilibrium = False
#                logger.info('{} property package does not support phase '
#                            'equilibrium has_phase_equilibrium set to '
#                            'False.'.format(self.name))
#
#        # If element balances, check properties for list of elements
#        if self.config.material_balance_type == 'element_total':
#            try:
#                add_object_ref(self, "element_list",
#                               prop_block.element_list)
#            except AttributeError:
#                raise AttributeError("{} Selected property package does not "
#                                     "support elemental mass balances"
#                                     .format(self.name))
