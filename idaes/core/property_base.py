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
This module contains classes for property blocks and property parameter blocks.
"""
from __future__ import division

# Import Python libraries
import inspect
import logging

# Import Pyomo libraries
from pyomo.common.config import ConfigValue, In

# Other third-party
import six

# Import IDAES cores
from idaes.core.process_block import ProcessBlock
from idaes.core import ProcessBlockData
from idaes.core import property_meta

# Some more information about this module
__author__ = "Andrew Lee"

__all__ = ['PropertyBlockDataBase', 'PropertyParameterBase']

# Set up logger
logger = logging.getLogger('idaes.core')


class PropertyParameterBase(ProcessBlockData,
                            property_meta.HasPropertyClassMetadata):
    """
        This is the base class for property parameter blocks. These are blocks
        that contain a set of parameters associated with a specific property
        package, and are linked to by all instances of that property package.
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("default_arguments", ConfigValue(
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
        # Get module reference and store on block
        frm = inspect.stack()[1]
        self.property_module = inspect.getmodule(frm[0])


class PropertyBlockBase(ProcessBlock):
    """
        This is the base class for property block objects. These are used when
        constructing the SimpleBlock or IndexedBlock which will contain the
        PropertyData objects, and contains methods that can be applied to
        multiple PropertyBlockData objects simultaneously.
    """
    def initialize(self, *args):
        """
        This is a default initialization routine for PropertyBlocks to ensure
        that a routine is present. All PropertyBlockData classes should
        overload this method with one suited to the particular property package

        This method prints a warning and then ends.

        Args:
            None

        Returns:
            None
        """
        logger.warning('{} Property package does not define an initialization'
                       'method. PropertyBlock has not been initialized'
                       .format(self.name))


class PropertyBlockDataBase(ProcessBlockData):
    """
        This is the base class for property block data objects. These are
        blocks that contain the Pyomo components associated with calculating a
        set of thermophysical, transport and reacion properties for a given
        material.
    """
    # Create Class ConfigBlock
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("parameters", ConfigValue(
            description="""A reference to an instance of the Property Parameter
                         Block associated with this property package."""))
    CONFIG.declare("has_sum_fractions", ConfigValue(
            default=False,
            domain=In([True, False]),
            description="Sum of material fractions flag",
            doc="""Flag indicating whether this property block should include
                constraints for the sum of material fractions (e.g. sum of mole
                fractions == 1) (default=False).
                """))
    CONFIG.declare("calculate_equilibrium_reactions", ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Equilibrium reaction constraint flag",
            doc="""Flag indicating whether equilibrium reaction constraints
                should be constructed in this property block (default=True).
                """))
    CONFIG.declare("calculate_phase_equilibrium", ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Phase equilibrium constraint flag",
            doc="""Flag indicating whether phase equilibrium constraints
                should be constructed in this property block (default=True).
                """))

    def build(self):
        """
        General build method for PropertyBlockDatas. Inheriting models should
        call super().build.

        Args:
            None

        Returns:
            None
        """
        pass

    def __getattr__(self, attr):
        """
        This method is used to avoid generating unnecessary property
        calculations in property blocks. __getattr__ is called whenever a
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
                raise ValueError("{} Trying to remove call {} from __getattr__"
                                 " call list, however this is not the most "
                                 "recent call in the list ({}). This indicates"
                                 " a bug in the __getattr__ calls. Please "
                                 "contact the IDAES developers with this bug."
                                 .format(self.name,
                                         attr,
                                         self.__getattrcalls[-1]))

        # Check that attr is not something we shouldn't touch
        if attr == "domain" or attr.startswith("_"):
            # Don't interfere with anything by getting attributes that are
            # none of my business
            raise AttributeError('{} {} does not exist, but is a protected '
                                 'attribute. Check the naming of your '
                                 'components to avoid any reserved names'
                                 .format(self.name, attr))

        # Check for recursive calls
        try:
            # Check to see if attr already appears in call list
            if attr in self.__getattrcalls:
                # If it does, indicates a recursive loop.
                if attr == self.__getattrcalls[-1]:
                    # attr method is calling itself
                    self.__getattrcalls.append(attr)
                    raise Exception('{} _{} made a recursive call to '
                                    'itself, indicating a potential '
                                    'recursive loop. This is generally '
                                    'caused by the {} method failing to '
                                    'create the {} component.'
                                    .format(self.name, attr, attr, attr))
                else:
                    self.__getattrcalls.append(attr)
                    raise Exception('{} a potential recursive loop has been '
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
        except AttributeError:
            # Creat a list of calls if one does not already exist
            self.__getattrcalls = [attr]

        # Get property information from get_supported_properties
        try:
            m = self.config.parameters.get_supported_properties()

            if m is None:
                raise ValueError('{} Property package get_supported_properties'
                                 ' method returned None when trying to create '
                                 '{}. Please contact the developer of the '
                                 'property package'.format(self.name, attr))
        except KeyError:
            # If attr not in get_supported_properties, assume package does not
            # support property
            clear_call_list(self, attr)
            raise AttributeError('{} {} is not supported by property package '
                                 '(property is not listed in '
                                 'get_supported_properties).'
                                 .format(self.name, attr, attr))

        # Get method name from get_supported_properties
        try:
            if m[attr]['method'] is None:
                # If method is none, property should be constructed
                # by property package, so raise AttributeError
                clear_call_list(self, attr)
                raise AttributeError(
                        '{} {} should be constructed automatically '
                        'by property package, but is not present. '
                        'This can be caused by methods being called '
                        'out of order.'.format(self.name, attr))
            elif m[attr]['method'] is False:
                # If method is False, package does not support property
                # Raise NotImplementedError
                clear_call_list(self, attr)
                raise NotImplementedError(
                        '{} {} is not supported by property package '
                        '(property method is listed as False in '
                        'get_supported_properties).'
                        .format(self.name, attr))
            elif isinstance(m[attr]['method'], str):
                # Try to get method name in from PropertyBlock object
                try:
                    f = getattr(self, m[attr]['method'])
                except AttributeError:
                    # If fails, method does not exist
                    clear_call_list(self, attr)
                    raise AttributeError(
                            '{} {} get_supported_properties method '
                            'returned a name that does not correspond'
                            ' to any method in the property package. '
                            'Please contact the developer of the '
                            'property package.'.format(self.name, attr))
            else:
                # Otherwise method name is invalid
                clear_call_list(self, attr)
                raise ValueError('{} {} get_supported_properties method '
                                 'returned invalid value for method name. '
                                 'Please contact the developer of the '
                                 'property package.'
                                 .format(self.name, attr))
        except KeyError:
            # No method key - raise Exception
            # Need to use an AttributeError so Pyomo.DAE will handle this
            clear_call_list(self, attr)
            raise AttributeError('{} get_supported_properties method '
                                 'does not contain a method for {}. '
                                 'Please contact the developer of the '
                                 'property package.'.format(self.name, attr))

        # Call attribute if it is callable
        # If this fails, it should return a meaningful error.
        if callable(f):
            try:
                f()
            except:
                # Clear call list and reraise error
                clear_call_list(self, attr)
                raise
        else:
            # If f is not callable, inform the user and clear call list
            clear_call_list(self, attr)
            raise AttributeError('{} has an attribute _{}, however it is not '
                                 'callable.'.format(self.name, attr))

        # Clear call list, and return
        comp = getattr(self, attr)
        clear_call_list(self, attr)
        return comp
