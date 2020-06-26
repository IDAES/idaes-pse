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
These classes handle the metadata aspects of classes representing
property packages.

Implementors of property packages need to do the following:

1. Create a new class that inherits from
   :class:`idaes.core.property_base.PhysicalParameterBlock`, which in turn
   inherits from :class:`HasPropertyClassMetadata`, in this module.

2. In that class, implement the `define_metadata()` method, inherited from
   :class:`HasPropertyClassMetadata`. This method is called
   automatically, once, when the `get_metadata()` method is first invoked.
   An empty metadata object (an instance of :class:`PropertyClassMetadata`)
   will be passed in, which the method should populate with information about
   properties and default units.

Example::

    from idaes.core.property_base import PhysicalParameterBlock

    class MyPropParams(PhysicalParameterBlock):

        @classmethod
        def define_metadata(cls, meta):
            meta.add_default_units({foo.U.TIME: 'fortnights',
                                   foo.U.MASS: 'stones'})
            meta.add_properties({'under_sea': {'units': 'leagues'},
                                'tentacle_size': {'units': 'yards'}})
            meta.add_required_properties({'under_sea': 'leagues',
                                'tentacle_size': 'yards'})

        # Also, of course, implement the non-metadata methods that
        # do the work of the class.

"""

from pyomo.core.base.units_container import _PyomoUnit

from idaes.core.util.exceptions import PropertyPackageError
import idaes.logger as idaeslog

__author__ = 'Dan Gunter <dkgunter@lbl.gov>, Andrew Lee'

_log = idaeslog.getLogger(__name__)


class HasPropertyClassMetadata(object):
    """Interface for classes that have PropertyClassMetadata.
    """
    _metadata = None

    @classmethod
    def get_metadata(cls):
        """Get property parameter metadata.

        If the metadata is not defined, this will instantiate a new
        metadata object and call `define_metadata()` to set it up.

        If the metadata is already defined, it will be simply returned.

        Returns:
            PropertyClassMetadata: The metadata
        """
        if cls._metadata is None:
            pcm = PropertyClassMetadata()
            cls.define_metadata(pcm)
            cls._metadata = pcm
        return cls._metadata

    @classmethod
    def define_metadata(cls, pcm):
        """Set all the metadata for properties and units.

        This method should be implemented by subclasses.
        In the implementation, they should set information into the
        object provided as an argument.

        Args:
            pcm (PropertyClassMetadata): Add metadata to this object.

        Returns:
            None
        """
        raise NotImplementedError()


class UnitNames(object):
    """Names for recognized units.
    """
    TM = TIME = 'time'
    LN = LENGTH = 'length'
    MA = MASS = 'mass'
    AM = AMOUNT = 'amount'
    TE = TEMPERATURE = 'temperature'
    CU = CURRENT = 'current'
    LI = LUMINOUS_INTENSITY = 'luminous intensity'


_base_units = {
    'time': None,
    'length': None,
    'mass': None,
    'amount': None,
    'temperature': None,
    'current': None,
    'luminous intensity': None}


class PropertyClassMetadata(object):
    """Container for metadata about the property class, which includes
       default units and properties.

    Example usage::

            foo = PropertyClassMetadata()
            foo.add_default_units({foo.U.TIME: 'fortnights',
                                   foo.U.MASS: 'stones'})
            foo.add_properties({'under_sea': {'units': 'leagues'},
                                'tentacle_size': {'units': 'yards'}})
            foo.add_required_properties({'under_sea': 'leagues',
                                        'tentacle_size': 'yards'})

    """
    #: Alias for class enumerating supported/known unit types
    U = UnitNames

    def __init__(self):
        self._default_units = _base_units.copy()
        self._properties = {}
        self._required_properties = {}

    @property
    def default_units(self):
        return self._default_units

    @property
    def properties(self):
        return self._properties

    @property
    def required_properties(self):
        return self._required_properties

    def add_default_units(self, u):
        """Add a dict with keys for the
        quantities used in the property package (as strings) and values of
        their default units as unit objects or strings.

        The quantities used by the framework are in constants
        defined in :class:`UnitNames`, aliased here in the class
        attribute `U`.

        Args:
            u (dict): Key=property, Value=units

        Returns:
            None
        """
        self._default_units.update(u)

        # Validate values. Pyomo units are all-or-nothing, so check to see that
        # this is the case
        _units = 0
        for q, u in self._default_units.items():
            if isinstance(u, _PyomoUnit):
                _units += 1
            elif u is None and (q == "luminous intensity" or q == "current"):
                # these units are infrequently used in PSE, so allow users
                # to skip these
                continue
            elif _units > 0:
                # Mix of units and non-unit objects
                raise PropertyPackageError(
                    "default_units ({}: {}): if using Pyomo Units objects, "
                    "all units must be defined using Units objects (not "
                    "compount units)."
                    .format(q, u))

        # Take opportunity to log a deprecation warning if units are not used
        # TODO : Add back in after June 2020 release
        # if _units == 0:
        #     _log.warning("DEPRECATED: IDAES is moving to using Pyomo Units "
        #                  "when defining default units, which are used "
        #                  "to automatically determine units of measurement "
        #                  "for quantities and convert where necessary. "
        #                  "Users are strongly encouraged to convert their "
        #                  "property packages to use Pyomo Units obejcts.")

    def add_properties(self, p):
        """Add properties to the metadata.

        For each property, the value should be another dict which may contain
        the following keys:

        - 'method': (required) the name of a method to construct the
                    property as a str, or None if the property will be
                    constructed by default.
        - 'units': (optional) units of measurement for the property.

        Args:
            p (dict): Key=property, Value=PropertyMetadata or equiv. dict

        Returns:
            None
        """
        for k, v in p.items():
            if not isinstance(v, PropertyMetadata):
                v = PropertyMetadata(name=k, **v)
            self._properties[k] = v

    def add_required_properties(self, p):
        """Add required properties to the metadata.

        For each property, the value should be the expected units of
        measurement for the property.

        Args:
            p (dict): Key=property, Value=units

        Returns:
            None
        """
        # Using the same PropertyMetadata class as for units, but 'method'
        # will always be none
        for k, v in p.items():
            if not isinstance(v, PropertyMetadata):
                v = PropertyMetadata(name=k, units=v)
            self._required_properties[k] = v


class PropertyMetadata(dict):
    """Container for property parameter metadata.

    Instances of this class are exactly dictionaries, with the
    only difference being some guidance on the values expected in the
    dictionary from the constructor.
    """

    def __init__(self, name=None, method=None, units=None):
        if name is None:
            raise TypeError('"name" is required')
        d = {'name': name, 'method': method}
        if units is not None:
            d['units'] = units
        else:
            # Adding a default "null" unit in case it is not provided by user
            d['units'] = "-"
        super(PropertyMetadata, self).__init__(d)
