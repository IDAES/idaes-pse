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
from pyomo.environ import units
from pyomo.core.base.units_container import _PyomoUnit, InconsistentUnitsError

from idaes.core.util.exceptions import PropertyPackageError
import idaes.logger as idaeslog

__author__ = "Dan Gunter <dkgunter@lbl.gov>, Andrew Lee"

_log = idaeslog.getLogger(__name__)


class HasPropertyClassMetadata(object):
    """Interface for classes that have PropertyClassMetadata."""

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


class UnitSet(object):
    """
    Object defining the set of recognised quantities in IDAES and their base units.

    Units of measurement are defined by setting units for the seven base SI quantities
    (amount, current, length, luminous intensity, mass, temperature and time), from which units
    for all other quantities are derived. The units of the seven base quantities must be provided
    when instantiating the UnitSet, otherwise base SI units are assumed.

    Units can be accesses by via either a property on the UnitSet (e.g., UnitSet.TIME) or
    via an index on the UnitSet (e.g., UnitSet["time"]).
    """

    _base_quantities = {
        "AMOUNT": units.mol,
        "CURRENT": units.watt,
        "LENGTH": units.meter,
        "LUMINOUS_INTENSITY": units.candela,
        "MASS": units.kilogram,
        "TEMPERATURE": units.kelvin,
        "TIME": units.seconds,
    }

    def __init__(
        self,
        amount: _PyomoUnit = units.mol,
        current: _PyomoUnit = units.watt,
        length: _PyomoUnit = units.meter,
        luminous_intensity: _PyomoUnit = units.candela,
        mass: _PyomoUnit = units.kilogram,
        temperature: _PyomoUnit = units.kelvin,
        time: _PyomoUnit = units.seconds,
    ):
        self._time = time
        self._length = length
        self._mass = mass
        self._amount = amount
        self._temperature = temperature
        self._current = current
        self._luminous_intensity = luminous_intensity

        # Check that valid units were assigned
        for q, expected_dim in self._base_quantities.items():
            u = getattr(self, q)
            if not isinstance(u, _PyomoUnit):
                # Check for non-unit inputs from user
                raise PropertyPackageError(
                    f"Unrecognized units of measurement for quantity {q} ({u})"
                )

            # Check for expected dimensionality
            try:
                # Try to convert user-input to SI units of expected dimensions
                units.convert(u, expected_dim)
            except InconsistentUnitsError:
                # An error indicates a mismatch in units or the units registry
                raise PropertyPackageError(
                    f"Invalid units of measurement for quantity {q} ({u}). "
                    "Please ensure units provided are valid for this quantity and "
                    "use the Pyomo unit registry."
                )

    def __getitem__(self, key):
        try:
            # Check to catch cases where luminous intensity has a space
            return getattr(self, key.upper().replace(" ", "_"))
        except AttributeError:
            raise PropertyPackageError(
                f"Unrecognised quantity {key}. Please check that this is a recognised quantity "
                "defined in idaes.core.base.property_meta.UnitSet."
            )

    def unitset_is_consistent(self, other):
        return all(getattr(self, q) is getattr(other, q) for q in self._base_quantities)

    @property
    def TIME(self):
        return self._time

    @property
    def LENGTH(self):
        return self._length

    @property
    def MASS(self):
        return self._mass

    @property
    def AMOUNT(self):
        return self._amount

    @property
    def TEMPERATURE(self):
        return self._temperature

    @property
    def CURRENT(self):
        return self._current

    @property
    def LUMINOUS_INTENSITY(self):
        return self._luminous_intensity

    # Length based
    @property
    def AREA(self):
        return self._length**2

    @property
    def VOLUME(self):
        return self._length**3

    # Flows
    @property
    def FLOW_MASS(self):
        return self._mass * self._time**-1

    @property
    def FLOW_MOLE(self):
        return self._amount * self._time**-1

    @property
    def FLOW_VOL(self):
        return self._length**3 * self._time**-1

    @property
    def FLUX_MASS(self):
        return self._mass * self._time**-1 * self._length**-2

    @property
    def FLUX_MOLE(self):
        return self._amount * self._time**-1 * self._length**-2

    @property
    def FLUX_ENERGY(self):
        return self._mass * self._time**-3

    # Velocity and Acceleration
    @property
    def VELOCITY(self):
        return self._length * self._time**-1

    @property
    def ACCELERATION(self):
        return self._length * self._time**-2

    # Pressures
    @property
    def PRESSURE(self):
        return self._mass * self._length**-1 * self._time**-2

    @property
    def GAS_CONSTANT(self):
        return (
            self._mass
            * self._length**2
            * self._time**-2
            * self._temperature**-1
            * self._amount**-1
        )

    # Densities
    @property
    def DENSITY_MASS(self):
        return self._mass * self._length**-3

    @property
    def DENSITY_MOLE(self):
        return self._amount * self._length**-3

    @property
    def MOLECULAR_WEIGHT(self):
        return self._mass / self._amount

    # Energy
    @property
    def ENERGY(self):
        return self._mass * self._length**2 * self._time**-2

    @property
    def ENERGY_MASS(self):
        return self._length**2 * self._time**-2

    @property
    def ENERGY_MOLE(self):
        return self._mass * self._length**2 * self._time**-2 * self._amount**-1

    @property
    def POWER(self):
        return self._mass * self._length**2 * self._time**-3

    # Heat Related
    @property
    def HEAT_CAPACITY_MASS(self):
        return self._length**2 * self._time**-2 * self._temperature**-1

    @property
    def HEAT_CAPACITY_MOLE(self):
        return (
            self._mass
            * self._length**2
            * self._time**-2
            * self._temperature**-1
            * self._amount**-1
        )

    @property
    def HEAT_TRANSFER_COEFFICIENT(self):
        return self._mass * self._time**-3 * self._temperature**-1

    # Entropy
    @property
    def ENTROPY(self):
        return (
            self._mass * self._length**2 * self._time**-2 * self._temperature**-1
        )

    @property
    def ENTROPY_MASS(self):
        return self._length**2 * self._time**-2 * self._temperature**-1

    @property
    def ENTROPY_MOLE(self):
        return (
            self._mass
            * self._length**2
            * self._time**-2
            * self._temperature**-1
            * self._amount**-1
        )

    # Transport Properties
    @property
    def DYNAMIC_VISCOSITY(self):
        return self._mass * self._length**-1 * self._time**-1

    @property
    def THERMAL_CONDUCTIVITY(self):
        return self._mass * self._length * self._time**-3 * self._temperature**-1


class PropertyClassMetadata(object):
    """Container for metadata about the property class, which includes
       default units and properties.

    Example usage::

            foo = PropertyClassMetadata()
            foo.add_default_units(time = pyo.units.fortnights,
                                  mass = pyo.units.stones)
            foo.add_properties({'under_sea': {'units': 'leagues'},
                                'tentacle_size': {'units': 'yards'}})
            foo.add_required_properties({'under_sea': 'leagues',
                                        'tentacle_size': 'yards'})

    """

    def __init__(self):
        # TODO: Deprecate in favour of common units property
        self._default_units = None
        self._properties = {}
        self._required_properties = {}

    @property
    def default_units(self):
        # TODO: Deprecate in favour of common units property
        return self._default_units

    @property
    def derived_units(self):
        # TODO: Deprecate in favour of common units property
        return self._default_units

    @property
    def properties(self):
        return self._properties

    @property
    def required_properties(self):
        return self._required_properties

    def add_default_units(self, u):
        """Add a dict with keys for the base quantities used in the
        property package (as strings) and values of their default units as Pyomo unit objects.

        If units are not provided for a quantity, it will be assumed to use base SI unites.

        Args:
            u (dict): Key=property, Value=units

        Returns:
            None
        """
        # TODO: Could look at replacing dict with defined arguments
        # This would be a big API change
        try:
            self._default_units = UnitSet(**u)
        except TypeError:
            raise TypeError(
                f"Unexpected argument for base quantities found when creating UnitSet. "
                "Please ensure that units are only defined for the seven base quantities."
            )

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

    def get_derived_units(self, units):
        # TODO: Deprecate in favour of common units property
        return self.derived_units[units]


class PropertyMetadata(dict):
    """Container for property parameter metadata.

    Instances of this class are exactly dictionaries, with the
    only difference being some guidance on the values expected in the
    dictionary from the constructor.
    """

    def __init__(self, name=None, method=None, units=None):
        if name is None:
            raise TypeError('"name" is required')
        d = {"name": name, "method": method}
        if units is not None:
            d["units"] = units
        else:
            # Adding a default "null" unit in case it is not provided by user
            d["units"] = "-"
        super(PropertyMetadata, self).__init__(d)
