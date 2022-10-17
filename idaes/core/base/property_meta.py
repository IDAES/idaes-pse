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

from pyomo.core.base.units_container import _PyomoUnit

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
    when instantiating the UnitSet - units are optional for current and luminous intensity.

    Units can be accesses by via either a property on the UnitSet (e.g., UnitSet.TIME) or
    via an index on the UnitSet (e.g., UnitSet["time"]).
    """

    # TODO: Test this in isolation
    def __init__(self, **kwds):
        self._time = kwds.pop("time", None)
        self._length = kwds.pop("length", None)
        self._mass = kwds.pop("mass", None)
        self._amount = kwds.pop("amount", None)
        self._temperature = kwds.pop("temperature", None)
        self._current = kwds.pop("current", None)
        self._luminous_intensity = kwds.pop("luminous_intensity", None)
        if kwds:
            raise PropertyPackageError(
                f"Unrecognized base quantities: {[i for i in kwds.keys()]}"
            )

        # Check that valid units were assigned
        for q in [
            "time",
            "length",
            "mass",
            "amount",
            "temperature",
            "current",
            "luminous_intensity",
        ]:
            u = getattr(self, "_" + q)
            if u is None:
                if q in ["luminous_intensity", "current"]:
                    # these units are infrequently used in PSE, so allow users
                    # to skip these
                    continue
                else:
                    raise PropertyPackageError(
                        f"Units of measurement not provided for base quantity {q}. Units must be provided "
                        "for all base quantities except for current and luminous intensity."
                    )
            elif not isinstance(u, _PyomoUnit):
                raise PropertyPackageError(
                    f"Unrecognized units of measurement for quantity {q} ({u})"
                )
        # TODO: Could add a check for dimensionality of base units (i.e., units for time are a measure of time)

    def __getitem__(self, key):
        # Check to catch cases where luminous intensity has a space
        if key == "luminous intensity":
            key = "luminous_intensity"

        try:
            return getattr(self, key.upper())
        except AttributeError:
            raise PropertyPackageError(
                f"Unrecognised quantity {key}. Please check that this is a recognised quantity "
                "defined in idaes.core.base.property_meta.UnitSet."
            )

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
        self._default_units = UnitSet(**u)

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
