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
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

from pyomo.environ import units
from pyomo.core.base.units_container import _PyomoUnit, InconsistentUnitsError
from pyomo.common.deprecation import deprecation_warning

from idaes.core.util.exceptions import PropertyPackageError
from idaes.core.base.property_set import StandardPropertySet, PropertySetBase
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

            # Check that the metadata was actually populated
            # Check requires looking at private attributes
            # pylint: disable-next=protected-access
            if pcm._properties is None or pcm._default_units is None:
                raise PropertyPackageError(
                    "Property package did not populate all expected metadata."
                )
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
        "CURRENT": units.ampere,
        "LENGTH": units.meter,
        "LUMINOUS_INTENSITY": units.candela,
        "MASS": units.kilogram,
        "TEMPERATURE": units.kelvin,
        "TIME": units.seconds,
    }

    def __init__(self):
        self._time = units.seconds
        self._length = units.meter
        self._mass = units.kilogram
        self._amount = units.mole
        self._temperature = units.kelvin
        self._current = units.ampere
        self._luminous_intensity = units.candela

    def set_units(
        self,
        amount: _PyomoUnit = units.mol,
        current: _PyomoUnit = units.ampere,
        length: _PyomoUnit = units.meter,
        luminous_intensity: _PyomoUnit = units.candela,
        mass: _PyomoUnit = units.kilogram,
        temperature: _PyomoUnit = units.kelvin,
        time: _PyomoUnit = units.seconds,
    ):
        """
        Set desired units of measurement for the seven base quantities.

        Args:
            amount: units for amount (default = moles)
            current: units for current (default = Amperes)
            length: units for length (default = meters)
            luminous_intensity: units for luminous intensity (default = candela)
            mass: units for mass (default = kilograms)
            temperature: units for temperature (default = Kelvins)
            time: units for time (default = seconds)

        Returns:
            None
        """
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

    def __getitem__(self, key: str):
        try:
            # Check to catch cases where luminous intensity has a space
            return getattr(self, key.upper().replace(" ", "_"))
        except AttributeError:
            raise PropertyPackageError(
                f"Unrecognised quantity {key}. Please check that this is a recognised quantity "
                "defined in idaes.core.base.property_meta.UnitSet."
            )

    def unitset_is_consistent(self, other: "UnitSet"):
        """
        Checks that defined units of measurement for base quantities are consistent with those
        in other UnitSet.

        Args:
            other: UnitSet to check for consistency with

        Returns:
            Bool indicating whether units are consistent
        """
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

    @property
    def MOLAR_VOLUME(self):
        return self._length**3 * self._amount**-1

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

    # Velocity, Acceleration and Force
    @property
    def VELOCITY(self):
        return self._length * self._time**-1

    @property
    def ACCELERATION(self):
        return self._length * self._time**-2

    @property
    def FORCE(self):
        return self._length * self._mass * self._time**-2

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

    # Densities & Concentrations
    @property
    def DENSITY_MASS(self):
        return self._mass * self._length**-3

    @property
    def DENSITY_MOLE(self):
        return self._amount * self._length**-3

    @property
    def MOLALITY(self):
        return self._amount * self._mass

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

    @property
    def VOLTAGE(self):
        return self._mass * self._length**2 * self._time**-3 * self._current**-1

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
    def DIFFUSIVITY(self):
        return self._length**2 * self._time**-1

    @property
    def DYNAMIC_VISCOSITY(self):
        return self._mass * self._length**-1 * self._time**-1

    @property
    def KINEMATIC_VISCOSITY(self):
        return self._length**2 * self._time**-1

    @property
    def SURFACE_TENSION(self):
        return self._mass * self._time**-2

    @property
    def THERMAL_CONDUCTIVITY(self):
        return self._mass * self._length * self._time**-3 * self._temperature**-1


class PropertyClassMetadata(object):
    """
    Container for metadata about the property class, which includes
    default units and properties.

    Example usage::

        foo = PropertyClassMetadata()
        foo.add_default_units(time = pyo.units.fortnights,
                              mass = pyo.units.stones)
        foo.add_properties({'under_sea': {'method': 'submarine', 'units': 'leagues', 'required': False, 'supported': True},
                            'tentacle_size': {'method': 'kraken', 'units': 'yards', 'required': True, 'supported': True}})

    """

    def __init__(self):
        # TODO: Deprecate in favour of common units property
        self._default_units = UnitSet()
        # Assume a default PropertySet to begin with. Property packages can replace this
        # with more specialized forms if required
        self._properties = StandardPropertySet(parent=self)

    def define_property_set(self, propset: PropertySetBase):
        """
        Define the type of property set to use for this package.

        Args:
            propset: PropertySet class (must derive from PropertySetBase)

        Returns:
            None
        """
        if not issubclass(propset, PropertySetBase):
            raise PropertyPackageError(
                f"{propset} does not derive from IDAES PropertySetBase class."
            )
        self._properties = propset(parent=self)

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

    def add_default_units(self, u: dict):
        """
        Set units of measurement for base quantities used in this property package. Units
        should be provided as a dict with keys being the seven base quantities and values
        being Pyomo unit expressions. These will be used to update the UnitSet associated
        with this property package.

        Args:
            u (dict): Key=property, Value=units

        Returns:
            None

        Raises:
            TypeError if definitions for unexpected quantities are found
        """
        # TODO: Could look at replacing dict with defined arguments
        # This would be a big API change
        try:
            self._default_units.set_units(**u)
        except TypeError:
            raise TypeError(
                "Unexpected argument for base quantities found when creating UnitSet. "
                "Please ensure that units are only defined for the seven base quantities."
            )

    def add_properties(self, p: dict):
        """Add properties to the metadata.

        For each property, the value should be another dict which may contain
        the following keys:

        - 'units': (optional) units of measurement for the property.
        - 'indices': (optional) list of sub-property indices for this property. If None, use default set, if False unindexed.
        - 'method': (optional, only if 'indices' is None or False) the name of a method to construct the
                    property as a str, or None if the property will be
                    constructed by default.
        - 'supported': (optional, only if 'indices' is None or False) bool indicating if this property is
                       supported by this package.
        - 'required': (optional, only if 'indices' is None or False) bool indicating if this property is
                      required by this package.
        - 'valid_range': (optional, only if 'indices' is None or False) 2-tuple containing range of validity for
                      property values (lower, upper).
        - 'initialize': (optional) dict indicating 'method', 'required', 'supported' and 'valid_range' values for sub-properties by index.

        Args:
            p (dict): Key=property, Value=dict

        Returns:
            None
        """
        # TODO: Deprecate in favour of directly updating or adding metadata
        for k, v in p.items():
            units = v.pop("units", None)
            try:
                try:
                    n, i = self._properties.get_name_and_index(k)
                except ValueError:
                    msg = (
                        f"The property name {k} in property metadata is not a recognized "
                        "standard property name defined in this PropertySet. Please refer "
                        "to IDAES standard names in the IDAES documentation. You can use "
                        "the define_custom_properties() rather than the add_properties() "
                        "method to define metadata for this property. You can also use a "
                        "different property set by calling the define_property_set() method."
                    )
                    deprecation_warning(
                        msg=msg, logger=_log, version="2.0.0", remove_in="3.0.0"
                    )
                    n = k
                    i = None
                getattr(self._properties, n)[i].update_property(**v)
            except AttributeError:
                # TODO: Deprecate this and make it raise an exception if an unknown property is encountered
                # Force users to explicitly declare new/custom properties
                self._properties.define_property(name=k, **v, units=units)

    def define_custom_properties(self, p: dict):
        """Add custom properties to the metadata.

        For each property, the value should be another dict which may contain
        the following keys:

        - 'units': (optional) units of measurement for the property.
        - 'indices': (optional) list of sub-property indices for this property. If None, use default set, if False unindexed.
        - 'method': (optional, only if 'indices' is None or False) the name of a method to construct the
                    property as a str, or None if the property will be
                    constructed by default.
        - 'supported': (optional, only if 'indices' is None or False) bool indicating if this property is
                       supported by this package.
        - 'required': (optional, only if 'indices' is None or False bool indicating if this property is
                      required by this package.
        - 'initialize': (optional) dict indicating 'method', 'required' and 'supported' values for sub-properties by index.

        Args:
            p (dict): Key=property, Value=dict

        Returns:
            None
        """
        for k, v in p.items():
            self._properties.define_property(name=k, **v)

    def add_required_properties(self, p: str):
        # TODO: Deprecate
        """Add required properties to the metadata.

        Update 'required' attribute of specified properties.
        Note that argument must be a dict for backwards compatibility.

        Args:
            p (dict): Key=property, Value=(ignored)

        Returns:
            None
        """
        for k in p.keys():
            try:
                self._properties[k].set_required(True)
            except KeyError:
                self._properties.define_property(name=k, supported=False, required=True)

    def get_derived_units(self, units: str):
        # TODO: Deprecate in favour of common units property
        return self.derived_units[units]
