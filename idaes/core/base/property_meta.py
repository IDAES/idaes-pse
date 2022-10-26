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

    def __init__(self):
        self._time = units.seconds
        self._length = units.meter
        self._mass = units.kilogram
        self._amount = units.mole
        self._temperature = units.kelvin
        self._current = units.watts
        self._luminous_intensity = units.candela

    def set_units(
        self,
        amount: _PyomoUnit = units.mol,
        current: _PyomoUnit = units.watt,
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
            current: units for current (default = Watts)
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


class PropertyMetadata(object):
    """
    Metadata object for defining a property.

    This object stores all the metadata associated with a single property, including:

        - standard name
        - units of measurement for this property (defined via associated UnitSet)
        - method that constructs this property and associated constraints (if build-on-demand)
        - whether property is supported by this package
        - whether this package expects this property to be provided by another package
    """

    # TODO: Add optional doc string to metadata objects
    def __init__(
        self, name=None, method=None, units=None, supported=False, required=False
    ):
        if name is None:
            raise TypeError('"name" is required')
        self._name = name
        self._method = method
        self._supported = supported
        self._required = required
        self._units = units  # TODO: Validate units are from UnitSet or dimensionless
        # TODO: For future, this would be the place to store default scaling information, etc.
        # TODO: Could also define default bounds, nominal values, etc.

    def __getitem__(self, key):
        try:
            return getattr(self, "_" + key)
        except AttributeError:
            # TODO: Real error message - needs to be a KeyError to work with getattr elsewhere
            raise KeyError()

    @property
    def name(self):
        """
        Standard name for property
        """
        return self._name

    @property
    def method(self):
        """
        Reference to method that can be called to construct this property and associated
        constraints if using build-on-demand approach.
        """
        return self._method

    @property
    def units(self):
        """
        Units of measurement for this property. This should be a reference to a quantity defined
        in the UnitSet associated with this property package.
        """
        return self._units

    @property
    def supported(self):
        """
        Bool indicating whether this property package supports calculation of this property.
        """
        return self._supported

    @property
    def required(self):
        """
        Bool indicating whether this property package requires calculation of this property
        by another property package.

        This is most commonly used by reaction packages which rely of thermophysical property
        packages to define other properties.
        """
        return self._required

    def set_method(self, meth):
        """
        Set method attribute of property.

        Args:
            meth: reference to method required to construct this property

        Returns:
            None
        """
        # TODO: Validate that meth is callable?
        self._method = meth

    def set_supported(self, supported=True):
        """
        Set supported attribute of property

        Args:
            supported: bool indicating whether package supports this property

        Returns:
            None
        """
        # TODO: Validate that supported is bool
        self._supported = supported

    def set_required(self, required=True):
        """
        Set required attribute of property

        Args:
            required: bool indicating whether package requires this property be defined by
            another property package

        Returns:
            None
        """
        # TODO: Validate that required is bool
        self._required = required

    def update_property(self, dict):
        """
        Update attributes of this property.

        Args:
            dict: containing desired attributes for this property. Supported keys are 'method',
            'required' and 'supported'.

        Returns:
            None

        Note that if not provided a value, 'supported` is assumed to be True.
        """
        if "method" in dict:
            self.set_method(dict["method"])
        if "required" in dict:
            self.set_required(dict["required"])
        if "supported" in dict:
            self.set_supported(dict["supported"])
        else:
            # Assume supported if not explicitly stated
            # TODO: Reconsider in the future, for now do this for backwards compatibility
            self.set_supported(True)


class PropertySet(object):
    """
    Metadata object which defines all the properties associated with a given property package.

    This object defines all the standard properties supported by IDAES, and also allows for
    definition of new properties if required.
    """

    # TODO: Add doc string

    def __init__(self, parent):
        self.__parent_block = parent

        self._define_standard_properties()

    def __getitem__(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            # TODO: Real error message - needs to be a KeyError to work with getattr elsewhere
            raise KeyError()

    def __iter__(self):
        for a in dir(self):
            aobj = getattr(self, a)
            if isinstance(aobj, PropertyMetadata):
                yield aobj

    def define_property(
        self, name=None, method=None, supported=True, required=False, units=None
    ):
        """
        Define a new property called `name`.

        Args:
            name: name of new property (required)
            method: reference to build-on-demand method for property (optional, default=None)
            supported: bool indicating if package supports this property (optional, default=True)
            required: bool indicating if package requires this property from another package (optional, default=False)
            units: quantity defined in associated UnitSet defining the units of measurement for property

        Returns:
            None
        """
        if hasattr(self, name):
            raise PropertyPackageError(
                f"A property with the name {name} already exists. Please use update_property "
                "method if you wish to update an existing property's metadata."
            )

        setattr(
            self,
            name,
            PropertyMetadata(
                name=name,
                method=method,
                supported=supported,
                required=required,
                units=units,
            ),
        )

    def check_required_properties(self, other):
        """
        Check that other property package supports all properties marked as required by this package.

        Args:
            other: PropertySet to check for supported properties

        Returns:
            list of properties required by this package which are not supported by other package
        """
        unsupported = []
        for a in dir(self):
            aobj = getattr(self, a)
            if isinstance(aobj, PropertyMetadata) and aobj.required:
                try:
                    if not other[aobj.name].supported:
                        unsupported.append(aobj.name)
                except KeyError:
                    unsupported.append(aobj.name)

        return unsupported

    # TODO: Define standard properties
    # TODO: Link units to UnitSet

    def list_supported_properties(self):
        """
        Return a list of properties supported by this package

        Returns:
            list
        """
        list = []
        for p in self:
            if p.supported:
                list.append(p)
        return list

    def list_unsupported_properties(self):
        """
        Return a list of properties not supported by this package

        Returns:
            list
        """
        list = []
        for p in self:
            if not p.supported:
                list.append(p)
        return list

    def list_required_properties(self):
        """
        Return a list of properties required by this package

        Returns:
            list
        """
        list = []
        for p in self:
            if p.required:
                list.append(p)
        return list

    @property
    def unitset(self):
        """
        UnitSet associated with this PropertySet (via the parent metadata object).
        """
        return self.__parent_block._default_units

    def _define_standard_properties(self):
        # Concrete definition of all the standard IDAES properties
        # TODO: Should we separate thermophysical and reaction properties?
        # AL: I am inclined to say no - define all of them, and state which are supported
        # This would allow for hybrid packages in the future
        self.define_property(
            name="act_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="act_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="act_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="act_coeff_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="act_coeff_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="act_coeff_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="compress_fact",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="compress_fact_phase",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="conc_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="conc_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="conc_mass_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="conc_mass_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="conc_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self.define_property(
            name="conc_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self.define_property(
            name="conc_mol_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self.define_property(
            name="conc_mol_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )

        self.define_property(
            name="cp_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cp_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cp_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cp_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cp_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self.define_property(
            name="cp_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self.define_property(
            name="cp_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self.define_property(
            name="cp_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )

        self.define_property(
            name="cv_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cv_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cv_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cv_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self.define_property(
            name="cv_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self.define_property(
            name="cv_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self.define_property(
            name="cv_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self.define_property(
            name="cv_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )

        self.define_property(
            name="dens_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="dens_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="dens_mass_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="dens_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self.define_property(
            name="dens_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self.define_property(
            name="dens_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self.define_property(
            name="dens_mol_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self.define_property(
            name="dens_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )

        self.define_property(
            name="diffus_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )
        self.define_property(
            name="diffus_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )
        self.define_property(
            name="diffus_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )
        self.define_property(
            name="diffus_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )

        self.define_property(
            name="energy_internal_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="energy_internal_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="energy_internal_mass_phase_como",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="energy_internal_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="energy_internal_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="energy_internal_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )

        self.define_property(
            name="enth_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="enth_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="enth_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="enth_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="enth_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="enth_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="enth_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )

        self.define_property(
            name="entr_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="entr_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="entr_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="entr_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="entr_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="entr_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="entr_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )

        self.define_property(
            name="flow_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self.define_property(
            name="flow_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self.define_property(
            name="flow_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self.define_property(
            name="flow_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self.define_property(
            name="flow_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self.define_property(
            name="flow_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self.define_property(
            name="flow_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self.define_property(
            name="flow_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self.define_property(
            name="flow_vol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )
        self.define_property(
            name="flow_vol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )
        self.define_property(
            name="flow_vol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )
        self.define_property(
            name="flow_vol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )

        self.define_property(
            name="fug_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="fug_coeff_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="heat_capacity_ratio",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="heat_capacity_ratio_phase",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="gibbs_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="gibbs_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="gibbs_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self.define_property(
            name="gibbs_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="gibbs_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="gibbs_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )

        self.define_property(
            name="isentropic_speed_sound_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.VELOCITY,
        )
        self.define_property(
            name="isothermal_speed_sound_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.VELOCITY,
        )

        self.define_property(
            name="henry",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
            # TODO: Units are an issue here, as there are multiple ways to define this
        )

        self.define_property(
            name="mass_frac_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="mass_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="mass_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="mass_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="mole_frac_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="mole_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="mole_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="mole_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="molality_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLALITY,
        )
        self.define_property(
            name="molality_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLALITY,
        )
        self.define_property(
            name="molality_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLALITY,
        )

        self.define_property(
            name="mw",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )
        self.define_property(
            name="mw_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )
        self.define_property(
            name="mw_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )
        self.define_property(
            name="mw_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )

        self.define_property(
            name="phase_frac",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="pressure",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_bubble",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_dew",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_osm_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_red",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="pressure_sat",  # TODO: Deprecate in favour of pressure_sat
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="pressure_sat_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self.define_property(
            name="surf_tens_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.SURFACE_TENSION,
        )

        self.define_property(
            name="temperature",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self.define_property(
            name="temperature_bubble",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self.define_property(
            name="temperature_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self.define_property(
            name="temperature_dew",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self.define_property(
            name="temperature_red",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="temperature_sat",  # TODO: Deprecate in favour of temperature_sat_comp?
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )

        self.define_property(
            name="therm_cond",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.THERMAL_CONDUCTIVITY,
        )
        self.define_property(
            name="therm_cond_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.THERMAL_CONDUCTIVITY,
        )

        self.define_property(
            name="visc_d",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DYNAMIC_VISCOSITY,
        )
        self.define_property(
            name="visc_d_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DYNAMIC_VISCOSITY,
        )
        self.define_property(
            name="visc_k_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.KINEMATIC_VISCOSITY,
        )
        self.define_property(
            name="visc_k",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.KINEMATIC_VISCOSITY,
        )

        self.define_property(
            name="vol_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLAR_VOLUME,
        )

        self.define_property(
            name="vol_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLAR_VOLUME,
        )

        # Log terms
        self.define_property(
            name="log_act_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_act_phase_solvents",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_act_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_act_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_conc_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_conc_mol_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_conc_mol_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="log_mass_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mass_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mass_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="log_molality_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_molality_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_molality_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="log_mole_frac_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_pbub",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_pdew",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_tbub",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_tdew",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_mole_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self.define_property(
            name="log_pressure_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_pressure_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_pressure_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        # Reaction Properties
        # TODO: Units are also problematic here - no single definition
        self.define_property(
            name="dh_rxn",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self.define_property(
            name="k_eq",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="log_k_eq",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="k_rxn",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self.define_property(
            name="reaction_rate",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )


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
        self._default_units = UnitSet()
        self._properties = PropertySet(parent=self)

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

    def add_default_units(self, u):
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
        - 'supported': (optional, default = True) bool indicating if this property is
                       supported by this package.
        - 'required': (optional, default = False) bool indicating if this property is
                      required by this package.

        Args:
            p (dict): Key=property, Value=dict

        Returns:
            None
        """
        for k, v in p.items():
            try:
                getattr(self._properties, k).update_property(v)
            except AttributeError:
                # TODO: Deprecate this and make it raise an exception if an unknown property is encountered
                # # Force users to explicitly declare new/custom properties
                # self._properties.define_property(name=k, **v)
                raise PropertyPackageError(
                    f"Encountered metadata for unknown property {k}. If you are trying to create a custom "
                    "property please use the define_custom_properties method instead."
                )

    def define_custom_properties(self, p):
        """Add custom properties to the metadata.

        For each property, the value should be another dict which may contain
        the following keys:

        - 'method': (required) the name of a method to construct the
                    property as a str, or None if the property will be
                    constructed by default.
        - 'units': (optional) units of measurement for the property.
        - 'supported': (optional, default = True) bool indicating if this property is
                       supported by this package.
        - 'required': (optional, default = False) bool indicating if this property is
                      required by this package.

        Args:
            p (dict): Key=property, Value=dict

        Returns:
            None
        """
        for k, v in p.items():
            self._properties.define_property(name=k, **v)

    def add_required_properties(self, p):
        # TODO: Update doc string
        """Add required properties to the metadata.

        Update 'required' attribute of specified properties.
        Note that argument must be a dict for backwards compatibility.

        Args:
            p (dict): Key=property, Value=(ignored)

        Returns:
            None
        """
        for k, v in p.items():
            try:
                self._properties[k].set_required(True)
            except KeyError:
                self._properties.define_property(name=k, supported=False, required=True)

    def get_derived_units(self, units):
        # TODO: Deprecate in favour of common units property
        return self.derived_units[units]
