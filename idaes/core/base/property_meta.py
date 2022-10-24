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


class PropertyMetadata(object):
    """Container for property parameter metadata.

    Instances of this class are exactly dictionaries, with the
    only difference being some guidance on the values expected in the
    dictionary from the constructor.
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
        return self._name

    @property
    def method(self):
        return self._method

    @property
    def units(self):
        return self._units

    @property
    def supported(self):
        return self._supported

    @property
    def required(self):
        return self._required

    def set_method(self, meth):
        # TODO: Validate that meth is callable?
        self._method = meth

    def set_supported(self, supported=True):
        # TODO: Validate that supported is bool
        self._supported = supported

    def set_required(self, required=True):
        # TODO: Validate that required is bool
        self._required = required

    def update_property(self, dict):
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
    # TODO: Add doc string

    def __init__(self, parent):
        self.__parent_block = parent

        self._define_standard_properties()

    def __getitem__(self, key):
        try:
            return getattr(self, "_" + key)
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
        # Method to define new, custom properties
        setattr(
            self,
            "_" + name,
            PropertyMetadata(
                name=name,
                method=method,
                supported=supported,
                required=required,
                units=units,
            ),
        )

    def check_required_properties(self, other):
        # Check that other package supports properties this package requires
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

    @property
    def unitset(self):
        return self.__parent_block._default_units

    def _define_standard_properties(self):
        # Concrete definition of all the standard IDAES properties
        # TODO: Should we separate thermophysical and reaction properties?
        # AL: I am inclined to say no - define all of them, and state which are supported
        # This would allow for hybrid packages in the future
        self.define_property(
            name="flow_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        # "flow_mol": {"method": "_flow_mol"},
        # "flow_vol": {"method": "_flow_vol"},
        # "flow_mass": {"method": "_flow_mass"},
        # "flow_mass_phase": {"method": "_flow_mass_phase"},
        # "flow_vol_phase": {"method": "_flow_vol_phase"},
        # "flow_mol_phase": {"method": "_flow_mol_phase"},
        # "flow_mass_comp": {"method": "_flow_mass_comp"},
        # "flow_mol_comp": {"method": "_flow_mol_comp"},
        # "flow_mass_phase_comp": {"method": "_flow_mass_phase_comp"},
        # "flow_mol_phase_comp": {"method": "_flow_mol_phase_comp"},
        # "mole_frac_comp": {"method": "_mole_frac_comp"},
        # "mole_frac_phase_comp": {"method": None},
        # "phase_frac": {"method": None},
        # "temperature": {"method": None},
        # "pressure": {"method": None},
        # "act_phase_comp": {"method": "_act_phase_comp"},
        # "act_phase_comp_true": {"method": "_act_phase_comp_true"},
        # "act_phase_comp_apparent": {"method": "_act_phase_comp_apparent"},
        # "act_coeff_phase_comp": {"method": "_act_coeff_phase_comp"},
        # "act_coeff_phase_comp_true": {"method": "_act_coeff_phase_comp_true"},
        # "act_coeff_phase_comp_apparent": {
        #     "method": "_act_coeff_phase_comp_apparent"
        # },
        # "compress_fact_phase": {"method": "_compress_fact_phase"},
        # "conc_mol_comp": {"method": "_conc_mol_comp"},
        # "conc_mol_phase_comp": {"method": "_conc_mol_phase_comp"},
        # "conc_mol_phase_comp_apparent": {
        #     "method": "_conc_mol_phase_comp_apparent"
        # },
        # "conc_mol_phase_comp_true": {"method": "_conc_mol_phase_comp_true"},
        # "cp_mol": {"method": "_cp_mol"},
        # "cp_mol_phase": {"method": "_cp_mol_phase"},
        # "cp_mol_phase_comp": {"method": "_cp_mol_phase_comp"},
        # "cv_mol": {"method": "_cv_mol"},
        # "cv_mol_phase": {"method": "_cv_mol_phase"},
        # "cv_mol_phase_comp": {"method": "_cv_mol_phase_comp"},
        # "diffus_phase_comp": {"method": "_diffus_phase_comp"},
        # "diffus_phase_comp_apparent": {"method": "_diffus_phase_comp_apparent"},
        # "diffus_phase_comp_true": {"method": "_diffus_phase_comp_true"},
        # "heat_capacity_ratio_phase": {"method": "_heat_capacity_ratio_phase"},
        # "dens_mass": {"method": "_dens_mass"},
        # "dens_mass_phase": {"method": "_dens_mass_phase"},
        # "dens_mol": {"method": "_dens_mol"},
        # "dens_mol_phase": {"method": "_dens_mol_phase"},
        # "energy_internal_mol": {"method": "_energy_internal_mol"},
        # "energy_internal_mol_phase": {"method": "_energy_internal_mol_phase"},
        # "energy_internal_mol_phase_comp": {
        #     "method": "_energy_internal_mol_phase_comp"
        # },
        # "enth_mol": {"method": "_enth_mol"},
        # "enth_mol_phase": {"method": "_enth_mol_phase"},
        # "enth_mol_phase_comp": {"method": "_enth_mol_phase_comp"},
        # "entr_mol": {"method": "_entr_mol"},
        # "entr_mol_phase": {"method": "_entr_mol_phase"},
        # "entr_mol_phase_comp": {"method": "_entr_mol_phase_comp"},
        # "fug_phase_comp": {"method": "_fug_phase_comp"},
        # "fug_coeff_phase_comp": {"method": "_fug_coeff_phase_comp"},
        # "gibbs_mol": {"method": "_gibbs_mol"},
        # "gibbs_mol_phase": {"method": "_gibbs_mol_phase"},
        # "gibbs_mol_phase_comp": {"method": "_gibbs_mol_phase_comp"},
        # "isentropic_speed_sound_phase": {
        #     "method": "_isentropic_speed_sound_phase"
        # },
        # "isothermal_speed_sound_phase": {
        #     "method": "_isothermal_speed_sound_phase"
        # },
        # "henry": {"method": "_henry"},
        # "mass_frac_phase_comp": {"method": "_mass_frac_phase_comp"},
        # "mass_frac_phase_comp_apparent": {
        #     "method": "_mass_frac_phase_comp_apparent"
        # },
        # "mass_frac_phase_comp_true": {"method": "_mass_frac_phase_comp_true"},
        # "molality_phase_comp": {"method": "_molality_phase_comp"},
        # "molality_phase_comp_apparent": {
        #     "method": "_molality_phase_comp_apparent"
        # },
        # "molality_phase_comp_true": {"method": "_molality_phase_comp_true"},
        # "mw": {"method": "_mw"},
        # "mw_comp": {"method": "_mw_comp"},
        # "mw_phase": {"method": "_mw_phase"},
        # "pressure_phase_comp": {"method": "_pressure_phase_comp"},
        # "pressure_phase_comp_true": {"method": "_pressure_phase_comp_true"},
        # "pressure_phase_comp_apparent": {
        #     "method": "_pressure_phase_comp_apparent"
        # },
        # "pressure_bubble": {"method": "_pressure_bubble"},
        # "pressure_dew": {"method": "_pressure_dew"},
        # "pressure_osm_phase": {"method": "_pressure_osm_phase"},
        # "pressure_sat_comp": {"method": "_pressure_sat_comp"},
        # "surf_tens_phase": {"method": "_surf_tens_phase"},
        # "temperature_bubble": {"method": "_temperature_bubble"},
        # "temperature_dew": {"method": "_temperature_dew"},
        # "therm_cond_phase": {"method": "_therm_cond_phase"},
        # "visc_d_phase": {"method": "_visc_d_phase"},
        # "vol_mol_phase": {"method": "_vol_mol_phase"},
        # "vol_mol_phase_comp": {"method": "_vol_mol_phase_comp"},
        # "dh_rxn": {"method": "_dh_rxn"},
        # "log_act_phase_comp": {"method": "_log_act_phase_comp"},
        # "log_act_phase_solvents": {"method": "_log_act_phase_solvents"},
        # "log_act_phase_comp_true": {"method": "_log_act_phase_comp_true"},
        # "log_act_phase_comp_apparent": {
        #     "method": "_log_act_phase_comp_apparent"
        # },
        # "log_conc_mol_phase_comp": {"method": "_log_conc_mol_phase_comp"},
        # "log_conc_mol_phase_comp_true": {
        #     "method": "_log_conc_mol_phase_comp_true"
        # },
        # "log_mass_frac_phase_comp": {"method": "_log_mass_frac_phase_comp"},
        # "log_mass_frac_phase_comp_apparent": {
        #     "method": "_log_mass_frac_phase_comp_apparent"
        # },
        # "log_mass_frac_phase_comp_true": {
        #     "method": "_log_mass_frac_phase_comp_true"
        # },
        # "log_molality_phase_comp": {"method": "_log_molality_phase_comp"},
        # "log_molality_phase_comp_apparent": {
        #     "method": "_log_molality_phase_comp_apparent"
        # },
        # "log_molality_phase_comp_true": {
        #     "method": "_log_molality_phase_comp_true"
        # },
        # "log_mole_frac_comp": {"method": "_log_mole_frac_comp"},
        # "log_mole_frac_tbub": {"method": "_log_mole_frac_tbub"},
        # "log_mole_frac_tdew": {"method": "_log_mole_frac_tdew"},
        # "log_mole_frac_pbub": {"method": "_log_mole_frac_pbub"},
        # "log_mole_frac_pdew": {"method": "_log_mole_frac_pdew"},
        # "log_mole_frac_phase_comp": {"method": "_log_mole_frac_phase_comp"},
        # "log_mole_frac_phase_comp_apparent": {
        #     "method": "_log_mole_frac_phase_comp_apparent"
        # },
        # "log_mole_frac_phase_comp_true": {
        #     "method": "_log_mole_frac_phase_comp_true"
        # },
        # "log_pressure_phase_comp": {"method": "_log_pressure_phase_comp"},
        # "log_pressure_phase_comp_apparent": {
        #     "method": "_log_pressure_phase_comp_apparent"
        # },
        # "log_pressure_phase_comp_true": {
        #     "method": "_log_pressure_phase_comp_true"
        # },
        # "log_k_eq": {"method": "_log_k_eq"},


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

        Args:
            p (dict): Key=property, Value=PropertyMetadata or equiv. dict

        Returns:
            None
        """
        for k, v in p.items():
            try:
                prop = getattr(self._properties, k).update_property(**v)
            except AttributeError:
                # TODO: Deprecate this and make it raise an exception if an unknown property is encountered?
                # Force users to explicitly declare new/custom properties
                self._properties.define_property(name=k, **v)

    def add_required_properties(self, p):
        # TODO: Update doc string
        """Add required properties to the metadata.

        For each property, the value should be the expected units of
        measurement for the property.

        Args:
            p (dict): Key=property, Value=units

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
