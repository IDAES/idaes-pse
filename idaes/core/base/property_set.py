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
This module defines the sets of properties used by IDAES.

There is a base set of standard properties which are applicable in most applications,
and it is intended that specialty applications can and will define their own
property sets as required (e.g. electrolyte systems).
"""
from pyomo.environ import units

from idaes.core.util.exceptions import PropertyPackageError
import idaes.logger as idaeslog

__author__ = "Dan Gunter <dkgunter@lbl.gov>, Andrew Lee"

_log = idaeslog.getLogger(__name__)


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


class PropertySetBase(object):
    """
    Base class for defining property sets.

    This defines the common methods expected of all PropertySets.
    """

    def __init__(self, parent):
        self.__parent_block = parent

    def __getitem__(self, key):
        try:
            return getattr(self, key.upper())
        except AttributeError:
            try:
                return getattr(self, key)
            except AttributeError:
                raise KeyError(f"Property {key} is not defined in this PropertySet.")

    def __iter__(self):
        for a in dir(self):
            if not a.startswith("_"):
                # Ignore anything starting with "_" as they will show up twice otherwise
                # This will also skip private attributes of object class
                aobj = getattr(self, a)
                if isinstance(aobj, PropertyMetadata):
                    yield aobj

    def _add_property_meta(
        self,
        local_name=None,
        name=None,
        method=None,
        supported=True,
        required=False,
        units=None,
    ):
        """
        Internal method for creating property metadata objects
        """
        if hasattr(self, name):
            raise PropertyPackageError(
                f"A property with the name {name} already exists. Please use update_property "
                "method if you wish to update an existing property's metadata."
            )

        setattr(
            self,
            local_name,
            PropertyMetadata(
                name=name,
                method=method,
                supported=supported,
                required=required,
                units=units,
            ),
        )

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
        self._add_property_meta(
            local_name=name,
            name=name,
            method=method,
            supported=supported,
            required=required,
            units=units,
        )

    def _define_standard_property(
        self, name=None, method=None, supported=True, required=False, units=None
    ):
        """
        Private method for defining standard properties (for which a @property will also be declared).

        This sets local_name = "_"+name
        """
        self._add_property_meta(
            local_name="_" + name,
            name=name,
            method=method,
            supported=supported,
            required=required,
            units=units,
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
        Reference to UnitSet associated with this PropertySet (via the parent metadata object).
        """
        return self.__parent_block._default_units


class StandardPropertySet(PropertySetBase):
    """
    This object defines all the standard properties supported by IDAES, and also allows for
    definition of new properties if required.
    """

    def __init__(self, parent):
        super().__init__(parent)

        # Concrete definition of all the standard IDAES properties
        # TODO: Should we separate thermophysical and reaction properties?
        # AL: I am inclined to say no - define all of them, and state which are supported
        # This would allow for hybrid packages in the future
        self._define_standard_property(
            name="act_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="act_coeff_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="compress_fact",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="compress_fact_phase",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="conc_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="conc_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="conc_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="conc_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="cp_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cp_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cp_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cp_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cp_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cp_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cp_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cp_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cv_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cv_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cv_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cv_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MASS,
        )
        self._define_standard_property(
            name="cv_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cv_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cv_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="cv_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.HEAT_CAPACITY_MOLE,
        )
        self._define_standard_property(
            name="dens_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="dens_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="dens_mass_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="dens_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="dens_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="dens_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="dens_mol_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="dens_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="diffus_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )
        self._define_standard_property(
            name="diffus_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )
        self._define_standard_property(
            name="energy_internal_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="energy_internal_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="energy_internal_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="energy_internal_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="energy_internal_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="energy_internal_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="enth_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="enth_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="enth_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="enth_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="enth_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="enth_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="enth_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="entr_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="entr_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="entr_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="entr_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="entr_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="entr_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="entr_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="flow_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self._define_standard_property(
            name="flow_mass_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self._define_standard_property(
            name="flow_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self._define_standard_property(
            name="flow_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MASS,
        )
        self._define_standard_property(
            name="flow_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self._define_standard_property(
            name="flow_mol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self._define_standard_property(
            name="flow_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self._define_standard_property(
            name="flow_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_MOLE,
        )
        self._define_standard_property(
            name="flow_vol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )
        self._define_standard_property(
            name="flow_vol_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )
        self._define_standard_property(
            name="flow_vol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )
        self._define_standard_property(
            name="flow_vol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.FLOW_VOL,
        )

        self._define_standard_property(
            name="fug_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="fug_coeff_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="heat_capacity_ratio",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="heat_capacity_ratio_phase",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="gibbs_mass",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="gibbs_mass_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="gibbs_mass_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MASS,
        )
        self._define_standard_property(
            name="gibbs_mol",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="gibbs_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="gibbs_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )

        self._define_standard_property(
            name="isentropic_speed_sound_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.VELOCITY,
        )
        self._define_standard_property(
            name="isothermal_speed_sound_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.VELOCITY,
        )

        self._define_standard_property(
            name="henry",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
            # TODO: Units are an issue here, as there are multiple ways to define this
        )

        self._define_standard_property(
            name="mass_frac_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="mass_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="mole_frac_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="mole_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="molality_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLALITY,
        )
        self._define_standard_property(
            name="mw",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )
        self._define_standard_property(
            name="mw_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )
        self._define_standard_property(
            name="mw_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )
        self._define_standard_property(
            name="mw_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLECULAR_WEIGHT,
        )

        self._define_standard_property(
            name="phase_frac",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="pressure",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_bubble",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_dew",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_osm_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_red",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="pressure_sat",  # TODO: Deprecate in favour of pressure_sat
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_sat_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="surf_tens_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.SURFACE_TENSION,
        )

        self._define_standard_property(
            name="temperature",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self._define_standard_property(
            name="temperature_bubble",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self._define_standard_property(
            name="temperature_crit",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self._define_standard_property(
            name="temperature_dew",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )
        self._define_standard_property(
            name="temperature_red",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="temperature_sat",  # TODO: Deprecate in favour of temperature_sat_comp?
            method=None,
            supported=False,
            required=False,
            units=self.unitset.TEMPERATURE,
        )

        self._define_standard_property(
            name="therm_cond",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.THERMAL_CONDUCTIVITY,
        )
        self._define_standard_property(
            name="therm_cond_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.THERMAL_CONDUCTIVITY,
        )

        self._define_standard_property(
            name="visc_d",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DYNAMIC_VISCOSITY,
        )
        self._define_standard_property(
            name="visc_d_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DYNAMIC_VISCOSITY,
        )
        self._define_standard_property(
            name="visc_k_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.KINEMATIC_VISCOSITY,
        )
        self._define_standard_property(
            name="visc_k",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.KINEMATIC_VISCOSITY,
        )

        self._define_standard_property(
            name="vol_mol_phase",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLAR_VOLUME,
        )

        self._define_standard_property(
            name="vol_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLAR_VOLUME,
        )

        # Log terms
        self._define_standard_property(
            name="log_act_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_conc_mol_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mass_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_molality_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_pbub",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_pdew",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_tbub",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_tdew",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_pressure_phase_comp",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        # Reaction Properties
        # TODO: Units are also problematic here - no single definition
        self._define_standard_property(
            name="dh_rxn",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.ENERGY_MOLE,
        )
        self._define_standard_property(
            name="k_eq",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_k_eq",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="k_rxn",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="reaction_rate",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

    @property
    def ACT_PHASE_COMP(self):
        return self._act_phase_comp

    @property
    def ACT_COEFF_PHASE_COMP(self):
        return self._act_coeff_phase_comp

    @property
    def COMPRESS_FACT(self):
        return self._compress_fact

    @property
    def COMPRESS_FACT_PHASE(self):
        return self._compress_fact_phase

    @property
    def CONC_MASS_COMP(self):
        return self._conc_mass_comp

    @property
    def CONC_MASS_PHASE_COMP(self):
        return self._conc_mass_phase_comp

    @property
    def CONC_MOL_COMP(self):
        return self._conc_mol_comp

    @property
    def CONC_MOL_PHASE_COMP(self):
        return self._conc_mol_phase_comp

    @property
    def CP_MASS(self):
        return self._cp_mass

    @property
    def CP_MASS_COMP(self):
        return self._cp_mass_comp

    @property
    def CP_MASS_PHASE(self):
        return self._cp_mass_phase

    @property
    def CP_MASS_PHASE_COMP(self):
        return self._cp_mass_phase_comp

    @property
    def CP_MOL(self):
        return self._cp_mol

    @property
    def CP_MOL_COMP(self):
        return self._cp_mol_comp

    @property
    def CP_MOL_PHASE(self):
        return self._cp_mol_phase

    @property
    def CP_MOL_PHASE_COMP(self):
        return self._cp_mol_phase_comp

    @property
    def CV_MASS(self):
        return self._cv_mass

    @property
    def CV_MASS_COMP(self):
        return self._cv_mass_comp

    @property
    def CV_MASS_PHASE(self):
        return self._cv_mass_phase

    @property
    def CV_MASS_PHASE_COMP(self):
        return self._cv_mass_phase_comp

    @property
    def CV_MOL(self):
        return self._cv_mol

    @property
    def CV_MOL_COMP(self):
        return self._cv_mol_comp

    @property
    def CV_MOL_PHASE(self):
        return self._cv_mol_phase

    @property
    def CV_MOL_PHASE_COMP(self):
        return self._cv_mol_phase_comp

    @property
    def DENS_MASS(self):
        return self._dens_mass

    @property
    def DENS_MASS_COMP(self):
        return self._dens_mass_comp

    @property
    def DENS_MASS_CRIT(self):
        return self._dens_mass_crit

    @property
    def DENS_MASS_PHASE(self):
        return self._dens_mass_phase

    @property
    def DENS_MOL(self):
        return self._dens_mol

    @property
    def DENS_MOL_COMP(self):
        return self._dens_mol_comp

    @property
    def DENS_MOL_CRIT(self):
        return self._dens_mol_crit

    @property
    def DENS_MOL_PHASE(self):
        return self._dens_mol_phase

    @property
    def DIFFUS_COMP(self):
        return self._diffus_comp

    @property
    def DIFFUS_PHASE_COMP(self):
        return self._diffus_phase_comp

    @property
    def ENERGY_INTERNAL_MASS(self):
        return self._energy_internal_mass

    @property
    def ENERGY_INTERNAL_MASS_PHASE(self):
        return self._energy_internal_mass_phase

    @property
    def ENERGY_INTERNAL_MASS_PHASE_COMP(self):
        return self._energy_internal_mass_phase_comp

    @property
    def ENERGY_INTERNAL_MOL(self):
        return self._energy_internal_mol

    @property
    def ENERGY_INTERNAL_MOL_PHASE(self):
        return self._energy_internal_mol_phase

    @property
    def ENERGY_INTERNAL_MOL_PHASE_COMP(self):
        return self._energy_internal_mol_phase_comp

    @property
    def ENTH_MASS(self):
        return self._enth_mass

    @property
    def ENTH_MASS_PHASE(self):
        return self._enth_mass_phase

    @property
    def ENTH_MASS_PHASE_COMP(self):
        return self._enth_mass_phase_comp

    @property
    def ENTH_MOL(self):
        return self._enth_mol

    @property
    def ENTH_MOL_COMP(self):
        return self._enth_mol_comp

    @property
    def ENTH_MOL_PHASE(self):
        return self._enth_mol_phase

    @property
    def ENTH_MOL_PHASE_COMP(self):
        return self._enth_mol_phase_comp

    @property
    def ENTR_MASS(self):
        return self._entr_mass

    @property
    def ENTR_MASS_PHASE(self):
        return self._entr_mass_phase

    @property
    def ENTR_MASS_PHASE_COMP(self):
        return self._entr_mass_phase_comp

    @property
    def ENTR_MOL(self):
        return self._entr_mol

    @property
    def ENTR_MOL_COMP(self):
        return self._entr_mol_comp

    @property
    def ENTR_MOL_PHASE(self):
        return self._entr_mol_phase

    @property
    def ENTR_MOL_PHASE_COMP(self):
        return self._entr_mol_phase_comp

    @property
    def FLOW_MASS(self):
        return self._flow_mass

    @property
    def FLOW_MASS_COMP(self):
        return self._flow_mass_comp

    @property
    def FLOW_MASS_PHASE(self):
        return self._flow_mass_phase

    @property
    def FLOW_MASS_PHASE_COMP(self):
        return self._flow_mass_phase_comp

    @property
    def FLOW_MOL(self):
        return self._flow_mol

    @property
    def FLOW_MOL_COMP(self):
        return self._flow_mol_comp

    @property
    def FLOW_MOL_PHASE(self):
        return self._flow_mol_phase

    @property
    def FLOW_MOL_PHASE_COMP(self):
        return self._flow_mol_phase_comp

    @property
    def FLOW_VOL(self):
        return self._flow_vol

    @property
    def FLOW_VOL_COMP(self):
        return self._flow_vol_comp

    @property
    def FLOW_VOL_PHASE(self):
        return self._flow_vol_phase

    @property
    def FLOW_VOL_PHASE_COMP(self):
        return self._flow_vol_phase_comp

    @property
    def FUG_PHASE_COMP(self):
        return self._fug_phase_comp

    @property
    def FUG_COEFF_PHASE_COMP(self):
        return self._fug_coeff_phase_comp

    @property
    def HEAT_CAPACITY_RATIO(self):
        return self._heat_capacity_ratio

    @property
    def HEAT_CAPACITY_RATIO_PHASE(self):
        return self._heat_capacity_ratio_phase

    @property
    def GIBBS_MASS(self):
        return self._gibbs_mass

    @property
    def GIBBS_MASS_PHASE(self):
        return self._gibbs_mass_phase

    @property
    def GIBBS_MASS_PHASE_COMP(self):
        return self._gibbs_mass_phase_comp

    @property
    def GIBBS_MOL(self):
        return self._gibbs_mol

    @property
    def GIBBS_MOL_PHASE(self):
        return self._gibbs_mol_phase

    @property
    def GIBBS_MOL_PHASE_COMP(self):
        return self._gibbs_mol_phase_comp

    @property
    def ISENTROPIC_SPEED_SOUND_PHASE(self):
        return self._isentropic_speed_sound_phase

    @property
    def ISOTHERMAL_SPEED_SOUND_PHASE(self):
        return self._isothermal_speed_sound_phase

    @property
    def HENRY(self):
        return self._henry  # TODO: Units issue

    @property
    def MASS_FRAC_COMP(self):
        return self._mass_frac_comp

    @property
    def MASS_FRAC_PHASE_COMP(self):
        return self._mass_frac_phase_comp

    @property
    def MOLE_FRAC_COMP(self):
        return self._mole_frac_comp

    @property
    def MOLE_FRAC_PHASE_COMP(self):
        return self._mole_frac_phase_comp

    @property
    def MOLALITY_PHASE_COMP(self):
        return self._molality_phase_comp

    @property
    def MW(self):
        return self._mw

    @property
    def MW_COMP(self):
        return self._mw_comp

    @property
    def MW_PHASE(self):
        return self._mw_phase

    @property
    def MW_PHASE_COMP(self):
        return self._mw_phase_comp  # TODO: THis should be redundant with _comp

    @property
    def PHASE_FRAC(self):
        return self._phase_frac

    @property
    def PRESSURE(self):
        return self._pressure

    @property
    def PRESSURE_PHASE_COMP(self):
        return self._pressure_phase_comp

    @property
    def PRESSURE_BUBBLE(self):
        return self._pressure_bubble

    @property
    def PRESSURE_CRIT(self):
        return self._pressure_crit

    @property
    def PRESSURE_DEW(self):
        return self._pressure_dew

    @property
    def PRESSURE_OSM_PHASE(self):
        return self._pressure_osm_phase

    @property
    def PRESSURE_RED(self):
        return self._pressure_red

    @property
    def PRESSURE_SAT(self):
        return self._pressure_sat  # TODO: Deprecate for _comp

    @property
    def PRESSURE_SAT_COMP(self):
        return self._pressure_sat_comp

    @property
    def SURF_TENS_PHASE(self):
        return self._surf_tens_phase

    @property
    def TEMPERATURE(self):
        return self._temperature

    @property
    def TEMPERATURE_BUBBLE(self):
        return self._temperature_bubble

    @property
    def TEMPERATURE_CRIT(self):
        return self._temperature_crit

    @property
    def TEMPERATURE_DEW(self):
        return self._temperature_dew

    @property
    def TEMPERATURE_RED(self):
        return self._temperature_red

    @property
    def TEMPERATURE_SAT(self):
        return self._temperature_sat  # TODO: Deprecate for _comp

    @property
    def THERM_COND(self):
        return self._therm_cond

    @property
    def THERM_COND_PHASE(self):
        return self._therm_cond_phase

    @property
    def VISC_D(self):
        return self._visc_d

    @property
    def VISC_D_PHASE(self):
        return self._visc_d_phase

    @property
    def VISC_K_PHASE(self):
        return self._visc_k_phase

    @property
    def VISC_K(self):
        return self._visc_k

    @property
    def VOL_MOL_PHASE(self):
        return self._vol_mol_phase

    @property
    def VOL_MOL_PHASE_COMP(self):
        return self._vol_mol_phase_comp

    @property
    def LOG_ACT_PHASE_COMP(self):
        return self._log_act_phase_comp

    @property
    def LOG_CONC_MOL_PHASE_COMP(self):
        return self._log_conc_mol_phase_comp

    @property
    def LOG_MASS_FRAC_PHASE_COMP(self):
        return self._log_mass_frac_phase_comp

    @property
    def LOG_MOLALITY_PHASE_COMP(self):
        return self._log_molality_phase_comp

    @property
    def LOG_MOLE_FRAC_COMP(self):
        return self._log_mole_frac_comp

    @property
    def LOG_MOLE_FRAC_PBUB(self):
        return self._log_mole_frac_pbub

    @property
    def LOG_MOLE_FRAC_PDEW(self):
        return self._log_mole_frac_pdew

    @property
    def LOG_MOLE_FRAC_TBUB(self):
        return self._log_mole_frac_tbub

    @property
    def LOG_MOLE_FRAC_TDEW(self):
        return self._log_mole_frac_tdew

    @property
    def LOG_MOLE_FRAC_PHASE_COMP(self):
        return self._log_mole_frac_phase_comp

    @property
    def LOG_PRESSURE_PHASE_COMP(self):
        return self._log_pressure_phase_comp

    @property
    def DH_RXN(self):
        return self._dh_rxn

    @property
    def K_EQ(self):
        return self._k_eq

    @property
    def LOG_K_EQ(self):
        return self._log_k_eq

    @property
    def K_RXN(self):
        return self._k_rxn

    @property
    def REACTION_RATE(self):
        return self._reaction_rate


class ElectrolytePropertySet(StandardPropertySet):
    """
    This object defines all the standard properties supported by IDAES, and also allows for
    definition of new properties if required.
    """

    def __init__(self, parent):
        super().__init__(parent)

        # Definition of additional properties required for electrolyte applications
        self._define_standard_property(
            name="act_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="act_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="act_coeff_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="act_coeff_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="conc_mass_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="conc_mass_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MASS,
        )
        self._define_standard_property(
            name="conc_mol_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )
        self._define_standard_property(
            name="conc_mol_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DENSITY_MOLE,
        )

        self._define_standard_property(
            name="diffus_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )
        self._define_standard_property(
            name="diffus_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.DIFFUSIVITY,
        )

        self._define_standard_property(
            name="mass_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="mass_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="mole_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="mole_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="molality_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLALITY,
        )
        self._define_standard_property(
            name="molality_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.MOLALITY,
        )

        self._define_standard_property(
            name="pressure_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )
        self._define_standard_property(
            name="pressure_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=self.unitset.PRESSURE,
        )

        # Log terms
        self._define_standard_property(
            name="log_act_phase_solvents",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_act_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_act_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_conc_mol_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_conc_mol_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="log_mass_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mass_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="log_molality_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_molality_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="log_mole_frac_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_mole_frac_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

        self._define_standard_property(
            name="log_pressure_phase_comp_apparent",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )
        self._define_standard_property(
            name="log_pressure_phase_comp_true",
            method=None,
            supported=False,
            required=False,
            units=units.dimensionless,
        )

    @property
    def ACT_PHASE_COMP_APPARENT(self):
        return self._act_phase_comp_apparent

    @property
    def ACT_PHASE_COMP_TRUE(self):
        return self._act_phase_comp_true

    @property
    def ACT_COEFF_PHASE_COMP_APPARENT(self):
        return self._act_coeff_phase_comp_apparent

    @property
    def ACT_COEFF_PHASE_COMP_TRUE(self):
        return self._act_coeff_phase_comp_true

    @property
    def CONC_MASS_PHASE_COMP_APPARENT(self):
        return self._conc_mass_phase_comp_apparent

    @property
    def CONC_MASS_PHASE_COMP_TRUE(self):
        return self._conc_mass_phase_comp_true

    @property
    def CONC_MOL_PHASE_COMP_APPARENT(self):
        return self._conc_mol_phase_comp_apparent

    @property
    def CONC_MOL_PHASE_COMP_TRUE(self):
        return self._conc_mol_phase_comp_true

    @property
    def DIFFUS_PHASE_COMP_APPARENT(self):
        return self._diffus_phase_comp_apparent

    @property
    def DIFFUS_PHASE_COMP_TRUE(self):
        return self._diffus_phase_comp_true

    @property
    def MASS_FRAC_PHASE_COMP_APPARENT(self):
        return self._mass_frac_phase_comp_apparent

    @property
    def MASS_FRAC_PHASE_COMP_TRUE(self):
        return self._mass_frac_phase_comp_true

    @property
    def MOLE_FRAC_PHASE_COMP_APPARENT(self):
        return self._mole_frac_phase_comp_apparent

    @property
    def MOLE_FRAC_PHASE_COMP_TRUE(self):
        return self._mole_frac_phase_comp_true

    @property
    def MOLALITY_PHASE_COMP_APPARENT(self):
        return self._molality_phase_comp_apparent

    @property
    def MOLALITY_PHASE_COMP_TRUE(self):
        return self._molality_phase_comp_true

    @property
    def PRESSURE_PHASE_COMP_APPARENT(self):
        return self._pressure_phase_comp_apparent

    @property
    def PRESSURE_PHASE_COMP_TRUE(self):
        return self._pressure_phase_comp_true

    @property
    def LOG_ACT_PHASE_SOLVENTS(self):
        return self._log_act_phase_solvents

    @property
    def LOG_ACT_PHASE_COMP_APPARENT(self):
        return self._log_act_phase_comp_apparent

    @property
    def LOG_ACT_PHASE_COMP_TRUE(self):
        return self._log_act_phase_comp_true

    @property
    def LOG_CONC_MOL_PHASE_COMP_APPARENT(self):
        return self._log_conc_mol_phase_comp_apparent

    @property
    def LOG_CONC_MOL_PHASE_COMP_TRUE(self):
        return self._log_conc_mol_phase_comp_true

    @property
    def LOG_MASS_FRAC_PHASE_COMP_APPARENT(self):
        return self._log_mass_frac_phase_comp_apparent

    @property
    def LOG_MASS_FRAC_PHASE_COMP_TRUE(self):
        return self._log_mass_frac_phase_comp_true

    @property
    def LOG_MOLALITY_PHASE_COMP_APPARENT(self):
        return self._log_molality_phase_comp_apparent

    @property
    def LOG_MOLALITY_PHASE_COMP_TRUE(self):
        return self._log_molality_phase_comp_true

    @property
    def LOG_MOLE_FRAC_PHASE_COMP_APPARENT(self):
        return self._log_mole_frac_phase_comp_apparent

    @property
    def LOG_MOLE_FRAC_PHASE_COMP_TRUE(self):
        return self._log_mole_frac_phase_comp_true

    @property
    def LOG_PRESSURE_PHASE_COMP(self):
        return self._log_pressure_phase_comp

    @property
    def LOG_PRESSURE_PHASE_COMP_APPARENT(self):
        return self._log_pressure_phase_comp_apparent

    @property
    def LOG_PRESSURE_PHASE_COMP_TRUE(self):
        return self._log_pressure_phase_comp_true
