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
from copy import copy

from pyomo.environ import units
from pyomo.core.base.units_container import _PyomoUnit

from idaes.core.util.exceptions import PropertyPackageError, BurntToast
import idaes.logger as idaeslog

__author__ = "Dan Gunter <dkgunter@lbl.gov>, Andrew Lee"

_log = idaeslog.getLogger(__name__)


class _PropertyMetadataIndex(object):
    def __init__(
        self,
        parent,
        method: str = None,
        units: _PyomoUnit = None,
        supported: bool = False,
        required: bool = False,
    ):
        super().__setattr__("_parent", parent)
        super().__setattr__("_method", method)
        super().__setattr__("_supported", supported)
        super().__setattr__("_required", required)

        # TODO: For future, this would be the place to store default scaling information, etc.
        # TODO: Could also define default bounds, nominal values, etc.

    def __setattr__(self, key, value):
        raise TypeError("Property metadata does not support assignment.")

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
        return self._parent.units

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

    def set_method(self, meth: str):
        """
        Set method attribute of property.

        Args:
            meth: reference to method required to construct this property

        Returns:
            None
        """
        super().__setattr__("_method", str(meth))

    def set_supported(self, supported: bool = True):
        """
        Set supported attribute of property

        Args:
            supported: bool indicating whether package supports this property

        Returns:
            None
        """
        super().__setattr__("_supported", bool(supported))

    def set_required(self, required: bool = True):
        """
        Set required attribute of property

        Args:
            required: bool indicating whether package requires this property be defined by
            another property package

        Returns:
            None
        """
        super().__setattr__("_required", bool(required))

    def update_property(
        self, method: str = None, required: bool = None, supported: bool = None
    ):
        """
        Update attributes of this property.

        Args:
            method: method name (str) to be assigned to construct property
            required : bool indicating whether this property package requires this property to be
                defined by another package
            supported : bool indicating whether this property package supports this property

        Returns:
            None

        Note that if not provided a value, 'supported' is assumed to be True.
        """
        if method is not None:
            self.set_method(method)
        if required is not None:
            self.set_required(required)
        if supported is not None:
            self.set_supported(supported)
        else:
            # Assume supported if not explicitly stated
            # TODO: Reconsider in the future, for now do this for backwards compatibility
            self.set_supported(True)


class PropertyMetadata(object):
    """
    Metadata object for defining a property.

    This object stores all the metadata associated with a single property, including:

        - name of property
        - documentation of property
        - units of measurement for this property (defined via associated UnitSet)
        - method that constructs this property and associated constraints (if build-on-demand)
        - whether property is supported by this package
        - whether this package expects this property to be provided by another package
    """

    # TODO: Add optional doc string to metadata objects
    def __init__(
        self,
        name: str = None,
        doc: str = None,
        method: str = None,
        units: _PyomoUnit = None,
        supported: bool = False,
        required: bool = False,
        indices: list = None,
    ):
        if name is None:
            raise TypeError('"name" is required')
        if doc is None:
            doc = name

        super().__setattr__("_name", name)
        super().__setattr__("_doc", doc)

        # TODO: Move to indexed elements
        super().__setattr__("_method", method)
        super().__setattr__("_supported", supported)
        super().__setattr__("_required", required)

        # TODO: Validate units are from UnitSet or dimensionless - needs deprecation
        super().__setattr__("_units", units)

        # Create entries for indexed sub-properties
        super().__setattr__("_none", _PropertyMetadataIndex(parent=self))
        if indices is not None:
            for i in indices:
                super().__setattr__("_" + i, _PropertyMetadataIndex(parent=self))

    def __getitem__(self, key: str):
        try:
            return getattr(self, "_" + key)
        except AttributeError:
            # TODO: Real error message - needs to be a KeyError to work with getattr elsewhere
            raise KeyError()

    def __repr__(self):
        return f"{self._doc} ({self._units}%s%s)" % (
            "supported" if self._supported else "",
            "required" if self._required else "",
        )

    def __setattr__(self, key, value):
        raise TypeError("Property metadata does not support assignment.")

    @property
    def name(self):
        """
        Doc string for property
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

    def set_method(self, meth: str):
        """
        Set method attribute of property.

        Args:
            meth: reference to method required to construct this property

        Returns:
            None
        """
        # TODO: Validate that meth is callable?
        super().__setattr__("_method", meth)

    def set_supported(self, supported: bool = True):
        """
        Set supported attribute of property

        Args:
            supported: bool indicating whether package supports this property

        Returns:
            None
        """
        # TODO: Validate that supported is bool
        super().__setattr__("_supported", supported)

    def set_required(self, required: bool = True):
        """
        Set required attribute of property

        Args:
            required: bool indicating whether package requires this property be defined by
            another property package

        Returns:
            None
        """
        # TODO: Validate that required is bool
        super().__setattr__("_required", required)

    def update_property(
        self, method: str = None, required: bool = None, supported: bool = None
    ):
        """
        Update attributes of this property.

        Args:
            method: method name (str) to be assigned to construct property
            required : bool indicating whether this property package requires this property to be
                defined by another package
            supported : bool indicating whether this property package supports this property

        Returns:
            None

        Note that if not provided a value, 'supported' is assumed to be True.
        """
        if method is not None:
            self.set_method(method)
        if required is not None:
            self.set_required(required)
        if supported is not None:
            self.set_supported(supported)
        else:
            # Assume supported if not explicitly stated
            # TODO: Reconsider in the future, for now do this for backwards compatibility
            self.set_supported(True)


class PropertySetBase(object):
    """
    Base class for defining property sets.

    This defines the common methods expected of all PropertySets.
    """

    # Define the standard indices for IDAES properties
    _defined_indices = ["comp", "phase", "phase_comp"]

    def __init__(self, parent):
        super().__setattr__("_parent_block", parent)
        super().__setattr__("_defined_properties", [])
        super().__setattr__("_defined_indices", copy(self.__class__._defined_indices))

        # Find stadanrd properties defined in class and create instance versions
        for i in dir(self.__class__):
            if not i.startswith("_"):
                # Ignore anything starting with "_"
                iobj = getattr(self.__class__, i)
                if isinstance(iobj, PropertyMetadata):
                    units = iobj.units
                    # Check if units is a placeholder string and update if required
                    if isinstance(units, str):
                        units = getattr(self._parent_block.default_units, units)
                    self._add_property(
                        name=i,
                        doc=iobj._doc,
                        method=None,
                        supported=None,
                        required=None,
                        units=units,
                    )

    def __setattr__(self, key, value):
        raise TypeError(
            "PropertySets do not support direct assignment. Please use define_property"
        )

    def __getitem__(self, key: str):
        try:
            return getattr(self, key)
        except AttributeError:
            raise KeyError(f"Property {key} is not defined in this PropertySet.")

    def __iter__(self):
        for a in self._defined_properties:
            yield a

    def define_property(
        self,
        name: str,
        doc: str = None,
        method: str = None,
        required: bool = False,
        supported: bool = True,
        units: _PyomoUnit = None,
    ):
        """
        Define a new property called `name`.

        Args:
            name: name of new property (required)
            doc: doc string for property
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

        self._add_property(
            name=name,
            doc=doc,
            method=method,
            supported=supported,
            required=required,
            units=units,
        )

    def _add_property(
        self,
        name: str,
        doc: str = None,
        method: str = None,
        required: bool = False,
        supported: bool = False,
        units: _PyomoUnit = None,
    ):
        super().__setattr__(
            name,
            PropertyMetadata(
                name=name,
                doc=doc,
                method=method,
                supported=supported,
                required=required,
                units=units,
                indices=self._defined_indices,
            ),
        )
        self._defined_properties.append(getattr(self, name))

    def check_required_properties(self, other: "PropertySetBase"):
        """
        Check that other property package supports all properties marked as required by this package.

        Args:
            other: PropertySet to check for supported properties

        Returns:
            list of properties required by this package which are not supported by other package
        """
        unsupported = []
        for a in self._defined_properties:
            if a.required:
                try:
                    if not other[a.name].supported:
                        unsupported.append(a.name)
                except KeyError:
                    unsupported.append(a.name)

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
        return self._parent_block.default_units

    def get_name_and_index(self, property: str):
        """
        Separates an indexed property name into the main property and index,

        This method is written assuming the standard indices, and checks for 'phase' and then
        'component' before checking for any custom indices. Developers of custom PropertySets
        may wish to overload this with custom search logic better suited to the set of defined
        indices.

        Returns:
            name, index: strings indicating the name of the base property and indexing set.
        """
        name = None
        index = None

        if property in self._defined_properties:
            name = property
        if "phase" in property:
            sep_point = property.index("phase") - 1
            name = property[:sep_point]
            index = property[sep_point + 1 :]
        elif "comp" in property:
            sep_point = property.index("comp") - 1
            name = property[:sep_point]
            index = property[sep_point + 1 :]
        else:
            for i in self._defined_indices:
                if i in property:
                    sep_point = property.index(i) - 1
                    name = property[:sep_point]
                    index = property[sep_point + 1 :]

        if name is None:
            raise ValueError(
                f"Unhandled property: {property}. This is mostly likely due to the "
                "property not being defined in this PropertySet."
            )

        return name, index


class StandardPropertySet(PropertySetBase):
    """
    This object defines all the standard properties supported by IDAES, and also allows for
    definition of new properties if required.
    """

    # Concrete definition of all the standard IDAES properties
    # TODO: Should we separate thermophysical and reaction properties?
    # AL: I am inclined to say no - define all of them, and state which are supported
    # This would allow for hybrid packages in the future
    act_phase_comp = PropertyMetadata(
        name="act_phase_comp",
        units=units.dimensionless,
    )
    act_coeff_phase_comp = PropertyMetadata(
        name="act_coeff_phase_comp",
        units=units.dimensionless,
    )
    compress_fact = PropertyMetadata(
        name="compress_fact",
        units=units.dimensionless,
    )
    compress_fact_phase = PropertyMetadata(
        name="compress_fact_phase",
        units=units.dimensionless,
    )
    conc_mass_comp = PropertyMetadata(
        name="conc_mass_comp",
        units="DENSITY_MASS",
    )
    conc_mass_phase_comp = PropertyMetadata(
        name="conc_mass_phase_comp",
        units="DENSITY_MASS",
    )
    conc_mol_comp = PropertyMetadata(
        name="conc_mol_comp",
        units="DENSITY_MOLE",
    )
    conc_mol_phase_comp = PropertyMetadata(
        name="conc_mol_phase_comp",
        units="DENSITY_MOLE",
    )
    cp_mass = PropertyMetadata(
        name="cp_mass",
        units="HEAT_CAPACITY_MASS",
    )
    cp_mass_comp = PropertyMetadata(
        name="cp_mass_comp",
        units="HEAT_CAPACITY_MASS",
    )
    cp_mass_phase = PropertyMetadata(
        name="cp_mass_phase",
        units="HEAT_CAPACITY_MASS",
    )
    cp_mass_phase_comp = PropertyMetadata(
        name="cp_mass_phase_comp",
        units="HEAT_CAPACITY_MASS",
    )
    cp_mol = PropertyMetadata(
        name="cp_mol",
        units="HEAT_CAPACITY_MOLE",
    )
    cp_mol_comp = PropertyMetadata(
        name="cp_mol_comp",
        units="HEAT_CAPACITY_MOLE",
    )
    cp_mol_phase = PropertyMetadata(
        name="cp_mol_phase",
        units="HEAT_CAPACITY_MOLE",
    )
    cp_mol_phase_comp = PropertyMetadata(
        name="cp_mol_phase_comp",
        units="HEAT_CAPACITY_MOLE",
    )
    cv_mass = PropertyMetadata(
        name="cv_mass",
        units="HEAT_CAPACITY_MASS",
    )
    cv_mass_comp = PropertyMetadata(
        name="cv_mass_comp",
        units="HEAT_CAPACITY_MASS",
    )
    cv_mass_phase = PropertyMetadata(
        name="cv_mass_phase",
        units="HEAT_CAPACITY_MASS",
    )
    cv_mass_phase_comp = PropertyMetadata(
        name="cv_mass_phase_comp",
        units="HEAT_CAPACITY_MASS",
    )
    cv_mol = PropertyMetadata(
        name="cv_mol",
        units="HEAT_CAPACITY_MOLE",
    )
    cv_mol_comp = PropertyMetadata(
        name="cv_mol_comp",
        units="HEAT_CAPACITY_MOLE",
    )
    cv_mol_phase = PropertyMetadata(
        name="cv_mol_phase",
        units="HEAT_CAPACITY_MOLE",
    )
    cv_mol_phase_comp = PropertyMetadata(
        name="cv_mol_phase_comp",
        units="HEAT_CAPACITY_MOLE",
    )
    dens_mass = PropertyMetadata(
        name="dens_mass",
        units="DENSITY_MASS",
    )
    dens_mass_comp = PropertyMetadata(
        name="dens_mass_comp",
        units="DENSITY_MASS",
    )
    dens_mass_crit = PropertyMetadata(
        name="dens_mass_crit",
        units="DENSITY_MASS",
    )
    dens_mass_phase = PropertyMetadata(
        name="dens_mass_phase",
        units="DENSITY_MASS",
    )
    dens_mol = PropertyMetadata(
        name="dens_mol",
        units="DENSITY_MOLE",
    )
    dens_mol_comp = PropertyMetadata(
        name="dens_mol_comp",
        units="DENSITY_MOLE",
    )
    dens_mol_crit = PropertyMetadata(
        name="dens_mol_crit",
        units="DENSITY_MOLE",
    )
    dens_mol_phase = PropertyMetadata(
        name="dens_mol_phase",
        units="DENSITY_MOLE",
    )
    diffus_comp = PropertyMetadata(
        name="diffus_comp",
        units="DIFFUSIVITY",
    )
    diffus_phase_comp = PropertyMetadata(
        name="diffus_phase_comp",
        units="DIFFUSIVITY",
    )
    energy_internal_mass = PropertyMetadata(
        name="energy_internal_mass",
        units="ENERGY_MASS",
    )
    energy_internal_mass_phase = PropertyMetadata(
        name="energy_internal_mass_phase",
        units="ENERGY_MASS",
    )
    energy_internal_mass_phase_comp = PropertyMetadata(
        name="energy_internal_mass_phase_comp",
        units="ENERGY_MASS",
    )
    energy_internal_mol = PropertyMetadata(
        name="energy_internal_mol",
        units="ENERGY_MOLE",
    )
    energy_internal_mol_phase = PropertyMetadata(
        name="energy_internal_mol_phase",
        units="ENERGY_MOLE",
    )
    energy_internal_mol_phase_comp = PropertyMetadata(
        name="energy_internal_mol_phase_comp",
        units="ENERGY_MOLE",
    )
    enth_mass = PropertyMetadata(
        name="enth_mass",
        units="ENERGY_MASS",
    )
    enth_mass_phase = PropertyMetadata(
        name="enth_mass_phase",
        units="ENERGY_MASS",
    )
    enth_mass_phase_comp = PropertyMetadata(
        name="enth_mass_phase_comp",
        units="ENERGY_MASS",
    )
    enth_mol = PropertyMetadata(
        name="enth_mol",
        units="ENERGY_MOLE",
    )
    enth_mol_comp = PropertyMetadata(
        name="enth_mol_comp",
        units="ENERGY_MOLE",
    )
    enth_mol_phase = PropertyMetadata(
        name="enth_mol_phase",
        units="ENERGY_MOLE",
    )
    enth_mol_phase_comp = PropertyMetadata(
        name="enth_mol_phase_comp",
        units="ENERGY_MOLE",
    )
    entr_mass = PropertyMetadata(
        name="entr_mass",
        units="ENERGY_MASS",
    )
    entr_mass_phase = PropertyMetadata(
        name="entr_mass_phase",
        units="ENERGY_MASS",
    )
    entr_mass_phase_comp = PropertyMetadata(
        name="entr_mass_phase_comp",
        units="ENERGY_MASS",
    )
    entr_mol = PropertyMetadata(
        name="entr_mol",
        units="ENERGY_MOLE",
    )
    entr_mol_comp = PropertyMetadata(
        name="entr_mol_comp",
        units="ENERGY_MOLE",
    )
    entr_mol_phase = PropertyMetadata(
        name="entr_mol_phase",
        units="ENERGY_MOLE",
    )
    entr_mol_phase_comp = PropertyMetadata(
        name="entr_mol_phase_comp",
        units="ENERGY_MOLE",
    )
    flow_mass = PropertyMetadata(
        name="flow_mass",
        units="FLOW_MASS",
    )
    flow_mass_comp = PropertyMetadata(
        name="flow_mass_comp",
        units="FLOW_MASS",
    )
    flow_mass_phase = PropertyMetadata(
        name="flow_mass_phase",
        units="FLOW_MASS",
    )
    flow_mass_phase_comp = PropertyMetadata(
        name="flow_mass_phase_comp",
        units="FLOW_MASS",
    )
    flow_mol = PropertyMetadata(
        name="flow_mol",
        units="FLOW_MOLE",
    )
    flow_mol_comp = PropertyMetadata(
        name="flow_mol_comp",
        units="FLOW_MOLE",
    )
    flow_mol_phase = PropertyMetadata(
        name="flow_mol_phase",
        units="FLOW_MOLE",
    )
    flow_mol_phase_comp = PropertyMetadata(
        name="flow_mol_phase_comp",
        units="FLOW_MOLE",
    )
    flow_vol = PropertyMetadata(
        name="flow_vol",
        units="FLOW_VOL",
    )
    flow_vol_comp = PropertyMetadata(
        name="flow_vol_comp",
        units="FLOW_VOL",
    )
    flow_vol_phase = PropertyMetadata(
        name="flow_vol_phase",
        units="FLOW_VOL",
    )
    flow_vol_phase_comp = PropertyMetadata(
        name="flow_vol_phase_comp",
        units="FLOW_VOL",
    )
    fug_phase_comp = PropertyMetadata(
        name="fug_phase_comp",
        units="PRESSURE",
    )
    fug_coeff_phase_comp = PropertyMetadata(
        name="fug_coeff_phase_comp",
        units=units.dimensionless,
    )
    heat_capacity_ratio = PropertyMetadata(
        name="heat_capacity_ratio",
        units=units.dimensionless,
    )
    heat_capacity_ratio_phase = PropertyMetadata(
        name="heat_capacity_ratio_phase",
        units=units.dimensionless,
    )
    gibbs_mass = PropertyMetadata(
        name="gibbs_mass",
        units="ENERGY_MASS",
    )
    gibbs_mass_phase = PropertyMetadata(
        name="gibbs_mass_phase",
        units="ENERGY_MASS",
    )
    gibbs_mass_phase_comp = PropertyMetadata(
        name="gibbs_mass_phase_comp",
        units="ENERGY_MASS",
    )
    gibbs_mol = PropertyMetadata(
        name="gibbs_mol",
        units="ENERGY_MOLE",
    )
    gibbs_mol_phase = PropertyMetadata(
        name="gibbs_mol_phase",
        units="ENERGY_MOLE",
    )
    gibbs_mol_phase_comp = PropertyMetadata(
        name="gibbs_mol_phase_comp",
        units="ENERGY_MOLE",
    )
    isentropic_speed_sound_phase = PropertyMetadata(
        name="isentropic_speed_sound_phase",
        units="VELOCITY",
    )
    isothermal_speed_sound_phase = PropertyMetadata(
        name="isothermal_speed_sound_phase",
        units="VELOCITY",
    )
    henry = PropertyMetadata(
        name="henry",
        units=units.dimensionless,
        # TODO: Units are an issue here, as there are multiple ways to define this
    )
    mass_frac_comp = PropertyMetadata(
        name="mass_frac_comp",
        units=units.dimensionless,
    )
    mass_frac_phase_comp = PropertyMetadata(
        name="mass_frac_phase_comp",
        units=units.dimensionless,
    )
    mole_frac_comp = PropertyMetadata(
        name="mole_frac_comp",
        units=units.dimensionless,
    )
    mole_frac_phase_comp = PropertyMetadata(
        name="mole_frac_phase_comp",
        units=units.dimensionless,
    )
    molality_phase_comp = PropertyMetadata(
        name="molality_phase_comp",
        units="MOLALITY",
    )
    mw = PropertyMetadata(
        name="mw",
        units="MOLECULAR_WEIGHT",
    )
    mw_comp = PropertyMetadata(
        name="mw_comp",
        units="MOLECULAR_WEIGHT",
    )
    mw_phase = PropertyMetadata(
        name="mw_phase",
        units="MOLECULAR_WEIGHT",
    )
    mw_phase_comp = PropertyMetadata(
        name="mw_phase_comp",
        units="MOLECULAR_WEIGHT",
    )
    phase_frac = PropertyMetadata(
        name="phase_frac",
        units=units.dimensionless,
    )
    pressure = PropertyMetadata(
        name="pressure",
        units="PRESSURE",
    )
    pressure_phase_comp = PropertyMetadata(
        name="pressure_phase_comp",
        units="PRESSURE",
    )
    pressure_bubble = PropertyMetadata(
        name="pressure_bubble",
        units="PRESSURE",
    )
    pressure_crit = PropertyMetadata(
        name="pressure_crit",
        units="PRESSURE",
    )
    pressure_dew = PropertyMetadata(
        name="pressure_dew",
        units="PRESSURE",
    )
    pressure_osm_phase = PropertyMetadata(
        name="pressure_osm_phase",
        units="PRESSURE",
    )
    pressure_red = PropertyMetadata(
        name="pressure_red",
        units=units.dimensionless,
    )
    pressure_sat = PropertyMetadata(
        name="pressure_sat",  # TODO: Deprecate in favour of pressure_sat_comp
        units="PRESSURE",
    )
    pressure_sat_comp = PropertyMetadata(
        name="pressure_sat_comp",
        units="PRESSURE",
    )
    surf_tens_phase = PropertyMetadata(
        name="surf_tens_phase",
        units="SURFACE_TENSION",
    )
    temperature = PropertyMetadata(
        name="temperature",
        units="TEMPERATURE",
    )
    temperature_bubble = PropertyMetadata(
        name="temperature_bubble",
        units="TEMPERATURE",
    )
    temperature_crit = PropertyMetadata(
        name="temperature_crit",
        units="TEMPERATURE",
    )
    temperature_dew = PropertyMetadata(
        name="temperature_dew",
        units="TEMPERATURE",
    )
    temperature_red = PropertyMetadata(
        name="temperature_red",
        units=units.dimensionless,
    )
    temperature_sat = PropertyMetadata(
        name="temperature_sat",  # TODO: Deprecate in favour of temperature_sat_comp?
        units="TEMPERATURE",
    )
    therm_cond = PropertyMetadata(
        name="therm_cond",
        units="THERMAL_CONDUCTIVITY",
    )
    therm_cond_phase = PropertyMetadata(
        name="therm_cond_phase",
        units="THERMAL_CONDUCTIVITY",
    )
    visc_d = PropertyMetadata(
        name="visc_d",
        units="DYNAMIC_VISCOSITY",
    )
    visc_d_phase = PropertyMetadata(
        name="visc_d_phase",
        units="DYNAMIC_VISCOSITY",
    )
    visc_k = PropertyMetadata(
        name="visc_k",
        units="KINEMATIC_VISCOSITY",
    )
    visc_k_phase = PropertyMetadata(
        name="visc_k_phase",
        units="KINEMATIC_VISCOSITY",
    )
    vol_mol_phase = PropertyMetadata(
        name="vol_mol_phase",
        units="MOLAR_VOLUME",
    )
    vol_mol_phase_comp = PropertyMetadata(
        name="vol_mol_phase_comp",
        units="MOLAR_VOLUME",
    )
    # Log terms
    log_act_phase_comp = PropertyMetadata(
        name="log_act_phase_comp",
        units=units.dimensionless,
    )
    log_conc_mol_phase_comp = PropertyMetadata(
        name="log_conc_mol_phase_comp",
        units=units.dimensionless,
    )
    log_mass_frac_phase_comp = PropertyMetadata(
        name="log_mass_frac_phase_comp",
        units=units.dimensionless,
    )
    log_molality_phase_comp = PropertyMetadata(
        name="log_molality_phase_comp",
        units=units.dimensionless,
    )
    log_mole_frac_comp = PropertyMetadata(
        name="log_mole_frac_comp",
        units=units.dimensionless,
    )
    log_mole_frac_pbub = PropertyMetadata(
        name="log_mole_frac_pbub",
        units=units.dimensionless,
    )
    log_mole_frac_pdew = PropertyMetadata(
        name="log_mole_frac_pdew",
        units=units.dimensionless,
    )
    log_mole_frac_tbub = PropertyMetadata(
        name="log_mole_frac_tbub",
        units=units.dimensionless,
    )
    log_mole_frac_tdew = PropertyMetadata(
        name="log_mole_frac_tdew",
        units=units.dimensionless,
    )
    log_mole_frac_phase_comp = PropertyMetadata(
        name="log_mole_frac_phase_comp",
        units=units.dimensionless,
    )
    log_pressure_phase_comp = PropertyMetadata(
        name="log_pressure_phase_comp",
        units=units.dimensionless,
    )

    # Reaction Properties
    # TODO: Units are also problematic here - no single definition
    dh_rxn = PropertyMetadata(
        name="dh_rxn",
        units="ENERGY_MOLE",
    )
    k_eq = PropertyMetadata(
        name="k_eq",
        units=units.dimensionless,
    )
    log_k_eq = PropertyMetadata(
        name="log_k_eq",
        units=units.dimensionless,
    )
    k_rxn = PropertyMetadata(
        name="k_rxn",
        units=units.dimensionless,
    )
    reaction_rate = PropertyMetadata(
        name="reaction_rate",
        units=units.dimensionless,
    )


class ElectrolytePropertySet(StandardPropertySet):
    """
    This object defines all the standard properties supported by IDAES, and also allows for
    definition of new properties if required.
    """

    _defined_indices = [
        "comp",
        "phase",
        "phase_comp",
        "phase_comp_apparent",
        "phase_comp_true",
    ]

    # Definition of additional properties required for electrolyte applications
    act_phase_comp_apparent = PropertyMetadata(
        name="act_phase_comp_apparent",
        units=units.dimensionless,
    )
    act_phase_comp_true = PropertyMetadata(
        name="act_phase_comp_true",
        units=units.dimensionless,
    )
    act_coeff_phase_comp_apparent = PropertyMetadata(
        name="act_coeff_phase_comp_apparent",
        units=units.dimensionless,
    )
    act_coeff_phase_comp_true = PropertyMetadata(
        name="act_coeff_phase_comp_true",
        units=units.dimensionless,
    )
    conc_mass_phase_comp_apparent = PropertyMetadata(
        name="conc_mass_phase_comp_apparent",
        units="DENSITY_MASS",
    )
    conc_mass_phase_comp_true = PropertyMetadata(
        name="conc_mass_phase_comp_true",
        units="DENSITY_MASS",
    )
    conc_mol_phase_comp_apparent = PropertyMetadata(
        name="conc_mol_phase_comp_apparent",
        units="DENSITY_MOLE",
    )
    conc_mol_phase_comp_true = PropertyMetadata(
        name="conc_mol_phase_comp_true",
        units="DENSITY_MOLE",
    )
    diffus_phase_comp_apparent = PropertyMetadata(
        name="diffus_phase_comp_apparent",
        units="DIFFUSIVITY",
    )
    diffus_phase_comp_true = PropertyMetadata(
        name="diffus_phase_comp_true",
        units="DIFFUSIVITY",
    )
    mass_frac_phase_comp_apparent = PropertyMetadata(
        name="mass_frac_phase_comp_apparent",
        units=units.dimensionless,
    )
    mass_frac_phase_comp_true = PropertyMetadata(
        name="mass_frac_phase_comp_true",
        units=units.dimensionless,
    )
    mole_frac_phase_comp_apparent = PropertyMetadata(
        name="mole_frac_phase_comp_apparent",
        units=units.dimensionless,
    )
    mole_frac_phase_comp_true = PropertyMetadata(
        name="mole_frac_phase_comp_true",
        units=units.dimensionless,
    )
    molality_phase_comp_apparent = PropertyMetadata(
        name="molality_phase_comp_apparent",
        units="MOLALITY",
    )
    molality_phase_comp_true = PropertyMetadata(
        name="molality_phase_comp_true",
        units="MOLALITY",
    )
    pressure_phase_comp_apparent = PropertyMetadata(
        name="pressure_phase_comp_apparent",
        units="PRESSURE",
    )
    pressure_phase_comp_true = PropertyMetadata(
        name="pressure_phase_comp_true",
        units="PRESSURE",
    )
    # Log terms
    log_act_phase_solvents = PropertyMetadata(
        name="log_act_phase_solvents",
        units=units.dimensionless,
    )
    log_act_phase_comp_apparent = PropertyMetadata(
        name="log_act_phase_comp_apparent",
        units=units.dimensionless,
    )
    log_act_phase_comp_true = PropertyMetadata(
        name="log_act_phase_comp_true",
        units=units.dimensionless,
    )
    log_conc_mol_phase_comp_apparent = PropertyMetadata(
        name="log_conc_mol_phase_comp_apparent",
        units=units.dimensionless,
    )
    log_conc_mol_phase_comp_true = PropertyMetadata(
        name="log_conc_mol_phase_comp_true",
        units=units.dimensionless,
    )
    log_mass_frac_phase_comp_apparent = PropertyMetadata(
        name="log_mass_frac_phase_comp_apparent",
        units=units.dimensionless,
    )
    log_mass_frac_phase_comp_true = PropertyMetadata(
        name="log_mass_frac_phase_comp_true",
        units=units.dimensionless,
    )
    log_molality_phase_comp_apparent = PropertyMetadata(
        name="log_molality_phase_comp_apparent",
        units=units.dimensionless,
    )
    log_molality_phase_comp_true = PropertyMetadata(
        name="log_molality_phase_comp_true",
        units=units.dimensionless,
    )
    log_mole_frac_phase_comp_apparent = PropertyMetadata(
        name="log_mole_frac_phase_comp_apparent",
        units=units.dimensionless,
    )
    log_mole_frac_phase_comp_true = PropertyMetadata(
        name="log_mole_frac_phase_comp_true",
        units=units.dimensionless,
    )
    log_pressure_phase_comp_apparent = PropertyMetadata(
        name="log_pressure_phase_comp_apparent",
        units=units.dimensionless,
    )
    log_pressure_phase_comp_true = PropertyMetadata(
        name="log_pressure_phase_comp_true",
        units=units.dimensionless,
    )
