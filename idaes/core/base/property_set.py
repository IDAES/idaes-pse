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
This module defines the sets of properties used by IDAES.

There is a base set of standard properties which are applicable in most applications,
and it is intended that specialty applications can and will define their own
property sets as required (e.g. electrolyte systems).
"""
from copy import copy

from pyomo.environ import units as pyunits
from pyomo.core.base.units_container import _PyomoUnit

from idaes.core.util.exceptions import PropertyPackageError
import idaes.logger as idaeslog

__author__ = "Dan Gunter <dkgunter@lbl.gov>, Andrew Lee"

_log = idaeslog.getLogger(__name__)


class _PropertyMetadataIndex:
    def __init__(
        self,
        parent,
        idx,
        method: str = None,
        supported: bool = False,
        required: bool = False,
        valid_range=None,
    ):
        self._parent = parent
        self._idx = idx
        self._method = method
        self._supported = supported
        self._required = required
        self._valid_range = None
        self._set_valid_range(
            valid_range
        )  # This method does some basic validation of input
        self._lock_setattr = True
        # TODO: For future, this would be the place to store default scaling information, etc.
        # TODO: Could also define default bounds, nominal values, etc.

    def __setattr__(self, key, value):
        if hasattr(self, "_lock_setattr") and self._lock_setattr is True:
            raise TypeError("Property metadata does not support assignment.")
        super().__setattr__(key, value)

    def __repr__(self):
        return f"{self._parent._doc} ({self._parent._units}%s%s)" % (
            "supported" if self._supported else "",
            "required" if self._required else "",
        )

    @property
    def name(self):
        """
        Name of sub-property as a string
        """
        suffix = ""
        if self._idx != "none":
            suffix = f"_{self._idx}"
        return f"{self._parent.name}{suffix}"

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

    @property
    def valid_range(self):
        """
        Tuple indicating valid range of values for this property based on the data and method used
        to calculate it. This should take the form of (lower range, upper range) and can be used to
        set and verify bounds on properties before and after solving in order to inform users of any
        cases where the property values are being extrapolated.

        It is strongly recommended that developers set valid ranges for all state variables in
        a property package, as well as any other key properties.
        """
        return self._valid_range

    def set_method(self, meth: str):
        """
        Set method attribute of property.

        Args:
            meth: reference to method required to construct this property

        Returns:
            None
        """
        if meth is None:
            super().__setattr__("_method", None)
        else:
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

    def _set_valid_range(self, valid_range: tuple):
        """
        Set valid_range attribute of property. WARNING: users should generally not change
        the valid_range of a property, as this is determined by the data and method used
        to fit the property correlation.

        Args:
            valid_range: 2-tuple with valid range of values for property (lower, upper)

        Returns:
            None
        """
        if valid_range is not None:
            # Do some basic verification
            if not isinstance(valid_range, tuple) or not len(valid_range) == 2:
                raise ValueError(
                    f"valid_range must be a tuple of length 2 (got {valid_range})"
                )
            elif valid_range[0] > valid_range[1]:
                raise ValueError(
                    f"valid_range must be a 2-tuple with form (lower, upper): first value "
                    f"was greater than second value: {valid_range}"
                )

        super().__setattr__("_valid_range", valid_range)

    def update_property(
        self,
        method: str = None,
        required: bool = None,
        supported: bool = None,
        valid_range: tuple = None,
    ):
        """
        Update attributes of this property.

        Args:
            method: method name (str) to be assigned to construct property
            required : bool indicating whether this property package requires this property to be
                defined by another package
            supported : bool indicating whether this property package supports this property
            valid_range: 2-tuple with valid range of values for property (lower, upper)

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
        if valid_range is not None:
            self._set_valid_range(valid_range)


class PropertyMetadata:
    """
    Metadata object for defining a property.

    This object stores all the metadata associated with a single property, including:

        - name of property
        - documentation of property
        - units of measurement for this property (defined via associated UnitSet)
        - indexing set for sub-properties (e.g. by component or by phase)
        - attributes for sub-properties identifying

            - method that constructs this property and associated constraints (if build-on-demand)
            - whether property is supported by this package
            - whether this package expects this property to be provided by another package

    """

    # TODO: Add optional doc string to metadata objects
    def __init__(
        self,
        name: str = None,
        doc: str = None,
        units: _PyomoUnit = None,
        indices: list = None,
        initialize: dict = None,
    ):
        """
        Construct a PropertyMetadata object.

        Args:
            name: name of property
            doc: doc string for property
            units: units of measurement for property
            indices: list of sub-property indices
            initialize: dict of metadata attributes to set by sub-property


        Supported metadata attributes are:

            - method: str representing a method name to call to construct this property
            - required: bool indicating whether a property package requires this property to be
                        supplied by another property package
            - supported: bool indicating whether a property package supports this property

        Returns:
            None
        """
        if name is None:
            raise TypeError('"name" is required')
        if doc is None:
            doc = name

        self._name = name
        self._doc = doc

        # TODO: Validate units are from UnitSet or dimensionless - needs deprecation
        self._units = units

        # Record indices - needed for building instance specific cases
        self._indices = indices

        self._lock_setattr = True

        # Create entries for indexed sub-properties
        if indices is None or indices is False:
            indices = ["none"]

        for i in indices:
            super().__setattr__("_" + i, _PropertyMetadataIndex(parent=self, idx=i))

        # Initialize values in sub-properties
        if initialize is not None:
            for k, v in initialize.items():
                self[k].update_property(**v)

    def __getitem__(self, key: str):
        if key is None:
            key = "none"
        try:
            return getattr(self, "_" + key)
        except AttributeError:
            raise KeyError(f"Property {self._name} does not have index {key}.")

    def __repr__(self):
        return f"{self._doc} ({self._units})"

    def __setattr__(self, key, value):
        if hasattr(self, "_lock_setattr") and self._lock_setattr is True:
            raise TypeError("Property metadata does not support assignment.")
        super().__setattr__(key, value)

    @property
    def name(self):
        """
        Doc string for property
        """
        return self._name

    @property
    def units(self):
        """
        Units of measurement for this property. This should be a reference to a quantity defined
        in the UnitSet associated with this property package.
        """
        return self._units

    # TODO: An overall update method across multiple indices?


class PropertySetBase:
    """
    Base class for defining property sets.

    This defines the common methods expected of all PropertySets.
    """

    # Define the standard indices for IDAES properties
    _defined_indices = ["comp", "phase", "phase_comp"]

    def __init__(self, parent):
        self._parent_block = parent
        self._defined_properties = []
        self._defined_indices = copy(self.__class__._defined_indices)
        self._lock_setattr = True
        # Find standard properties defined in class and create instance versions
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
                        units=units,
                        indices=iobj._indices,  # get desired indices from parent
                    )

    def __setattr__(self, key, value):
        if hasattr(self, "_lock_setattr") and self._lock_setattr is True:
            raise TypeError(
                "PropertySets do not support direct assignment. Please use define_property"
            )
        super().__setattr__(key, value)

    def __getitem__(self, key: str):
        n, i = self.get_name_and_index(key)
        try:
            return getattr(self, n)[i]
        except AttributeError:
            raise KeyError(f"Property {key} is not defined in this PropertySet.")

    def __iter__(self):
        for a in self._defined_properties:
            yield getattr(self, a)

    def define_property(
        self,
        name: str,
        doc: str = None,
        units: _PyomoUnit = None,
        indices=None,
        method: str = None,
        required: bool = False,
        supported: bool = True,
        initialize: dict = None,
    ):
        """
        Define a new property called `name`.

        Args:
            name: name of new property (required)
            doc: doc string for property
            units: quantity defined in associated UnitSet defining the units of measurement for property
            indices: list, bool or None. Indicates what sub-property indices should be added. None = use default set,
                     False = unindexed (only 'none' entry)
            method: reference to build-on-demand method for property (optional,
                    only if indices is None or False, default=None)
            supported: bool indicating if package supports this property (optional,
                       only if indices is None or False, default=True)
            required: bool indicating if package requires this property from another package (optional,
                      only if indices is None or False, default=False)
            initialize: dict containing values for method, required and supported by sub-property indices
                        (optional). Cannot combine with method, required or supported arguments.

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
            indices=indices,
            initialize=initialize,
        )

    def _add_property(
        self,
        name: str,
        doc: str = None,
        method: str = None,
        required: bool = False,
        supported: bool = False,
        units: _PyomoUnit = None,
        indices: list = None,
        initialize: dict = None,
    ):
        if indices is not None:
            # Cannot assign method, required or supported as we don't know what index they apply to.
            # Raise exception if these are not default values.
            if required or supported or method is not None:
                raise ValueError(
                    "Cannot assign method, required or supported attributes directly to indexed "
                    "properties. Please use the initialize argument instead."
                )
        else:
            # Otherwise, take global arguments and set to initialize "none" index
            if initialize is None:
                initialize = {
                    "none": {
                        "method": method,
                        "required": required,
                        "supported": supported,
                    },
                }
            else:
                if required or supported or method is not None:
                    raise ValueError(
                        "Cannot provide values for initialize and any of method, required or supported "
                        "arguments. Please use the one approach or the other."
                    )

        if indices is None:
            indices = list(self._defined_indices)
            indices.append("none")
        elif indices is False:
            indices = ["none"]

        super().__setattr__(
            name,
            PropertyMetadata(
                name=name,
                doc=doc,
                units=units,
                indices=indices,
                initialize=initialize,
            ),
        )
        self._defined_properties.append(name)

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
            aobj = getattr(self, a)
            # TODO: Should refactor parent so this is not private
            # pylint: disable-next=protected-access
            for i in aobj._indices:
                if aobj[i].required:
                    try:
                        if not getattr(other, a)[i].supported:
                            unsupported.append(a + "_" + i)
                    except AttributeError:
                        unsupported.append(a + "_" + i)

        return unsupported

    def list_supported_properties(self):
        """
        Return a list of properties supported by this package

        Returns:
            list
        """
        slist = []
        for p in self:
            # TODO: Should refactor parent so this is not private
            # pylint: disable-next=protected-access
            for i in p._indices:
                if p[i].supported:
                    slist.append(p[i])
        return slist

    def list_unsupported_properties(self):
        """
        Return a list of properties not supported by this package

        Returns:
            list
        """
        ulist = []
        for p in self:
            # TODO: Should refactor parent so this is not private
            # pylint: disable-next=protected-access
            for i in p._indices:
                if not p[i].supported:
                    ulist.append(p[i])
        return ulist

    def list_required_properties(self):
        """
        Return a list of properties required by this package

        Returns:
            list
        """
        rlist = []
        for p in self:
            # TODO: Should refactor parent so this is not private
            # pylint: disable-next=protected-access
            for i in p._indices:
                if p[i].required:
                    rlist.append(p[i])
        return rlist

    @property
    def unitset(self):
        """
        Reference to UnitSet associated with this PropertySet (via the parent metadata object).
        """
        return self._parent_block.default_units

    def get_name_and_index(self, property_name: str):
        """
        Separates an indexed property name into the main property and index,

        This method is written assuming the standard indices, and checks for 'phase' and then
        'component' before checking for any custom indices. Developers of custom PropertySets
        may wish to overload this with custom search logic better suited to the set of defined
        indices.

        Returns:
            name, index: strings indicating the name of the base property and indexing set.
        """
        root_name = property_name
        index_name = None
        _defined_indices = {"phase_comp", "phase", "comp"} | set(self._defined_indices)

        _defined_indices = list(reversed(sorted(_defined_indices, key=len)))

        if property_name in self._defined_properties:
            return property_name, None
        else:
            for i in _defined_indices:
                if property_name.endswith("_" + i):
                    nchar = len(i) + 1
                    root_name = property_name[:-nchar]
                    index_name = i
                    break

        if root_name not in self._defined_properties:
            raise ValueError(
                f"Unhandled property: {property_name}. This is mostly likely due to the property not being defined in this PropertySet."
            )

        return root_name, index_name


class StandardPropertySet(PropertySetBase):
    """
    This object defines all the standard properties supported by IDAES, and also allows for
    definition of new properties if required.
    """

    # Concrete definition of all the standard IDAES properties
    # TODO: Should we separate thermophysical and reaction properties?
    # AL: I am inclined to say no - define all of them, and state which are supported
    # This would allow for hybrid packages in the future
    act = PropertyMetadata(
        name="act",
        doc="Chemical Activity",
        units=pyunits.dimensionless,
    )
    act_coeff = PropertyMetadata(
        name="act_coeff",
        doc="Chemical Activity Coefficient",
        units=pyunits.dimensionless,
    )
    compress_fact = PropertyMetadata(
        name="compress_fact",
        doc="Compressiblity Factor",
        units=pyunits.dimensionless,
    )
    conc_mass = PropertyMetadata(
        name="conc_mass",
        doc="Concentration on a Mass Basis",
        units="DENSITY_MASS",
    )
    conc_mol = PropertyMetadata(
        name="conc_mol",
        doc="Concentration on a Molar Basis",
        units="DENSITY_MOLE",
    )
    cp_mass = PropertyMetadata(
        name="cp_mass",
        doc="Constant Pressure Specific Heat Capacity (Mass Basis)",
        units="HEAT_CAPACITY_MASS",
    )
    cp_mol = PropertyMetadata(
        name="cp_mol",
        doc="Constant Pressure Specific Heat Capacity (Molar Basis)",
        units="HEAT_CAPACITY_MOLE",
    )
    cv_mass = PropertyMetadata(
        name="cv_mass",
        doc="Constant Volume Specific Heat Capacity (Mass Basis)",
        units="HEAT_CAPACITY_MASS",
    )
    cv_mol = PropertyMetadata(
        name="cv_mol",
        doc="Constant Volume Specific Heat Capacity (Molar Basis)",
        units="HEAT_CAPACITY_MOLE",
    )
    dens_mass = PropertyMetadata(
        name="dens_mass",
        doc="Density (Mass Basis)",
        units="DENSITY_MASS",
    )
    dens_mass_crit = PropertyMetadata(
        name="dens_mass_crit",
        doc="Mass Density at Critical Point",
        units="DENSITY_MASS",
    )
    dens_mol = PropertyMetadata(
        name="dens_mol",
        doc="Density (Molar Basis)",
        units="DENSITY_MOLE",
    )
    dens_mol_crit = PropertyMetadata(
        name="dens_mol_crit",
        doc="Molar Density at Critical Point",
        units="DENSITY_MOLE",
    )
    diffus = PropertyMetadata(
        name="diffus",
        doc="Diffusivity Coefficient",
        units="DIFFUSIVITY",
    )
    energy_internal_mass = PropertyMetadata(
        name="energy_internal_mass",
        doc="Specific Internal Energy (Mass Basis)",
        units="ENERGY_MASS",
    )
    energy_internal_mol = PropertyMetadata(
        name="energy_internal_mol",
        doc="Specific Internal Energy (Molar Basis)",
        units="ENERGY_MOLE",
    )
    enth_mass = PropertyMetadata(
        name="enth_mass",
        doc="Specific Enthalpy (Mass Basis)",
        units="ENERGY_MASS",
    )
    enth_mol = PropertyMetadata(
        name="enth_mol",
        doc="Specific Enthalpy (Molar Basis)",
        units="ENERGY_MOLE",
    )
    entr_mass = PropertyMetadata(
        name="entr_mass",
        doc="Specific Entropy (Mass Basis)",
        units="ENERGY_MASS",
    )
    entr_mol = PropertyMetadata(
        name="entr_mol",
        doc="Specific Entropy (Molar Basis)",
        units="ENERGY_MOLE",
    )
    flow_mass = PropertyMetadata(
        name="flow_mass",
        doc="Mass Flow Rate",
        units="FLOW_MASS",
    )
    flow_mol = PropertyMetadata(
        name="flow_mol",
        doc="Molar Flow Rate",
        units="FLOW_MOLE",
    )
    flow_vol = PropertyMetadata(
        name="flow_vol",
        doc="Volumetric Flow Rate",
        units="FLOW_VOL",
    )
    fug = PropertyMetadata(
        name="fug",
        doc="Fugacity",
        units="PRESSURE",
    )
    fug_coeff = PropertyMetadata(
        name="fug_coeff",
        doc="Fugacity Coefficient",
        units=pyunits.dimensionless,
    )
    heat_capacity_ratio = PropertyMetadata(
        name="heat_capacity_ratio",
        doc="Heat Capacity Ration",
        units=pyunits.dimensionless,
    )
    gibbs_mass = PropertyMetadata(
        name="gibbs_mass",
        doc="Specific Gibbs Energy (Mass Basis)",
        units="ENERGY_MASS",
    )
    gibbs_mol = PropertyMetadata(
        name="gibbs_mol",
        doc="Specific Gibbs Energy (Molar Basis)",
        units="ENERGY_MOLE",
    )
    isentropic_speed_sound = PropertyMetadata(
        name="isentropic_speed_sound",
        doc="Isentropic Speed of Sound",
        units="VELOCITY",
    )
    isothermal_speed_sound = PropertyMetadata(
        name="isothermal_speed_sound",
        doc="Isothermal Speed of Sound",
        units="VELOCITY",
    )
    henry = PropertyMetadata(
        name="henry",
        doc="Henry Constant",
        units=pyunits.dimensionless,
        # TODO: Units are an issue here, as there are multiple ways to define this
    )
    mass_frac = PropertyMetadata(
        name="mass_frac",
        doc="Mass Fraction",
        units=pyunits.dimensionless,
    )
    mole_frac = PropertyMetadata(
        name="mole_frac",
        doc="Mole Fraction",
        units=pyunits.dimensionless,
    )
    molality = PropertyMetadata(
        name="molality",
        doc="Molality",
        units="MOLALITY",
    )
    mw = PropertyMetadata(
        name="mw",
        doc="Molecular Weight",
        units="MOLECULAR_WEIGHT",
    )
    phase_frac = PropertyMetadata(
        name="phase_frac",
        doc="Phase Fraction",
        units=pyunits.dimensionless,
    )
    prandtl_number = PropertyMetadata(
        name="prandtl_number",
        doc="Prandtl Number",
        units=pyunits.dimensionless,
    )
    pressure = PropertyMetadata(
        name="pressure",
        doc="Pressure",
        units="PRESSURE",
    )
    pressure_bubble = PropertyMetadata(
        name="pressure_bubble",
        doc="Bubble Point Pressure",
        units="PRESSURE",
    )
    pressure_crit = PropertyMetadata(
        name="pressure_crit",
        doc="Pressure at Critical Point",
        units="PRESSURE",
    )
    pressure_dew = PropertyMetadata(
        name="pressure_dew",
        doc="Dew point Pressure",
        units="PRESSURE",
    )
    pressure_osm = PropertyMetadata(
        name="pressure_osm",
        doc="Osmotic Pressure",
        units="PRESSURE",
    )
    pressure_red = PropertyMetadata(
        name="pressure_red",
        doc="Reduced Pressure",
        units=pyunits.dimensionless,
    )
    pressure_sat = PropertyMetadata(
        name="pressure_sat",
        doc="Saturation Pressure",
        units="PRESSURE",
    )
    surf_tens = PropertyMetadata(
        name="surf_tens",
        doc="Surface Tension",
        units="SURFACE_TENSION",
    )
    temperature = PropertyMetadata(
        name="temperature",
        doc="Temperature",
        units="TEMPERATURE",
        indices=False,
    )
    temperature_bubble = PropertyMetadata(
        name="temperature_bubble",
        doc="Bubble Point Temperature",
        units="TEMPERATURE",
        indices=False,
    )
    temperature_crit = PropertyMetadata(
        name="temperature_crit",
        doc="Temperature at Critical Point",
        units="TEMPERATURE",
    )
    temperature_dew = PropertyMetadata(
        name="temperature_dew",
        doc="Dew Point Temperature",
        units="TEMPERATURE",
        indices=False,
    )
    temperature_red = PropertyMetadata(
        name="temperature_red",
        doc="Reduced Temperature",
        units=pyunits.dimensionless,
        indices=False,
    )
    temperature_sat = PropertyMetadata(
        name="temperature_sat",
        doc="Saturation Temperature",
        units="TEMPERATURE",
    )
    therm_cond = PropertyMetadata(
        name="therm_cond",
        doc="Thermal Conductivity",
        units="THERMAL_CONDUCTIVITY",
    )
    visc_d = PropertyMetadata(
        name="visc_d",
        doc="Dynamic Viscosity",
        units="DYNAMIC_VISCOSITY",
    )
    visc_k = PropertyMetadata(
        name="visc_k",
        doc="Kinematic Viscosity",
        units="KINEMATIC_VISCOSITY",
    )
    vol_mol = PropertyMetadata(
        name="vol_mol",
        doc="Molar Volume",
        units="MOLAR_VOLUME",
    )
    vol_mol_crit = PropertyMetadata(
        name="vol_mol",
        doc="Molar Volume at Critical Point",
        units="MOLAR_VOLUME",
    )
    # Log terms
    log_act = PropertyMetadata(
        name="log_act",
        doc="Log of Activity",
        units=pyunits.dimensionless,
    )
    log_act_coeff = PropertyMetadata(
        name="log_act_coeff",
        doc="Log of Activity Coefficient",
        units=pyunits.dimensionless,
    )
    log_conc_mol = PropertyMetadata(
        name="log_conc_mol",
        doc="Log of Molar Concentration",
        units=pyunits.dimensionless,
    )
    log_mass_frac = PropertyMetadata(
        name="log_mass_frac",
        doc="Log of Mass Fractions",
        units=pyunits.dimensionless,
    )
    log_molality = PropertyMetadata(
        name="log_molality",
        doc="Log of Molality",
        units=pyunits.dimensionless,
    )
    log_mole_frac = PropertyMetadata(
        name="log_mole_frac",
        doc="Log of Mole Fractions",
        units=pyunits.dimensionless,
    )
    log_mole_frac_pbub = PropertyMetadata(
        name="log_mole_frac_pbub",
        doc="Log of Mole Fractions at Bubble Point Pressure",
        units=pyunits.dimensionless,
    )
    log_mole_frac_pdew = PropertyMetadata(
        name="log_mole_frac_pdew",
        doc="Log of Mole Fractions at Dew Point Pressure",
        units=pyunits.dimensionless,
    )
    log_mole_frac_tbub = PropertyMetadata(
        name="log_mole_frac_tbub",
        doc="Log of Mole Fractions at Bubble Point Temperature",
        units=pyunits.dimensionless,
    )
    log_mole_frac_tdew = PropertyMetadata(
        name="log_mole_frac_tdew",
        doc="Log of Mole Fractions at Dew Point Temperature",
        units=pyunits.dimensionless,
    )
    log_pressure = PropertyMetadata(
        name="log_pressure",
        doc="Log of Pressure",
        units=pyunits.dimensionless,
    )

    # Reaction Properties
    # TODO: Units are also problematic here - no single definition
    dh_rxn = PropertyMetadata(
        name="dh_rxn",
        doc="Heat of Reaction",
        units="ENERGY_MOLE",
    )
    k_eq = PropertyMetadata(
        name="k_eq",
        doc="Equilibrium Coefficient",
        units=pyunits.dimensionless,
    )
    log_k_eq = PropertyMetadata(
        name="log_k_eq",
        doc="Log of Equilibrium Coefficient",
        units=pyunits.dimensionless,
    )
    k_rxn = PropertyMetadata(
        name="k_rxn",
        doc="Rate Constant",
        units=pyunits.dimensionless,
    )
    reaction_rate = PropertyMetadata(
        name="reaction_rate",
        doc="Reaction Rate",
        units=pyunits.dimensionless,
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
    # Log terms
    log_act_phase_solvents = PropertyMetadata(
        name="log_act_phase_solvents",
        doc="Log of Activity of Solvents",
        units=pyunits.dimensionless,
    )
    ionic_strength = PropertyMetadata(
        name="ionic_strength",
        doc="Ionic Strength (Molality Basis)",
        units="MOLALITY",
    )
    pH = PropertyMetadata(
        name="pH",
        doc="pH, -log10(a_H+)",
        units=pyunits.dimensionless,
    )
    pK = PropertyMetadata(
        name="pK",
        doc="pK, -log10(equilibrium constant)",
        units=pyunits.dimensionless,
    )
    pOH = PropertyMetadata(
        name="pOH",
        doc="pOH, -log10(a_OH-)",
        units=pyunits.dimensionless,
    )
    log10_act_coeff = PropertyMetadata(
        name="log10_act_coeff",
        doc="Log of Activity Coefficient",
        units=pyunits.dimensionless,
    )
    log10_molality = PropertyMetadata(
        name="log10_molality",
        doc="Log of Molality",
        units=pyunits.dimensionless,
    )
    log10_k_eq = PropertyMetadata(
        name="log10_k_eq",
        doc="Log of equilibrium coefficient",
        units=pyunits.dimensionless,
    )
    saturation_index = PropertyMetadata(
        name="saturation_index",
        doc="Saturation Index = log(IAP/Ksp)",
        units=pyunits.dimensionless,
    )
