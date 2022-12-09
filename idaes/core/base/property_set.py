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
        idx,
        method: str = None,
        units: _PyomoUnit = None,
        supported: bool = False,
        required: bool = False,
    ):
        super().__setattr__("_parent", parent)
        super().__setattr__("_idx", idx)
        super().__setattr__("_method", method)
        super().__setattr__("_supported", supported)
        super().__setattr__("_required", required)

        # TODO: For future, this would be the place to store default scaling information, etc.
        # TODO: Could also define default bounds, nominal values, etc.

    def __setattr__(self, key, value):
        raise TypeError("Property metadata does not support assignment.")

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
        if self._idx is not "none":
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
        if name is None:
            raise TypeError('"name" is required')
        if doc is None:
            doc = name

        super().__setattr__("_name", name)
        super().__setattr__("_doc", doc)

        # TODO: Validate units are from UnitSet or dimensionless - needs deprecation
        super().__setattr__("_units", units)

        # Record indices - needed for building instance specific cases
        super().__setattr__("_indices", indices)

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
        raise TypeError("Property metadata does not support assignment.")

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
        raise TypeError(
            "PropertySets do not support direct assignment. Please use define_property"
        )

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
        method: str = None,
        required: bool = False,
        supported: bool = True,
        units: _PyomoUnit = None,
        indices=None,
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
            indices: list, bool or None. Indicates what sub-property indices should be added. None = use default set,
                     False = unindexed (only 'none' entry)/

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
    ):
        # TODO: Need to support defining custom indexed properties and setting metadata
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
            ),
        )
        self._defined_properties.append(name)

        # Set metadata on _none index
        if "none" in indices:
            getattr(self, name)._none.update_property(
                method=method,
                supported=supported,
                required=required,
            )

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
        list = []
        for p in self:
            for i in p._indices:
                if p[i].supported:
                    list.append(p[i])
        return list

    def list_unsupported_properties(self):
        """
        Return a list of properties not supported by this package

        Returns:
            list
        """
        list = []
        for p in self:
            for i in p._indices:
                if not p[i].supported:
                    list.append(p[i])
        return list

    def list_required_properties(self):
        """
        Return a list of properties required by this package

        Returns:
            list
        """
        list = []
        for p in self:
            for i in p._indices:
                if p[i].required:
                    list.append(p[i])
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

        sep_point = None
        if property in self._defined_properties:
            sep_point = 0
        elif "phase" in property and not property.startswith("phase"):
            if "phase_comp" in property:
                sep_point = property.rindex("phase_comp")
            else:
                sep_point = property.rindex("phase")
        elif "comp" in property and not property.startswith("comp"):
            if "phase_comp" in property:
                sep_point = property.rindex("phase_comp")
            else:
                sep_point = property.rindex("comp")
        else:
            for i in self._defined_indices:
                if i in property:
                    sep_point = property.rindex(i)

        if sep_point is None:
            raise ValueError(
                f"Unhandled property: {property}. This is mostly likely due to the "
                "property not being defined in this PropertySet."
            )
        elif sep_point > 0:
            name = property[: sep_point - 1]
            index = property[sep_point:]
        else:
            name = property
            index = None

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
    act = PropertyMetadata(
        name="act",
        doc="Chemical Activity",
        units=units.dimensionless,
    )
    compress_fact = PropertyMetadata(
        name="compress_fact",
        doc="Compressiblity Factor",
        units=units.dimensionless,
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
    diffus_comp = PropertyMetadata(
        name="diffus_comp",
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
        units=units.dimensionless,
    )
    heat_capacity_ratio = PropertyMetadata(
        name="heat_capacity_ratio",
        doc="Heat Capacity Ration",
        units=units.dimensionless,
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
    isentropic_speed_sound_phase = PropertyMetadata(
        name="isentropic_speed_sound_phase",
        doc="Isentropic Speed of Sound",
        units="VELOCITY",
    )
    isothermal_speed_sound_phase = PropertyMetadata(
        name="isothermal_speed_sound_phase",
        doc="Isothermal Speed of Sound",
        units="VELOCITY",
    )
    henry = PropertyMetadata(
        name="henry",
        doc="Henry Constant",
        units=units.dimensionless,
        # TODO: Units are an issue here, as there are multiple ways to define this
    )
    mass_frac = PropertyMetadata(
        name="mass_frac",
        doc="Mass Fraction",
        units=units.dimensionless,
    )
    mole_frac = PropertyMetadata(
        name="mole_frac",
        doc="Mole Fraction",
        units=units.dimensionless,
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
        units=units.dimensionless,
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
        units=units.dimensionless,
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
        units=units.dimensionless,
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
        units=units.dimensionless,
    )
    log_conc_mol = PropertyMetadata(
        name="log_conc_mol",
        doc="Log of Molar Concentration",
        units=units.dimensionless,
    )
    log_mass_frac = PropertyMetadata(
        name="log_mass_frac_phase_comp",
        doc="Log of Mass Fractions",
        units=units.dimensionless,
    )
    log_molality = PropertyMetadata(
        name="log_molality",
        doc="Log of Molality",
        units=units.dimensionless,
    )
    log_mole_frac = PropertyMetadata(
        name="log_mole_frac",
        doc="Log of Mole Fractions",
        units=units.dimensionless,
    )
    log_mole_frac_pbub = PropertyMetadata(
        name="log_mole_frac_pbub",
        doc="Log of Mole Fractions at Bubble Point Pressure",
        units=units.dimensionless,
    )
    log_mole_frac_pdew = PropertyMetadata(
        name="log_mole_frac_pdew",
        doc="Log of Mole Fractions at Dew Point Pressure",
        units=units.dimensionless,
    )
    log_mole_frac_tbub = PropertyMetadata(
        name="log_mole_frac_tbub",
        doc="Log of Mole Fractions at Bubble Point Temperature",
        units=units.dimensionless,
    )
    log_mole_frac_tdew = PropertyMetadata(
        name="log_mole_frac_tdew",
        doc="Log of Mole Fractions at Dew Point Temperature",
        units=units.dimensionless,
    )
    log_pressure = PropertyMetadata(
        name="log_pressure",
        doc="Log of Pressure",
        units=units.dimensionless,
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
        units=units.dimensionless,
    )
    log_k_eq = PropertyMetadata(
        name="log_k_eq",
        doc="Log of Equilibrium Coefficient",
        units=units.dimensionless,
    )
    k_rxn = PropertyMetadata(
        name="k_rxn",
        doc="Rate Constant",
        units=units.dimensionless,
    )
    reaction_rate = PropertyMetadata(
        name="reaction_rate",
        doc="Reaction Rate",
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
    # Log terms
    log_act_phase_solvents = PropertyMetadata(
        name="log_act_phase_solvents",
        doc="Log of Activity of Solvents",
        units=units.dimensionless,
    )
