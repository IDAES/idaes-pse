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
"""This module provides functions to register and retrieve component information
"""

_components = {}


class _ComponentStruct(object):
    """Component registry entry structure"""

    def __init__(
        self,
        viscosity=False,
        thermal_conductivity=False,
        surface_tension=False,
        eos_ref=None,
        viscosity_ref=None,
        thermal_conductivity_ref=None,
        surface_tension_ref=None,
    ):
        """Create a component registry object.

        Args:
            viscosity (bool): True if viscosity model available
            thermal_conductivity (bool): True if thermal conductivity model available
            surface_tension (bool): True if surface tension model available
            eos_ref (str|list|None): Reference for equation of state model
            viscosity_ref (str|list|None): Reference for viscosity model
            thermal_conductivity_ref (str|list|None): Reference for thermal conductivity model
            surface_tension_ref (str|list|None): Reference for surface tension model

        Returns:
            _ComponentStruct"""
        # these are true if external functions to calculate are available
        self.viscosity = viscosity
        self.thermal_conductivity = thermal_conductivity
        self.surface_tension = surface_tension
        if isinstance(eos_ref, str) or eos_ref is None:
            self.eos_ref = eos_ref
        else:
            self.eos_ref = "\n".join(eos_ref)
        if isinstance(viscosity_ref, str) or viscosity_ref is None:
            self.viscosity_ref = viscosity_ref
        else:
            self.viscosity_ref = "\n".join(viscosity_ref)
        if (
            isinstance(thermal_conductivity_ref, str)
            or thermal_conductivity_ref is None
        ):
            self.thermal_conductivity_ref = thermal_conductivity_ref
        else:
            self.thermal_conductivity_ref = "\n".join(thermal_conductivity_ref)
        if isinstance(surface_tension_ref, str) or surface_tension_ref is None:
            self.surface_tension_ref = surface_tension_ref
        else:
            self.surface_tension_ref = "\n".join(surface_tension_ref)


def register_helmholtz_component(
    comp_str,
    viscosity=False,
    thermal_conductivity=False,
    surface_tension=False,
    eos_ref=None,
    viscosity_ref=None,
    thermal_conductivity_ref=None,
    surface_tension_ref=None,
):
    """Add a component to the registry

    Args:
        viscosity (bool): True if viscosity model available
        thermal_conductivity (bool): True if thermal conductivity model available
        surface_tension (bool): True if surface tension model available
        eos_ref (str|list|None): Reference for equation of state model
        viscosity_ref (str|list|None): Reference for viscosity model
        thermal_conductivity_ref (str|list|None): Reference for thermal conductivity model
        surface_tension_ref (str|list|None): Reference for surface tension model

    Returns:
        None
    """
    comp_str = comp_str.lower()
    _components[comp_str] = _ComponentStruct(
        viscosity=viscosity,
        thermal_conductivity=thermal_conductivity,
        surface_tension=surface_tension,
        eos_ref=eos_ref,
        viscosity_ref=viscosity_ref,
        thermal_conductivity_ref=thermal_conductivity_ref,
        surface_tension_ref=surface_tension_ref,
    )


def remove_component(comp_str):
    """Remove a component from the registry. If the component isn't
    registered, raise KeyError.

    Args:
        comp_str (str): Component to remove

    Returns:
        None

    Raises:
        KeyError: if the component isn't registered raise a KeyError
    """
    comp_str = comp_str.lower()
    del _components[comp_str]


def clear_component_registry():
    """Remove all components from registry"""
    _components.clear()


def registered_components():
    """Return a list of registered components"""
    return list(_components.keys())


def viscosity_available(comp_str):
    """Return whether a viscosity model is available for a component.

    Args:
        comp_str (str): component string

    Returns:
        bool: True if viscosity model is available, False if not
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].viscosity
    return False


def thermal_conductivity_available(comp_str):
    """Return whether a thermal conductivity model is available for a component.

    Args:
        comp_str (str): component string

    Returns:
        bool: True if thermal conductivity model is available, False if not
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].thermal_conductivity
    return False


def surface_tension_available(comp_str):
    """Return whether a surface tension model is available for a component.

    Args:
        comp_str (str): component string

    Returns:
        bool: True if surface tension model is available, False if not
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].surface_tension
    return False


def component_registered(comp_str):
    """Return whether a component is registered.

    Args:
        comp_str (str): component string

    Returns:
        bool: True if component is available, False if not
    """
    comp_str = comp_str.lower()
    return comp_str in _components


def eos_reference(comp_str):
    """Return the equation of state reference or None if
    not available

    Args:
        comp_str (str): component string

    Returns:
        (str|None): Equation of state reference
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].eos_ref
    return None


def viscosity_reference(comp_str):
    """Return the viscosity reference or None if
    not available

    Args:
        comp_str (str): component string

    Returns:
        (str|None): Viscosity reference
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].viscosity_ref
    return None


def thermal_conductivity_reference(comp_str):
    """Return the thermal conductivity reference or None if
    not available

    Args:
        comp_str (str): component string

    Returns:
        (str|None): Thermal conductivity reference
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].thermal_conductivity_ref
    return None


def surface_tension_reference(comp_str):
    """Return the surface tension reference or None if
    not available

    Args:
        comp_str (str): component string

    Returns:
        (str|None): Surface tension reference
    """
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].surface_tension_ref
    return None
