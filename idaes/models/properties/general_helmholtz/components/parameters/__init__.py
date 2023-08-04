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
"""This module gets the path of parameters files
"""

__author__ = "John Eslick"

import os
import re
import json

from pyomo.common.fileutils import this_file_dir

import idaes
from idaes.models.properties.general_helmholtz.components import (
    clear_component_registry,
    register_helmholtz_component,
)


def get_parameter_path():
    """Get the parameter file path

    Args:
        None

    Returns:
        str: path for parameter files
    """
    pth = idaes.cfg.properties.helmholtz.parameter_file_path
    if pth is None:
        pth = this_file_dir()
    return pth


def set_parameter_path(path):
    """Set the parameter file path, and register components found there.

    Args:
        path (str): Parameter file path

    Returns:
        None
    """
    idaes.properties.helmholtz.parameter_file_path = path


def auto_register():
    """Search the parameter file path for available components and register them"""
    pth = get_parameter_path()
    clear_component_registry()
    lst = os.listdir(pth)
    eos_set = set()  # list of components with EOS
    st_set = set()  # list of components with surface tension
    tcx_set = set()  # list of components with thermal conductivity
    visc_set = set()  # list of components with viscosity
    for fname in lst:
        rematch = re.fullmatch(r"(.*)_expressions_(.*).nl", fname)
        if rematch:
            comp = rematch.groups()[0]
            expr = rematch.groups()[1]
            if expr == "eos":
                eos_set.add(rematch.groups()[0])
            elif expr == "st":
                st_set.add(rematch.groups()[0])
            elif expr == "tcx":
                tcx_set.add(rematch.groups()[0])
            elif expr == "visc":
                visc_set.add(rematch.groups()[0])
        else:
            rematch = re.fullmatch(r"(.*).json", fname)
            eos_ref = None
            thermal_conductivity_ref = None
            viscosity_ref = None
            surface_tension_ref = None
            if rematch:
                with open(os.path.join(pth, fname), "r", encoding="utf-8") as fptr:
                    dct = json.load(fptr)
                try:
                    eos_ref = dct["eos"]["reference"]
                    if isinstance(eos_ref, list):
                        eos_ref = "\n".join(eos_ref)
                except KeyError:
                    eos_ref = None
                try:
                    thermal_conductivity_ref = dct["transport"]["thermal_conductivity"][
                        "reference"
                    ]
                    if isinstance(thermal_conductivity_ref, list):
                        thermal_conductivity_ref = "\n".join(eos_ref)
                except KeyError:
                    thermal_conductivity_ref = None
                try:
                    thermal_conductivity_ref = dct["transport"]["viscosity"][
                        "reference"
                    ]
                    if isinstance(viscosity_ref, list):
                        viscosity_ref = "\n".join(eos_ref)
                except KeyError:
                    viscosity_ref = None
                try:
                    surface_tension_ref = dct["transport"]["surface_tension"][
                        "reference"
                    ]
                    if isinstance(surface_tension_ref, list):
                        surface_tension_ref = "\n".join(eos_ref)
                except KeyError:
                    surface_tension_ref = None
    for comp in eos_set:
        register_helmholtz_component(
            comp,
            viscosity=comp in visc_set,
            thermal_conductivity=comp in tcx_set,
            surface_tension=comp in st_set,
            viscosity_ref=viscosity_ref,
            thermal_conductivity_ref=thermal_conductivity_ref,
            surface_tension_ref=surface_tension_ref,
        )
