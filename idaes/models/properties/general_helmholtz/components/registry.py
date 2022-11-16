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
"""This module provides functions to register and retrieve component information
"""

import os
import idaes

_default_data_dir = os.path.join(idaes.bin_directory, "helm_data")
_default_data_dir = os.path.join(_default_data_dir, "")

_components = {}


class _ComponentStruct(object):
    def __init__(self, transport_module=None):
        self.transport_module = transport_module


def register_helmholtz_component(comp_str, transport_module=None):
    _components[comp_str] = _ComponentStruct(transport_module=transport_module)


def get_transport_module(comp_str):
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].transport_module
    return None


def component_registered(comp_str):
    comp_str = comp_str.lower()
    return comp_str in _components
