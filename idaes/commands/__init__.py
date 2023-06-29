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
# TODO: Missing doc strings
# pylint: disable=missing-module-docstring

import time

# Need to time imports, so must start timer before importing
# pylint: disable=wrong-import-position
_command_import_start_time = time.time()

import pkgutil
from idaes.commands.base import command_base as cb

# import all the commands
for loader, dotted_module_name, is_pkg in pkgutil.walk_packages(__path__):
    module_name_parts = dotted_module_name.split(".")
    module_name = module_name_parts[-1]
    is_test_module = module_name.startswith("test_")
    if dotted_module_name == "base":
        pass
    elif is_test_module:
        pass
    else:
        loader.find_module(module_name).load_module(module_name)

_command_import_total_time = time.time() - _command_import_start_time
