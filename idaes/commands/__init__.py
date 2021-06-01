###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
#
###############################################################################

import time

_command_import_start_time = time.time()

from idaes.commands.base import command_base as cb

import pkgutil
# import all the commands
for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):
    if module_name == "base":
        pass
    elif module_name.startswith("test"):
        pass
    else:
        loader.find_module(module_name).load_module(module_name)

_command_import_total_time = time.time() - _command_import_start_time
