##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import pkgutil

from idaes.commands.base import command_base as cb

# import all the commands
for loader, module_name, is_pkg in pkgutil.walk_packages(__path__):
    if module_name == "base":
        pass
    elif module_name.startswith("test"):
        pass
    else:
        loader.find_module(module_name).load_module(module_name)
