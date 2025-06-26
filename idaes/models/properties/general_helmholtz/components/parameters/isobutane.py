#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""Generate parameter and expression files for isobutane
"""

# Look at the documentation at
# https://github.com/IDAES/idaes-ext/blob/f60f28be990755749906567ac73fb1e0a188fa42/src/general_helmholtz/doc/documentation.pdf
# to see how to make new parameter files.
# The only thing you actually need is the .json file.
# and then you can find all the parameters online, e.g from the referenced papers
# or other papers that include helmholtz formulations for a component.
# the biggest difficulty is making sure the equations match the equations in helmholtz,
# otherwise you've got to define them yourself.

# Once you've made the .json file, you can use this script to generate the
# .nl files (as explained in the documentation.)

import os
from pyomo.common.fileutils import this_file_dir
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        WriteParameters
    """
    main_param_file = os.path.join(this_file_dir(), "isobutane.json")
    we = WriteParameters(parameters=main_param_file)
    we.write(dry_run=dry_run)
    return we


if __name__ == "__main__":
    main()
