#!/usr/bin/env python

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
import math
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir
from pyomo.common.fileutils import find_library
from idaes.core.util.math import smooth_max
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
    # we.add(
    #     {
    #         # Add extra things here
    #     }
    # )
    we.write(dry_run=dry_run)
    return we


if __name__ == "__main__":
    main()
