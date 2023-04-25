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
"""Generate parameter and expression files for r227ea
"""

__author__ = "John Eslick"

from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def main():
    """Generate parameter and expression files."""
    we = WriteParameters(parameters="r227ea.json")
    we.write()


if __name__ == "__main__":
    main()
