###############################################################################
# ** Copyright Notice **
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
###############################################################################
"""
This is a sample visualization script. This script builds a simple model and 
calls visualize("sample_visualization") in order to pop up a webpage with 
a sample visualization.
"""
import sys
import time
from idaes.generic_models.flowsheets.demo_flowsheet import build_flowsheet


def main():
    m = build_flowsheet()
    m.fs.visualize("sample_visualization")
    print("Starting server. Press Control-C to stop.")
    try:
        while 1:
            time.sleep(1)
    except KeyboardInterrupt:
        print("Stopped.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
