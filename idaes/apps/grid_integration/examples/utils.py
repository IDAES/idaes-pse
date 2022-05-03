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

try:
    # Pyton 3.8+
    from importlib import resources
except ImportError:
    # Python 3.7
    import importlib_resources as resources
from pathlib import Path
import os

import pandas as pd


with resources.path("idaes.tests.prescient.5bus", "__init__.py") as pkg_file:
    prescient_5bus = Path(pkg_file).parent

# rts_gmlc_generator_dataframe = pd.read_csv("gen.csv")
rts_gmlc_generator_dataframe = pd.read_csv(os.path.join(prescient_5bus, "gen.csv"))
rts_gmlc_bus_dataframe = pd.read_csv(os.path.join(prescient_5bus, "bus.csv"))
