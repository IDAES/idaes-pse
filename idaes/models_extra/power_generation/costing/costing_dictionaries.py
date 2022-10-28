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
"""
Python script to read costing components
This script reads the library of costing components (scaled cost, reference
parameters, costing exponents, etc.) from the json files.
First, open json file, then create a python dictionary that gets imported into
power_plant_capcost.py

Two python dictionaries that are loaded:
* BB_costing_data
* sCO2_costing_params

"""
__author__ = "Costing Team (A. Noring, B. Paul, D. Caballero, and M. Zamarripa)"
__version__ = "1.0.0"

import os
import json
from pyomo.common.fileutils import this_file_dir

directory = this_file_dir()


def load_BB_costing_dictionary():
    """
    The costing data dictionary contains information from the BBR4 COE
    spreadsheet and from the QGESS on capital cost scaling methodology
    (DOE/NETL-2019/1784). Specifically it includes scaling exponents,
    valid ranges for the scaled parameter, and units for those ranges.
    It is important to note the units only apply to the ranges and are
    not neccessarily the units that the reference parameter value will
    be given in.. It includes the total plant cost (TPC), reference
    parameter value, and units for that value.

    Some accounts are costed using two different reference parameters, these
    accounts have been divided into two separate accounts following the naming
    convention x.x.a and x.x.b.

    This dictionary is nested with the following structure:
    tech type --> CCS --> account --> property name --> property values
    """

    # assuming the dictionary exists, load it so it is importable when called
    with open(os.path.join(directory, "BB_costing_data.json"), "r") as file:
        BB_costing_params = json.load(file)
    return BB_costing_params


def load_sCO2_costing_dictionary():
    with open(os.path.join(directory, "sCO2_costing_parameters.json"), "r") as file:
        sCO2_costing_params = json.load(file)
    return sCO2_costing_params
