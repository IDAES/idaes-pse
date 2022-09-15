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
power_plant_costing.py

Three python dictionaries that are loaded:
* BB_costing_exponents
* BB_costing_params
* sCO2_costing_params

"""
__author__ = "Costing Team (A. Noring, M. Zamarripa, B. Paul, D. Caballero)"
__version__ = "1.0.0"

import os
import json
from pyomo.common.fileutils import this_file_dir

directory = this_file_dir()


"""
The costing exponents dictionary contains information from the QGESS on
capital cost scaling methodology (DOE/NETL-2019/1784). Specifically it includes
scaling exponents, valid ranges for the scaled parameter, and units for those
ranges. It is important to note the units only apply to the ranges and are not
neccessarily the units that the reference parameter value will be given in.
This dictionary is nested with the following structure:

tech type --> account --> property name --> property value
"""

with open(os.path.join(directory, "BB_costing_exponents.json"), "r") as file:
    BB_costing_exponents = json.load(file)


"""
The costing params dictionary contains information from the BBR4 COE
spreadsheet. It includes the total plant cost (TPC), reference parameter value,
and units for that value.

Some accounts are costed using two different reference parameters, these
accounts have been divided into two separate accounts following the naming
convention x.x.a and x.x.b.

This dictionary is nested with the following structure:
tech type --> CCS --> account --> property name --> property values
"""
with open(os.path.join(directory, "BB_costing_parameters.json"), "r") as file:
    BB_costing_params = json.load(file)


with open(os.path.join(directory, "sCO2_costing_parameters.json"), "r") as file:
    sCO2_costing_params = json.load(file)

if not os.path.exists(
    os.path.join(directory, "BB_costing_data.json")
):  # make the dictionary
    # remove this section later, this is just a way to "zip together" the two BB_costing files

    BB_costing_data = BB_costing_params

    for tech in BB_costing_data.keys():  # do one technology at a time
        for ccs in ["A", "B"]:
            if ccs in BB_costing_data[tech]:  # check if CCS = A, for indexing
                accounts_dict = BB_costing_data[tech][ccs]  # shorter alias
                for account in accounts_dict.keys():  # do one account at a time
                    accounts_dict[account][
                        "BEC_units"
                    ] = "K$2018"  # add BEC units as thousands of 2018 USD
                    for accountkey in BB_costing_exponents[tech][
                        account
                    ].keys():  # get one " exponents"account property at a time
                        accounts_dict[account][accountkey] = BB_costing_exponents[tech][
                            account
                        ][accountkey]
                    sorted_accountkeys = sorted(
                        accounts_dict[account]
                    )  # now, sort the accountkeys alphabetically within each account
                    accounts_dict[account] = {
                        key: accounts_dict[account][key] for key in sorted_accountkeys
                    }  # re-add the keys in alphabetical order
                BB_costing_data[tech][
                    ccs
                ] = accounts_dict  # use the alias to update the original dictionary

    with open(os.path.join(directory, "BB_costing_data.json"), "w") as outfile:
        json.dump(BB_costing_data, outfile)
    print("Success! New costing dictionary generated.")

# assuming the dictionary exists, load it so it is importable when called
with open(os.path.join(directory, "BB_costing_data.json"), "r") as file:
    BB_costing_params = json.load(file)
