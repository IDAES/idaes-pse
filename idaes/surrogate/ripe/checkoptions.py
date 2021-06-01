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
import sys


def checkoptions(kwargs, sharedata, debug):
    # This subroutine checks kwargs for ripe-specific parameters
    # the subroutine checks for valid arguments, but only a subset
    # are check to initialize sharedata
    # Inputs:
    # kwargs     - kwargs from ripemodel()
    # sharedata  - dictionary shared across subroutines
    # debug      - Additional shared dictionary
    # Outputs:
    # sharedata  - Originally defined in shared.py
    #    this dictionary contains information pertaining to RIPE's function

    inkeys = list(kwargs)  # kwargs.keys() python 2
    for key in inkeys:
        if key not in sharedata["kwargsin"] + sharedata["ivars"]:
            sys.exit(
                "Keyword argument "
                + key
                + " is not recognized. Consult the documnetation for appropriate arguments"
            )
        elif key == "minlp_path" or key == "alamo_path":
            sharedata[key] = kwargs[key]
        elif key == "keepfiles":
            sharedata[key] = kwargs[key]
        elif key == "showpyomo":
            sharedata[key] = kwargs[key]
        elif key == "Tref" or key == "Tr" or key == "tr" or key == "tref":
            sharedata["Tref"] = kwargs[key]
            sharedata["Tr"] = kwargs[key]
        elif key == "return_model":
            sharedata[key] = kwargs[key]
        elif key == "mechanisms" or key == "mech" or key == "mechs":
            kwargs["mech"] = kwargs[key]
        elif key == "stoichiometry" or key == "stoich" or key == "stoichs":
            kwargs["stoich"] = kwargs[key]
        elif key == "hide_output":
            sharedata[key] = kwargs[key]
        elif key == "deltaterm":
            sharedata[key] = kwargs[key]
        elif key == "zscale":
            sharedata[key] = kwargs[key]
        elif key == "ascale":
            sharedata[key] = kwargs[key]
        elif key == "onemechper":
            sharedata[key] = kwargs[key]
        elif key == "expand_output":
            sharedata[key] = kwargs[key]
        elif key == "time" or key == "t":
            kwargs["t"] = kwargs[key]
        elif (
            key == "Temp"
            or key == "Temperature"
            or key == "temp"
            or key == "temperature"
        ):
            kwargs["T"] = kwargs[key]

    return sharedata, kwargs
