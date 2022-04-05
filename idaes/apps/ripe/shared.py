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
Set the alamo and gams paths here.
"""
# import collections as col


# Input-specific options and mechanisms
# Mechanisms used in CLC models are enumerated as strings are parsed later
# dlist old clcmlist
dlist = [
    "c1",
    "c2",
    "c3",
    "c4",
    "c5",
    "c6",
    "c7",
    "c8",
    "c9",
    "c10",
    "c11",
    "c12",
    "c13",
    "c14",
    "c15",
    "c16",
    "c17",
    "c18",
    "c19",
    "c20",
    "c21",
]
rspace = {}
rspace["mechanisms"] = list(["iT", "ma", "rev", "homcat", "clc", "ads"] + dlist)
rspace["clc"] = list(dlist)
rspace["clcnames"] = list(
    [
        "Random Nucleation",
        "Power law order 2/3",
        "Power law order 2",
        "Power law order 4",
        "Power law order 3",
        "Power law order 1",
        "Avrami order 2",
        "Avrami order 3",
        "Avrami order 0.5",
        "Avrami order 1.5",
        "Avrami order 4",
        "Prout Tompkins",
        "Jander 3d diff",
        "Anti Jander 3d diff",
        "Valensi",
        "Parabolic",
        "Ginstling-Brountstein 3d diff",
        "zhuralev lesohkin tempelman",
        "Grain Model",
    ]
)
rspace["clcforms"] = list(
    [
        "1-X",
        "(n=2/3) nX^(n-1/n)",
        "(n=2) nX^(n-1/n)",
        "(n=4) nX^(n-1/n)",
        "(n=3) nX^(n-1/n)",
        "(n=1) nX^(n-1/n)",
        "(n=2) n(1-X)(-ln(1-X))^(n-1/n)",
        "(n=3) n(1-X)(-ln(1-X))^(n-1/n)",
        "(n=0.5) n(1-X)(-ln(1-X))^(n-1/n)",
        "(n=1.5) n(1-X)(-ln(1-X))^(n-1/n)",
        "(n=4) n(1-X)(-ln(1-X))^(n-1/n)",
        "X(1-X)",
        "3(1-X)^(1/3)(1/(1+X)^((-1/3)-1))",
        "3/2(1-X)^(2/3)(1/(1+X)^((-1/3)-1))",
        "1/(-ln(1-X))",
        "1/(2X)",
        "(3/2)(1-X)^(4/3)((1-X)^((-1/3)-1))^-1",
        "(3/2)((1-X)^((-1/3)-1))^-1",
        "(1-X)^(2/3)",
    ]
)
# rspace['pclabels']=list(['T','v','x0','flow'])

# Initialize additional dictionaries, and initialize options that will be passed between subroutines
sharedata = {}
dlist = list(["profiles", "otherdata", "pwalamo"])
for item in dlist:
    sharedata[item] = {}
dlist = list(["profiles", "tlo", "tup", "etrack", "err", "nmin"])
for item in dlist:
    sharedata["pwalamo"][item] = {}
# Specify options pertaining to applicaiton of pw-alamo
sharedata["pwalamo"]["err"] = 0.99
sharedata["pwalamo"]["nmin"] = 3
sharedata["otherdata"]["eobj"] = 2
sharedata["gasconst"] = 0.008314  # J/mol K
sharedata["bounds"] = {}
sharedata["bounds"]["k"] = {}
sharedata["bounds"]["e"] = {}
sharedata["bounds"]["k"]["min"] = 0.0
sharedata["bounds"]["k"]["max"] = 500.0
sharedata["bounds"]["e"]["min"] = 0.0
sharedata["bounds"]["e"]["max"] = 200.0
sharedata["df"] = 1
sharedata["npc"] = 1
sharedata["ratemthd"] = "deriv"
# possible arguments for kwargs in ripemodel
sharedata["kwargsin"] = [
    "mechanisms",
    "mech",
    "mechs",
    "stoich",
    "stoichiometry",
    "stoichs",
    "minlp_path",
    "showpyomo",
    "alamo_path",
    "nlp_path",
    "keepfiles",
    "return_model",
    "Tref",
    "sigma",
    "hide_output",
    "ccon",
    "wls",
    "Tr",
    "deltaterm",
    "ascale",
    "zscale",
    "onemechper",
    "expand_output",
    "temp",
    "time",
    "tr",
    "tref",
    "temperature",
    "Temperature",
    "Temp",
]
# Set paths here
sharedata["minlp_path"] = "baron"
sharedata["nlp_path"] = "ipopt"
sharedata["alamo_path"] = "alamo"
sharedata["showpyomo"] = False
sharedata["keepfiles"] = False
sharedata["return_model"] = False
sharedata["hide_output"] = False
sharedata["wls"] = False
sharedata["adddepcons"] = False
sharedata["addfamcons"] = False
sharedata["maxmiptime"] = 1000.0
sharedata["deltaterm"] = 0
sharedata["zscale"] = False
sharedata["ascale"] = False
sharedata["onemechper"] = True
sharedata["expand_output"] = False

# independent variables considered by atermconstruct
sharedata["ivars"] = ["t", "T", "x0", "flow", "vol", "other"]


#
debug = {}
debug["savescratch"] = False
debug["true"] = {}
debug["bignum"] = 10**20
debug["smallnum"] = 10**-10
debug["true"]["use"] = False
debug["true"]["fix"] = False
