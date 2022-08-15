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
Shared (utility) code for UI tests.
"""


def dict_diff(d1, d2, result=[], pfx=""):
    if isinstance(d1, list) and isinstance(d2, list):
        if len(d1) != len(d2):
            result.append(
                f"Array length at {pfx} first({len(d1)}) " f"!= second({len(d2)})"
            )
        else:
            pass  # good enough
    elif not isinstance(d1, dict) or not isinstance(d2, dict):
        if type(d1) == type(d2):
            same = None
            try:
                same = d1 == d2
            except:  # cannot compare them
                pass  # good enough
            if same is False:
                result.append(f"value at {pfx} first != second")
        else:
            result.append(
                f"type of value at {pfx} first({type(d1)}" f" != second({type(d2)}"
            )
    else:
        if set(d1.keys()) != set(d2.keys()):
            for k in d1:
                if k not in d2:
                    result.append(f"{pfx}.{k} in first, not in second")
            for k in d2:
                if k not in d1:
                    result.append(f"{pfx}.{k} in second, not in first")
        for k in d1:
            if k in d2:
                dict_diff(d1[k], d2[k], result=result, pfx=f"{pfx}.{k}")
    return result
