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


def almpickle(data, vargs):
    import pickle
    import sys

    # remove lambda function from results dictionary
    if len(vargs) == 1:
        fname = vargs[0]
    elif len(vargs) == 2:
        fname = vargs[0]
        dopick = vargs[1]
    else:
        fname = "alamo.pick"
        dopick = True

    if not isinstance(data, list):
        if "f(model)" in data.keys():
            data["f(model)"] = ()
        else:
            sys.stdout.write("Results dictionary does not contain lambda function")

    if dopick:
        pickle.dump(data, open(fname, "wb"))


def postpickle(data):
    # Recompile lambda function from model string and labels
    # import sympy
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import symbols, lambdify

    if isinstance(data["model"], type({})):
        data["f(model)"] = {}
        for olab in data["zlabels"]:
            model_str = data["model"].split("=")[1].replace("^", "**")
            data["f(model)"][olab] = lambdify(
                [symbols(data["xlabels"])], parse_expr(model_str), "numpy"
            )
    else:
        model_str = data["model"].split("=")[1].replace("^", "**")
        data["f(model)"] = lambdify(
            [symbols(data["xlabels"])], parse_expr(model_str), "numpy"
        )

    return data


def almunpickle(fname):
    # unpickle and relambdify
    from idaes.apps import alamopy_depr as alamopy

    # from alamopy import postpickle
    import pickle

    res = pickle.load(open(fname, "rb"))
    if "f(model)" in res[0].keys():
        res = alamopy.postpickle(res)
    return res
