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


def almconfidence(data, *vargs):
    """
    This function calculates a covariance matrix
    # and confidence intervals of the estimated alamo regression coefficients
    """
    import numpy as np
    from scipy.stats import t  # 2.7
    from sympy.parsing.sympy_parser import parse_expr
    from sympy import symbols, lambdify

    # Initialize additional metric dictionaries
    data["covariance"] = {}
    data["conf_inv"] = {}

    if data.get("xdata", None) is None:
        xdata = vargs[0]
        # zdata = vargs[1]
    else:
        xdata = data["xdata"]
        # zdata = data['zdata']

    ndata = np.shape(xdata)[0]
    if isinstance(data["model"], type({})):
        for okey in data["model"].keys():
            model = data["model"][okey]
            model = model.split("=")[1]

            # split the model on +/- to isolate each linear term
            # This section is not currently in compliance with custom basis functions
            model = model.replace(" - ", " + ").split(" + ")
            if model[0] == " " or model[0] == "":
                model = model[1:]

            nlinterms = len(model)
            sensmat = np.zeros([ndata, nlinterms])
            covar = np.zeros([nlinterms, nlinterms])
            coeffs = np.zeros([nlinterms])

            for j in range(nlinterms):
                thisterm = model[j].split(" * ")
                if thisterm[0].isdigit():
                    coeffs[j] = float(eval(thisterm[0]))
                else:
                    coeffs[j] = 1
                thisterm = thisterm[-1]
                thislam = lambdify(
                    [symbols(data["xlabels"])],
                    parse_expr(thisterm.replace("^", "**")),
                    "numpy",
                )
                for i in range(ndata):
                    sensmat[i, j] = thislam(xdata[i])

            sigma = float(data["ssr"][okey]) / (float(ndata) - float(nlinterms))
            ci = np.zeros([nlinterms])
            covar = sigma * np.linalg.inv(np.matmul(np.transpose(sensmat), sensmat))
            for i in range(0, nlinterms):
                ci[i] = t.ppf(1 - 0.025, int(ndata) - nlinterms) * np.sqrt(covar[i][i])

            data["covariance"][okey] = covar
            data["conf_inv"][okey] = list()
            for j in range(nlinterms):
                data["conf_inv"][okey].append(
                    "B" + str(j + 1) + " : " + str(coeffs[j]) + "+/-" + str(ci[j])
                )
    else:
        model = data["model"]
        model.split("=")[1]
        model = model.split("=")[1]
        # split the model on +/- to isolate each linear term
        # This section is not currently in compliance with custom basis functions
        model = model.replace(" - ", " + ").split(" + ")
        while "" in model:
            model.remove("")
        while " " in model:
            model.remove(" ")
        nlinterms = len(model)
        sensmat = np.zeros([ndata, nlinterms])
        covar = np.zeros([nlinterms, nlinterms])
        coeffs = np.zeros([nlinterms])
        for j in range(nlinterms):
            if " * " in model[j]:
                thisterm = model[j].split(" * ")
                if thisterm[0].isdigit():
                    coeffs[j] = float(eval(thisterm[0]))
                else:
                    coeffs[j] = 1
                thisterm = thisterm[-1]
            else:
                thisterm = model[j]
                coeffs[j] = float(eval(thisterm))
            thislam = lambdify(
                [symbols(data["xlabels"])],
                parse_expr(thisterm.replace("^", "**")),
                "numpy",
            )
            for i in range(ndata):
                sensmat[i, j] = thislam(xdata[i])

        sigma = float(data["ssr"]) / (float(ndata) - float(nlinterms))

        covar = sigma * np.linalg.inv(np.matmul(np.transpose(sensmat), sensmat))
        ci = np.zeros([nlinterms])
        for i in range(0, nlinterms):
            ci[i] = t.ppf(1 - 0.025, int(ndata) - nlinterms) * np.sqrt(
                covar[i][i]
            )  # 2.7

        data["covariance"] = covar
        data["conf_inv"] = list()
        for j in range(nlinterms):
            data["conf_inv"].append(
                "B" + str(j + 1) + " : " + str(coeffs[j]) + "+/-" + str(ci[j])
            )
    return data
