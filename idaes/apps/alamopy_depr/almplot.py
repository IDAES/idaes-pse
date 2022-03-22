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
Plot doalamo output including confidence intervals if they are calculated.
"""

from idaes.apps.alamopy_depr import almerror


def almplot(res, show=True):
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        raise almerror.AlamoInputError(
            "Cannot plot, possibly missing matplotlib package"
        )

    model = res["model"]

    if type(model) is dict:
        for key in list(res["model"].keys()):
            model = res["model"][key].replace(" - ", " + -")

            model = model.split("=")[1]
            model = model.split(" + ")
            if (
                model[0] == " "
            ):  # if there are more than one terms, the first split is ' '
                model = model[1:]

            ndp = 100
            t = np.linspace(0.08, 1.7, ndp)
            out = np.zeros([3, ndp])
            clo = np.zeros(ndp)
            chi = np.zeros(ndp)
            coeff = np.zeros(ndp)

            for i in range(len(model)):
                thisterm = model[i].split(" * ")
                if thisterm[0].isdigit():
                    coeff[i] = float(model[i].split(" * ")[0])
                else:
                    coeff[i] = 1
                if "conf_inv" in res.keys():
                    if type(res["conf_inv"]) is dict:
                        # key = list(res['conf_inv'].keys())[0]
                        clo[i] = coeff[i] - float(
                            res["conf_inv"][key][i].split("+/-")[1]
                        )
                        chi[i] = coeff[i] + float(
                            res["conf_inv"][key][i].split("+/-")[1]
                        )
                    else:
                        clo[i] = coeff[i] - float(res["conf_inv"][i].split("+/-")[1])
                        chi[i] = coeff[i] + float(res["conf_inv"][i].split("+/-")[1])

            for i in range(ndp):
                out[0, i] = (
                    float(coeff[0]) * t[i] ** 1.2
                    - float(coeff[1]) * t[i] ** 2
                    - float(coeff[2])
                )
                if "conf_inv" in res.keys():  # If confidence intervals exist
                    out[1, i] = clo[0] * t[i] ** 1.2 - chi[1] * t[i] ** 2 - chi[2]
                    out[2, i] = chi[0] * t[i] ** 1.2 - clo[1] * t[i] ** 2 - clo[2]

            plt.plot(t, out[0], "b-")
            if "conf_inv" in res.keys():
                plt.plot(t, out[1], "r--", t, out[2], "r--")

            plt.xlabel(res["xlabels"][0])
            plt.ylabel(key)

            if show:
                plt.show()
    else:
        model = res["model"].replace(" - ", " + -")

        model = model.split("=")[1]
        model = model.split(" + ")
        if model[0] == " ":  # if there are more than one terms, the first split is ' '
            model = model[1:]
        ndp = 100
        t = np.linspace(0.08, 1.7, ndp)
        out = np.zeros([3, ndp])
        clo = np.zeros(ndp)
        chi = np.zeros(ndp)
        coeff = np.zeros(ndp)

        for i in range(len(model)):
            thisterm = model[i].split(" * ")
            if thisterm[0].isdigit():
                coeff[i] = float(model[i].split(" * ")[0])
            else:
                coeff[i] = 1
            if "conf_inv" in res.keys():
                if type(res["conf_inv"]) is dict:
                    firstKey = list(res["conf_inv"].keys())[0]
                    clo[i] = coeff[i] - float(
                        res["conf_inv"][firstKey][i].split("+/-")[1]
                    )
                    chi[i] = coeff[i] + float(
                        res["conf_inv"][firstKey][i].split("+/-")[1]
                    )
                else:
                    clo[i] = coeff[i] - float(res["conf_inv"][i].split("+/-")[1])
                    chi[i] = coeff[i] + float(res["conf_inv"][i].split("+/-")[1])

        for i in range(ndp):
            out[0, i] = (
                float(coeff[0]) * t[i] ** 1.2
                - float(coeff[1]) * t[i] ** 2
                - float(coeff[2])
            )
            if "conf_inv" in res.keys():  # If confidence intervals exist
                out[1, i] = clo[0] * t[i] ** 1.2 - chi[1] * t[i] ** 2 - chi[2]
                out[2, i] = chi[0] * t[i] ** 1.2 - clo[1] * t[i] ** 2 - clo[2]

        plt.plot(t, out[0], "b-")
        if "conf_inv" in res.keys():
            plt.plot(t, out[1], "r--", t, out[2], "r--")

        plt.xlabel(res["xlabels"][0])
        plt.ylabel(res["zlabels"][0])
        if show:
            plt.show()
