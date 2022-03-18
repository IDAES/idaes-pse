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
from idaes.apps import ripe
import numpy as np

# import itertools as it
import sys


def ripemodel(data, **kwargs):
    # import pyomo.environ

    #    from ripe import mechs, kinforms, setshared
    # import os

    sharedata, debug = ripe.sharedata, ripe.debug
    sharedata, kwargs = ripe.checkoptions(kwargs, sharedata, debug)

    # input to ripe specifies knietic mechanisms for reaction stoichiometries
    # Get mechs formats or creates these kinetic mechanisms
    stoich, mechanisms, fixarray, ncons, mechlist = ripe.mechs.getmechs(kwargs)

    # Find dependent stoichiometries

    if False:
        #    if np.linalg.matrix_rank(np.array(stoich)) < len(stoich):
        sharedata["adddepcons"] = True
        dep_small, dep_large, ncon, fam_dep, nfam = ripe.bounds.stoich_cons(stoich)
        if len(ncon) > 0:
            sharedata["adddepcons"] = True
        if len(nfam) > 0:
            sharedata["addfamcons"] = True
    else:
        dep_small, dep_large, ncon, fam_dep, nfam = [False, False, False, False, False]
        sharedata["adddepcons"] = False
        sharedata["addfamcons"] = False
    sharedata["sdep"] = [dep_small, dep_large, ncon, fam_dep, nfam]

    # Construction of the activity term is handled here
    # i.e. - rate = k exp(-E/RT) * Aterm
    aterm, fdata, pc, data, scales = ripe.atermconstruct.makeaterm(
        data, stoich, mechanisms, kwargs, ncons, mechlist, fixarray, sharedata
    )
    # remove stoichiometric-dependent flag from mechs
    mechanisms = mechanisms[0]

    # Aterm is a 4-dimensional array where
    # n = number of process conditions (t included for dynamic data)
    # ns = number of species
    # nm = number of mechanisms
    # nh = number of stoichiometries
    n, ns, nm, nh = np.shape(aterm)

    # Need to calculate target values
    targets, s_targets, savetargets, sharedata = ripe.targets.gentargets(
        data, kwargs, fdata, pc, sharedata
    )

    # ptype determines arrhenious vs linear kinetic models

    if "T" in kwargs.keys():
        ptype = "arr"
    elif ns > 1:
        ptype = "iT"
    else:
        ptype = "simple"

    # Initialize sigma then calculate variance of residuals
    if "sigma" in kwargs.keys():
        if isinstance(kwargs["sigma"], type(0.0)) or isinstance(
            kwargs["sigma"], type(0)
        ):
            kwargs["sigma"] = np.array(kwargs["sigma"])
        sigma = np.squeeze(kwargs["sigma"])
        dshape = np.shape(sigma)
        if len(dshape) == 1:
            # First determine if the sigma specified is a single value or vector
            if dshape[0] == 1:
                sigma = np.multiply(np.ones([n, ns]), sigma[0])
                sharedata["wls"] = False
            elif dshape[0] == ns:
                sigma = np.multiply(np.ones([n, ns]), sigma[:])
                sharedata["wls"] = False
            else:
                sys.exit("Invalid shape specified for input argument sigma")
        elif len(dshape) == 2:
            sigma = np.reshape(sigma, (n, ns))
            sharedata["wls"] = True

        # Replace specified zero sigma with small number
        for i in range(dshape[0]):
            for j in range(dshape[1]):
                if sigma[i, j] < 0.0:
                    sys.exit("Negative varaince defined by user")
                elif sigma[i, j] == 0.0:
                    sigma[i, j] = 10**-8
        # Initial call of ripeomo to get refined bounds
        sharedata = ripe.bounds.get_bounds(
            [
                [aterm, targets, sigma, sharedata["bounds"]],
                ptype,
                fixarray,
                ncons,
                0,
                pc,
                sharedata,
            ]
        )
        # Need to claculate refined variable bounds here
    else:
        sigma = [[1.0] * int(ns)] * int(n)
        ccon = np.min([ns, ncons, n, nm, nh])
        sys.stdout.write(
            "Calcuating residual variance with data provided with cc = "
            + str(ccon)
            + "\n"
        )
        sigresults = ripe.genpyomo.ripeomo(
            [aterm, targets, sigma, sharedata["bounds"]],
            ptype,
            fixarray,
            ccon,
            0,
            pc,
            sharedata,
        )
        sigma = [sigresults["sigma"] * int(ns)] * int(n)
        if sharedata["wls"]:
            sigma = sigresults["wlsmat"]
    sigma = np.array(sigma)
    dshape = np.shape(sigma)

    for i in range(dshape[0]):
        for j in range(dshape[1]):
            if sigma[i, j] < 10**-8:
                sigma[i, j] = 10**-8

    # w_norm = sum(sum(sigma))
    savesig = sigma
    if sharedata["zscale"]:
        for i in range(n):
            for j in range(ns):
                sigma[i, j] = sigma[i, j] / (s_targets**2)
            # if sigma[i,j] < debug['smallnum'] or savesig

    if "ccon" in kwargs.keys():
        ccon = int(kwargs["ccon"][0][0])
        sys.stdout.write(
            "Solving RIPE model with cardinality constraint = " + str(ccon) + "\n"
        )
        results = ripe.genpyomo.ripeomo(
            [aterm, targets, sigma, sharedata["bounds"]],
            ptype,
            fixarray,
            ccon,
            1,
            pc,
            sharedata,
        )
    else:
        # ccon is not specified, use BIC to size model
        if ptype == "arr":
            ccon_list = range(1, 1 + np.min([nh, ns, int(np.floor((n - 1) / 2))]))
        else:
            ccon_list = range(1, 1 + np.min([nh, ns, int(n - 1)]))
        bic = []
        # Overwriting results variable here
        d_results = dict.fromkeys([0] + list(ccon_list))
        sys.stdout.write(
            "   ---- Calculating null values for model selection ----    \n"
        )
        # Getting sigma estimates from ~unbiased estimates
        d_results[0] = ripe.genpyomo.ripeomo(
            [aterm, targets, sigma, sharedata["bounds"]],
            ptype,
            fixarray,
            0,
            1,
            pc,
            sharedata,
        )
        bic.append(d_results[0]["OBJ"])
        sys.stdout.write(" - Null model BIC = " + str(bic[-1]) + "\n")
        for ccon in ccon_list:
            sys.stdout.write(
                " - Solving RIPE model with cardinality constraint = "
                + str(ccon)
                + " - \n"
            )
            d_results[ccon] = ripe.genpyomo.ripeomo(
                [aterm, targets, sigma, sharedata["bounds"]],
                ptype,
                fixarray,
                ccon,
                1,
                pc,
                sharedata,
            )
            bic.append(float(d_results[ccon]["OBJ"]) + np.log(n * ns) * ccon)
            sys.stdout.write(
                " - " + str(ccon) + "-term model BIC = " + str(bic[-1]) + "\n"
            )
            if bic[-1] > bic[-2]:
                ccon = ccon - 1
                break

        results = d_results[ccon]

    if "sigma" not in kwargs.keys():
        sigma = np.array(results["wlsmat"])
    else:
        sigma = savesig
    kres, eres, stoichres, mechres, r2 = ripe.confinv(
        aterm, targets, results, ccon, ptype, pc, sharedata, sigma
    )  # ,scales,s_targets)
    #    kres, eres, stoichres, mechres,r2 = ripe.confinv(scales['aterm_u'],savetargets,results,ccon,ptype,pc,sharedata,sigma,scales,s_targets)
    for i in range(len(mechres)):
        in1 = int(stoichres[i]) - 1
        in2 = int(mechres[i]) - 1
        stoichres[i] = stoich[in1]
        mechres[i] = mechlist[in2]
        if sharedata["ascale"]:
            kres[i][0] = kres[i][0] / float(scales["nfactor"][in2, in1])
            kres[i][1] = kres[i][1] / float(scales["nfactor"][in2, in1])
        if sharedata["zscale"]:
            kres[i][0] = s_targets * kres[i][0]
            kres[i][1] = s_targets * kres[i][1]

    s_kres, s_eres, ci_k, ci_e = print_results(
        kres, eres, stoichres, mechres, r2, sharedata
    )

    main_return = {}
    # Additional information is included if expand_outputs is true
    if sharedata["expand_output"]:
        main_return["wls"] = results["wlsmat"]
        if "sigma" in kwargs.keys():
            main_return["sigma"] = sigma
        else:
            main_return["sigma"] = results["sigma"]

    main_return["k"] = s_kres
    main_return["conf_inv"] = ci_k
    main_return["stoichiometry"] = stoichres
    main_return["mechanisms"] = mechres
    if eres != []:
        main_return["E"] = s_eres
        main_return["conf_inv"].append([ci_e])
    return main_return


def ripewrite(sd, string):
    import sys

    with open("riperesults.txt", "a") as a:
        a.write(string + "\n")
    if not sd["hide_output"]:
        sys.stdout.write(string + "\n")


def print_results(kres, eres, stoichres, mechres, r2, sd):
    import sys

    cc = len(mechres)

    s_kres = []
    s_eres = []
    ci_k = []
    ci_e = []
    tstring = (
        "  ----  RIPE selected "
        + str(cc)
        + " models in the optimal reaction network  ----  \n"
    )
    sys.stdout.write(tstring)
    with open("riperesults.txt", "w") as a:
        a.write(tstring)
    ripewrite(sd, " -+-  Overall R2 of selected model : " + str(r2))
    for i in range(cc):
        ripewrite(sd, " - Stochiometry selected for reaction " + str(i + 1))
        ripewrite(sd, "         " + str(stoichres[i]))
        ripewrite(sd, "   Mechanism selected for reaction " + str(i + 1))
        ripewrite(sd, "         " + str(mechres[i]))
        ripewrite(sd, "   Kinetic rate constant & 95% confidence interval")
        ripewrite(sd, "         " + str(kres[i][0]) + " +/- " + str(kres[i][1]))
        ci_k.append(kres[i][1])
        s_kres.append(kres[i][0])
        if eres != []:
            ripewrite(sd, "   Activation energy & 95% confidence interval")
            ripewrite(sd, "        " + str(eres[i][0]) + " +/- " + str(eres[i][1]))
            s_eres.append(eres[i][1])
            ci_e.append(eres[i][1])
        ripewrite(sd, "")

    return [s_kres, s_eres, ci_k, ci_e]
