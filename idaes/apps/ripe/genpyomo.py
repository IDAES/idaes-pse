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


def ripeomo(problem_data, ptype, fixarray, cc_int, task, pc, sharedata):
    # This subroutine creates a pyomo model to solve the RIPE MINLP
    # Inputs:
    # problem_data  - process data from ripemodel()
    # ptype         - (
    # fixarray      - array denoting valid stoichiometry/mechanism pairings
    # cc_int        - maximum cardinality
    # task          - 0 for presolve, 1 for RIPE call
    # pc            - process condition array from ripemodel()
    # sharedata     - shared data dctionary
    # Outputs:
    # riperes       - dioctionary containing parameter estimates

    aterm, target, sigma, bounds = problem_data
    # import re
    import pyomo.environ as pyo

    # from pyomo.opt import SolverFactory
    # import ripe
    import numpy as np

    gasc = sharedata["gasconst"]
    n, s, r, h = np.shape(aterm)
    # Import max bounds for parameter
    elb, eub, klb, kub = [
        bounds["e"]["min"],
        bounds["e"]["max"],
        bounds["k"]["min"],
        bounds["k"]["max"],
    ]
    if not isinstance(elb, type([])):
        elb = [elb] * r
    if not isinstance(eub, type([])):
        eub = [eub] * r
    if not isinstance(klb, type([])):
        klb = [klb] * r
    if not isinstance(kub, type([])):
        kub = [kub] * r

    # Initialize models
    opt = pyo.SolverFactory("baron")
    model = pyo.ConcreteModel()

    # Set i over observations
    model.i = pyo.RangeSet(n)
    # Set r over reaction mechanisms
    model.r = pyo.RangeSet(r)
    # Set s over species
    model.s = pyo.RangeSet(s)
    # Set h over stoichiometries
    model.h = pyo.RangeSet(h)

    # Define pre-exponential factor and binary variables
    model.k = pyo.Var(model.r, model.h, domain=pyo.NonNegativeReals)
    model.y = pyo.Var(model.r, model.h, domain=pyo.Binary)

    # Add set parameter and subset for elimination of dependent reactions
    # famcons involve reaction stoichiometries of the same order
    if sharedata["addfamcons"]:
        s_small, s_big, ncon, s_fam, nfam = sharedata["sdep"]
        model.fr = pyo.RangeSet(len(nfam))
        model.fcon = pyo.Param(
            model.fr, initialize=dict(((fr), nfam[fr - 1]) for fr in model.fr)
        )
        model.depfam = pyo.Set(
            model.fr, model.h, initialize=dict(((r), s_fam[r - 1]) for r in model.fr)
        )

        def famcon(model, fr):
            return (
                sum(sum(model.y[r, j] for r in model.r) for j in model.depfam[fr])
                <= model.fcon[fr]
            )

        model.famcons = pyo.Constraint(model.fr, rule=famcon)

    # depcons involve reaction stoiciometries of different orders
    if sharedata["adddepcons"]:
        s_small, s_big, ncon, s_fam, nfam = sharedata["sdep"]
        model.dr = pyo.RangeSet(len(ncon))
        model.scon = pyo.Param(
            model.dr, initialize=dict(((dr), ncon[dr - 1]) for dr in model.dr)
        )
        model.depsmall = pyo.Set(
            model.dr,
            within=model.h,
            initialize=dict(((r), s_small[r - 1]) for r in model.dr),
        )
        model.deplarge = pyo.Set(
            model.dr,
            within=model.h,
            initialize=dict(((r), s_big[r - 1]) for r in model.dr),
        )

        def depcon(model, dr):
            return (
                sum(sum(model.y[r, j] for r in model.r) for j in model.depsmall[dr])
                + model.scon[dr]
                * sum(sum(model.y[r, j] for r in model.r) for j in model.deplarge[dr])
                <= model.scon[dr]
            )

        model.depcons = pyo.Constraint(model.dr, rule=depcon)

    # ptype == 'arr' if temperature is included in kwargs and arrhenious relationships will be fit
    if ptype == "arr":
        model.E = pyo.Var(model.r, model.h, domain=pyo.Reals)
        model.T = pyo.Param(
            model.i, initialize=dict(((i), pc["T"][i - 1][0]) for i in model.i)
        )
        model.eub = pyo.Param(
            model.r, initialize=dict(((r), eub[r - 1]) for r in model.r)
        )
        model.elb = pyo.Param(
            model.r, initialize=dict(((r), elb[r - 1]) for r in model.r)
        )
        if "Tref" in pc.keys():
            model.Tref = pyo.Param(rule=pc["Tref"])

    model.aterm = pyo.Param(
        model.i,
        model.s,
        model.r,
        model.h,
        initialize=dict(
            ((i, s, r, h), aterm[i - 1, s - 1, r - 1, h - 1])
            for i in model.i
            for s in model.s
            for r in model.r
            for h in model.h
        ),
    )
    model.target = pyo.Param(
        model.i,
        model.s,
        initialize=dict(
            ((i, s), target[i - 1, s - 1]) for i in model.i for s in model.s
        ),
    )
    model.sigma = pyo.Param(
        model.i,
        model.s,
        initialize=dict(
            ((i, s), sigma[i - 1][s - 1]) for s in model.s for i in model.i
        ),
    )
    model.kub = pyo.Param(model.r, initialize=dict(((r), kub[r - 1]) for r in model.r))
    model.klb = pyo.Param(model.r, initialize=dict(((r), klb[r - 1]) for r in model.r))

    # Define isothermal objective function (no Arrhenious)
    def obj_expression(model):
        return sum(
            sum(
                (
                    model.target[i, s]
                    - sum(
                        sum(model.k[r, h] * model.aterm[i, s, r, h] for r in model.r)
                        for h in model.h
                    )
                )
                ** 2
                for i in model.i
            )
            / model.sigma[s]
            for s in model.s
        )

    # Define arrhenious objective function
    def arr_expression(model):
        return sum(
            sum(
                (
                    model.target[i, s]
                    - sum(
                        sum(
                            model.k[r, h]
                            * pyo.exp(-1.0 * model.E[r, h] / (gasc * model.T[i]))
                            * model.aterm[i, s, r, h]
                            for r in model.r
                        )
                        for h in model.h
                    )
                )
                ** 2
                for i in model.i
            )
            / model.sigma[s]
            for s in model.s
        )

    def refarr_expression(model):
        return sum(
            sum(
                (
                    model.target[i, s]
                    - sum(
                        sum(
                            model.k[r, h]
                            * pyo.exp(
                                -1.0
                                * (model.E[r, h] / gasc)
                                * (1 / (model.T[i]) - 1 / model.Tref)
                            )
                            * model.aterm[i, s, r, h]
                            for r in model.r
                        )
                        for h in model.h
                    )
                )
                ** 2
                for i in model.i
            )
            / model.sigma[s]
            for s in model.s
        )

    # Define additional objectives for wls, these could be combined with a simplified approach but are left separate for clarity for now
    def wls_obj_expression(model):
        return sum(
            sum(
                (
                    model.target[i, s]
                    - sum(
                        sum(model.k[r, h] * model.aterm[i, s, r, h] for r in model.r)
                        for h in model.h
                    )
                )
                ** 2
                / model.sigma[i, s]
                for i in model.i
            )
            for s in model.s
        )

    # Define arrhenious objective function
    def wls_arr_expression(model):
        return sum(
            sum(
                (
                    (
                        model.target[i, s]
                        - sum(
                            sum(
                                model.k[r, h]
                                * pyo.exp(-model.E[r, h] / (gasc * model.T[i]))
                                * (model.aterm[i, s, r, h])
                                for r in model.r
                            )
                            for h in model.h
                        )
                    )
                    ** 2
                    / model.sigma[i, s]
                )
                for i in model.i
            )
            for s in model.s
        )

    def wls_refarr_expression(model):
        #        return sum(sum((model.target[i,s] - sum(sum(model.k[r,h] * pyo.exp( -1.0 * model.e[r,h] * ( 1.0 / (8.314*model.T[i])))*model.aterm[i,s,r,h] for r in model.r) for h in model.h))**(2) / model.sigma[i,s] for i in model.i) for s in model.s)
        return sum(
            sum(
                (
                    model.target[i, s]
                    - sum(
                        sum(
                            model.k[r, h]
                            * pyo.exp(
                                -1.0
                                * (model.E[r, h] / gasc)
                                * (1 / (model.T[i]) - 1 / model.Tref)
                            )
                            * model.aterm[i, s, r, h]
                            for r in model.r
                        )
                        for h in model.h
                    )
                )
                ** 2
                / model.sigma[i, s]
                for i in model.i
            )
            for s in model.s
        )

    # Check for temperature dependencies
    if ptype == "arr":
        if "Tref" in pc.keys():
            model.OBJ = pyo.Objective(rule=wls_refarr_expression)
        else:
            model.OBJ = pyo.Objective(rule=wls_arr_expression)
    else:
        model.OBJ = pyo.Objective(rule=wls_obj_expression)

    # Define big-M constraints for k and possibly E
    def bigMub(model, r, h):
        return model.k[r, h] <= model.kub[r] * model.y[r, h]

    def bigMlb(model, r, h):
        return model.k[r, h] >= model.klb[r] * model.y[r, h]

    def bigEub(model, r, h):
        return model.E[r, h] <= model.eub[r] * model.y[r, h]

    def bigElb(model, r, h):
        return model.E[r, h] >= model.elb[r] * model.y[r, h]

    model.ubcon = pyo.Constraint(model.r, model.h, rule=bigMub)
    model.lbcon = pyo.Constraint(model.r, model.h, rule=bigMlb)
    if ptype == "arr":
        model.ubecon = pyo.Constraint(model.r, model.h, rule=bigEub)
        model.lbecon = pyo.Constraint(model.r, model.h, rule=bigElb)

    # Constriant defining one mechanisms per reaction stoichiometry
    def onemech(model, h):
        return sum(model.y[r, h] for r in model.r) <= 1.0

    if sharedata["onemechper"]:
        model.mechcon = pyo.Constraint(model.h, rule=onemech)

    # Cardinality constraint
    def ccon(model):
        return sum(sum(model.y[r, h] for h in model.h) for r in model.r) == cc_int

    # relax equality constraint
    def ccon_relax(model):
        return sum(sum(model.y[r, h] for h in model.h) for r in model.r) <= cc_int

    # hard cardinality vs soft for task 1 vs 0
    if task == 0:
        model.cardcon = pyo.Constraint(rule=ccon_relax)
    else:
        model.cardcon = pyo.Constraint(rule=ccon)

    # Fix array is used to fix binary variables of reactions/mechanisms that should not be considered
    nst, nm = np.shape(fixarray)
    for i in range(nst):
        for j in range(nm):
            if fixarray[i, j] == 0:
                model.y[j + 1, i + 1] = 0
                model.y[j + 1, i + 1].fixed = True

    # initialize results dictionary
    riperes = {}

    if task == 0:
        results = opt.solve(
            model,
            tee=sharedata["showpyomo"],
            options={"DeltaTerm": 1},
            keepfiles=sharedata["keepfiles"],
        )
    else:
        results = opt.solve(
            model,
            tee=sharedata["showpyomo"],
            keepfiles=sharedata["keepfiles"],
            options={
                "MaxTime": sharedata["maxmiptime"],
                "DeltaTerm": sharedata["deltaterm"],
            },
        )
    # Find residual values
    if ptype == "arr":
        if "Tref" in pc.keys():
            sr_list = [
                [
                    (
                        model.target[i, s]
                        - sum(
                            sum(
                                model.k[r, h].value
                                * pyo.exp(
                                    -1.0
                                    * (model.E[r, h].value / gasc)
                                    * (1 / (model.Tref) - 1 / (model.T[i]))
                                )
                                * model.aterm[i, s, r, h]
                                for r in model.r
                            )
                            for h in model.h
                        )
                    )
                    ** 2
                    for s in model.s
                ]
                for i in model.i
            ]
            ssr_list = ssr_list = [
                sum(
                    (
                        model.target[i, s]
                        - sum(
                            sum(
                                model.k[r, h].value
                                * pyo.exp(
                                    -1.0
                                    * (model.E[r, h].value / gasc)
                                    * (1 / (model.Tref) - 1 / (model.T[i]))
                                )
                                * model.aterm[i, s, r, h]
                                for r in model.r
                            )
                            for h in model.h
                        )
                    )
                    ** 2
                    for i in model.i
                )
                for s in model.s
            ]
        else:
            ssr_list = [
                sum(
                    (
                        model.target[i, s]
                        - sum(
                            sum(
                                model.k[r, h].value
                                * pyo.exp(
                                    -1.0 * model.E[r, h].value / (gasc * model.T[i])
                                )
                                * model.aterm[i, s, r, h]
                                for r in model.r
                            )
                            for h in model.h
                        )
                    )
                    ** 2
                    for i in model.i
                )
                for s in model.s
            ]
            sr_list = [
                [
                    (
                        model.target[i, s]
                        - sum(
                            sum(
                                model.k[r, h].value
                                * pyo.exp(
                                    -1.0 * model.E[r, h].value / (gasc * model.T[i])
                                )
                                * model.aterm[i, s, r, h]
                                for r in model.r
                            )
                            for h in model.h
                        )
                    )
                    ** 2
                    for s in model.s
                ]
                for i in model.i
            ]
    else:
        sr_list = [
            [
                (
                    model.target[i, s]
                    - sum(
                        sum(
                            model.k[r, h].value * model.aterm[i, s, r, h]
                            for r in model.r
                        )
                        for h in model.h
                    )
                )
                ** 2
                for s in model.s
            ]
            for i in model.i
        ]
        ssr_list = [
            sum(
                (
                    model.target[i, s]
                    - sum(
                        sum(
                            model.k[r, h].value * model.aterm[i, s, r, h]
                            for r in model.r
                        )
                        for h in model.h
                    )
                )
                ** 2
                for i in model.i
            )
            for s in model.s
        ]

    # Use results to estimate sigma values
    riperes["sigma"] = []
    for ssr in ssr_list:
        riperes["sigma"].append(ssr / np.max((1, float(n - cc_int))))
    riperes["wlsmat"] = []
    for sr in sr_list:
        riperes["wlsmat"].append([1 / np.max([0.000001, s]) for s in sr])
    model.solutions.store_to(results)

    # Populate results dictionary for post processing
    for key in ["ind", "k", "OBJ"]:
        riperes[key] = []
    if ptype == "arr":
        riperes["E"] = []

    # smallnum = 10 ** -6
    riperes["OBJ"] = results.Solution.Objective["OBJ"]["Value"]

    for v in results.Solution.Variable.keys():
        # find active binaries
        if "y" in v:
            # double test currently needed for obtained results
            if (
                results.Solution.Variable[v]["Value"] > 0.99
                and results.Solution.Variable[v]["Value"] < 1.1
            ):
                line = v.split(",")
                line = [line[0].split("[")[-1], line[1].split("]")[0]]
                line = [line[0], line[1]]
                riperes["ind"].append(line)
    # results may be in any order so we need to loop again and check the appropriate indicies
    if task == 0:
        riperes["maxk"] = 0.0
        riperes["maxe"] = 0.0

    for v in results.Solution.Variable.keys():
        line = v.split(",")
        line = [line[0].split("[")[-1], line[1].split("]")[0]]
        line = [line[0], line[1]]
        for inds in riperes["ind"]:
            if inds[0] == line[0] and inds[1] == line[1]:
                if "k" in v:
                    riperes["k"].append(
                        [
                            float(results.Solution.Variable[v]["Value"]),
                            [inds[0], inds[1]],
                        ]
                    )
                    if task == 0:
                        riperes["maxk"] = max(
                            riperes["maxk"],
                            float(results.Solution.Variable[v]["Value"]),
                        )
                if "E" in v:
                    riperes["E"].append(
                        [
                            float(results.Solution.Variable[v]["Value"]),
                            [inds[0], inds[1]],
                        ]
                    )
                    if task == 0:
                        riperes["maxe"] = max(
                            riperes["maxe"],
                            float(results.Solution.Variable[v]["Value"]),
                        )
    return riperes
