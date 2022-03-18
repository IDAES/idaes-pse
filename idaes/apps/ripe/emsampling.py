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
# This file contains two subroutines
#  - construcutmodel(results, kwargs) : constructs a steady-state pyomo reactor model for use in error maximization sampling
#  - ems(riperes,sim,lb,ub,nspec,kwargs) : perform error maximization sampling
import numpy as np
import rbfopt
import os

# import sys
# import random
import datetime
import time

# pkg
from idaes.apps import ripe


def constructmodel(riperes, **kwargs):
    # Inputs:
    # riperes - output of ripemodel()
    # kwargs  - kwargs supplied to ripemodel()
    # Outputs:
    # ripesim - pyomo model

    # sharedata must be supplied explicitly
    if "sharedata" in kwargs.keys():
        sharedata = kwargs["sharedata"]
    else:
        sharedata = ripe.sharedata

    # Define gas constant for clarity later
    gc = sharedata["gasconst"]

    # determine type of kinetic form used
    if "E" in riperes.keys():
        acte = riperes["E"]
        if "Tref" in kwargs.keys():
            Tref = sharedata["Tref"]
            ptype = 3
        else:
            ptype = 2
    else:
        ptype = 1
    prek = riperes["k"]
    mechs = riperes["mechanisms"]
    stoichs = riperes["stoichiometry"]
    ns = len(stoichs[0])
    nr = len(mechs)
    # replace 'massact' with species specific mechanisms (required of rbfopt)
    for i in range(nr):
        if mechs[i] == "massact":
            mechs[i] = ripe.mechs.mechperstoich(mechs[i], stoichs[i])

    # Define a steady-state reactor model in pyomo using ripe resutls
    def ripesim(data):
        import pyomo.environ as pyo

        # ripesim expects input in a particular order
        model = pyo.ConcreteModel()
        # s is index over species
        model.s = pyo.RangeSet(ns)
        # r over reactions
        model.r = pyo.RangeSet(nr)
        # initialize variable from data
        model.conc0 = pyo.Param(
            model.s, initialize=dict(((s), data[s - 1]) for s in model.s)
        )
        if ptype > 1:
            model.T = pyo.Param(initialize=data[ns])
            model.E = pyo.Param(
                model.r, initialize=dict(((r), acte[r - 1]) for r in model.r)
            )

        model.nu = pyo.Param(
            model.r,
            model.s,
            initialize=dict(
                ((r, s), stoichs[r - 1][s - 1]) for r in model.r for s in model.s
            ),
        )
        model.conc = pyo.Var(model.s, domain=pyo.NonNegativeReals, initialize=1.0)
        try:
            model.flow = pyo.Param(
                model.s, initialize=dict(((s), data[ns + s]) for s in model.s)
            )
        except Exception:
            model.flow = pyo.Param(model.s, initialize=1.0)
        try:
            model.vol = pyo.Param(initialize=float(data[2 * ns + 1]))
        except Exception:
            model.vol = pyo.Param(initialize=1.0)
        model.k = pyo.Param(
            model.r, initialize=dict(((r), prek[r - 1]) for r in model.r)
        )
        model.rate = pyo.Var(model.r, initialize=0.0)
        model.dum = pyo.Var(model.s, initialize=0.0)

        # define predicted rates of generation
        # different problem types require different rate forms
        def prates_1(model, r):
            return model.rate[r] == model.k[r] * mechs[r - 1](*model.conc[:])

        def prates_2(model, r):
            indata = [model.conc[i] for i in model.s] + [data[ns]]
            return model.rate[r] == model.k[r] * pyo.exp(
                -(model.E[r] / (gc * model.T))
            ) * mechs[r - 1](*indata)

        def prates_3(model, r):
            indata = [model.conc[i] for i in model.s] + [data[ns]]
            return model.rate[r] == model.k[r] * pyo.exp(
                -(model.E[r] / (gc)) * (1 / (model.T) - 1 / (Tref))
            ) * mechs[r - 1](*indata)

        if ptype == 1:
            trule = prates_1
        elif ptype == 2:
            trule = prates_2
        else:
            trule = prates_3
        model.rc = pyo.Constraint(model.r, rule=trule)

        def balance(model, s):
            return 100.0 * model.dum[s] == model.flow[s] * (1.0 / model.vol) * (
                model.conc0[s] - model.conc[s]
            ) + sum([model.nu[i, s] * model.rate[i] for i in model.r])

        model.bal = pyo.Constraint(model.s, rule=balance)

        def obj(model):
            return sum([model.dum[i] ** 2 for i in model.s])

        model.OBJ = pyo.Objective(rule=obj)
        opt = pyo.SolverFactory(sharedata["minlp_path"])
        results = opt.solve(model, tee=sharedata["showpyomo"])
        model.solutions.store_to(results)
        conres = []
        for i in range(ns):
            conres.append(
                results.Solution.Variable["conc[" + str(i + 1) + "]"]["Value"]
            )
        return conres

    return ripesim


def ems(riperes, sim, lb, ub, nspec, **kwargs):
    # This subroutine performs error maximization sampling
    # Inputs:
    # riperes - results from ripemodel()
    # sim     - Original black-box simulator (callable in python)
    # lb/ub   - bounds for each independent variable
    # nspec   - number of species (can be different than #lb/ub)
    # Outputs:
    # x       - next best input point
    # errs    - absolute errors obtained on predicted point

    # Non-default options can be provided explicitly through kwargs
    if "sharedata" in kwargs.keys():
        sharedata = kwargs["sharedata"]
    else:
        sharedata = ripe.sharedata
    # poskeys = sharedata["ivars"] + ["x"]
    inkeys = list(set(kwargs.keys()) - set(["sharedata"]))
    ndim = len(lb)
    ns = nspec

    # Call atermconstruct.formatinputs to get input data
    # in a convenient form
    if "x" in inkeys:
        conc = kwargs["x"]
        dflag = True
    else:
        conc = [1] * ns
        dflag = False

    if "frac" in kwargs.keys():
        dofrac = True
        # conc[i] = kwargs["frac"](conc[i])  # i not defined?
    else:
        dofrac = False
    data, kwargs, fdata, pc, alldata = ripe.atermconstruct.formatinputs(conc, kwargs)

    # Handle inkeys exceptions and (possibly) generate a pyomo model
    if "Tref" in inkeys:
        Tref = sharedata["Tref"]
    if "res_sim" in inkeys:
        res_sim = kwargs["res_sim"]
    else:
        if "Tref" in inkeys:
            res_sim = constructmodel(riperes, sharedata=sharedata, Tref=Tref)
        else:
            res_sim = constructmodel(riperes, sharedata=sharedata)

    ndim = len(lb)
    nd, ns = np.shape(fdata)

    if "T" in inkeys:
        params = [[], []]
        params[0] = riperes["k"]
        params[1] = riperes["E"]
    else:
        params = riperes["k"]

    # Check for multiple requested points
    # if "nreq" in inkeys:
    #     nreq = kwargs["nreq"]
    # else:
    #     nreq = 1

    # subroutine for using fractional arguments (mole fracs or partial pressure)
    def apply_frac(x, doit=True):
        if doit:
            return kwargs["frac"](x[:])
        else:
            return x

    def check_fun(x):
        # Check fun is used to evaluate errors of the current model on provided data
        x[:] = apply_frac(x[:], dofrac)
        return -1.0 * np.sum(
            np.divide(
                np.power(np.subtract(x[:ns], res_sim(x[ns:])), 2), riperes["sigma"]
            )
        )

    def max_err(x):
        # max_err determines the absolute errors (max called later)
        x[:] = apply_frac(x[:], dofrac)
        errs = np.absolute(np.subtract(sim(x), res_sim(x)))[0]
        return errs

    def error_fun(x):
        # error fun is the black box objective
        x[:] = apply_frac(x[:], dofrac)
        sim_conc = sim(x[:])
        res_conc = res_sim(x[:])
        return -1.0 * np.sum(
            np.power(np.divide(np.subtract(sim_conc, res_conc), sim_conc), 2)
        )

    t_targets = np.zeros([nd, 1])
    if dflag:
        for i in range(nd):
            t_targets[i] = check_fun(alldata[i][:])

    bb = rbfopt.RbfoptUserBlackBox(
        ndim, np.array(lb), np.array(ub), np.array(["R"] * ndim), error_fun
    )
    # Request that any point returned be no closer than the other closest point
    distance = 10.0
    for i in range(nd):
        for j in range(nd):
            if i != j:
                distance = max(
                    distance,
                    np.linalg.norm(
                        alldata[i, ns : ns + ndim] - alldata[j, ns : ns + ndim]
                    ),
                )
    settings = []
    t = datetime.datetime.now()
    d_mult = 1.0
    settings = rbfopt.RbfoptSettings(
        nlp_solver_path=sharedata["nlp_path"],
        minlp_solver_path=sharedata["minlp_path"],
        print_solver_output=False,
        min_dist=d_mult * distance,
        algorithm="Gutmann",
        rand_seed=int(time.mktime(t.timetuple())),
    )  # random seed required for RBFopt functionality
    # Alternative settings used in testing saved for posterity
    #    settings = rbfopt.RbfoptSettings( do_infstep=True,nlp_solver_path='/usr/local/ipopt/3.10.2/ipopt/ipopt.pc', minlp_solver_path='~/baron/baron',print_solver_output=False, algorithm='MSRM',global_search_method = 'solver',rand_seed = int(time.mktime(t.timetuple())), min_dist = d_mult*distance)
    # ),#init_strategy='all_corners',
    # num_global_searches = 0,)
    # provide initialization data if flag is tripped
    if not dflag:
        #        print 'check conc : ',alldata[:,ns:ns+ndim],t_targets
        alg = rbfopt.RbfoptAlgorithm(
            settings,
            bb,
            do_init_strategy=False,
            init_node_pos=alldata[:, ns : ns + ndim],
            init_node_val=t_targets,
        )  # ,num_nodes_at_restart=nd)
    else:
        alg = rbfopt.RbfoptAlgorithm(settings, bb)

    # Call to rbfopt in a manner that terminal is not polluted
    f = open(os.devnull, "w")
    alg.set_output_stream(f)
    val, new_x, itercount, evalcount, fast_evalcount = alg.optimize(pause_after_iters=1)
    f.close()
    if len(new_x) < ndim:
        x = list(new_x) + [0] * (ndim - len(new_x))
    else:
        x = list(new_x)

    # Calculate relevant metrics and identify output responsible for error violation
    errs = max_err(x)
    # maxv = np.max(errs)
    # loc = np.argmax(errs)

    if "T" in inkeys:
        # prop_x = x[:-1]
        prop_t = x[-1]
    # else:
    #     prop_x = x
    #    print 'Maximum error on new proposed point : ', maxv,' on species # ',loc+1
    #    print 'Proposed initial concentrations', prop_x
    if "T" in inkeys:
        print("Proposed Temperature", prop_t)

    return [x, errs]
