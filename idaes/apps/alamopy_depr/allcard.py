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

from importlib import reload

cvalsim = None  # could cause problems


def allcard(xdata, zdata, xval, zval, **kwargs):
    # enumerate all model cardinalities via ccmiqp and
    # use validation/cross-validaiton to determine
    from idaes.apps import alamopy_depr as alamopy

    # PYLINT-TODO-FIX alamopy.writethis.writethis doesn't seem to exist
    # pylint: disable=import-error
    from alamopy.writethis import writethis

    # pylint: enable=import-error
    import numpy as np
    import math

    # import sympy
    import random
    from random import shuffle

    # from scipy.optimize import curve_fit #2.7
    # from scipy.optimize import minimize #2.7
    from scipy.linalg import lstsq  # 2.7
    import os

    # from alamopy import almlsq as almfun
    # from alamopy import almlsqjac as almjac
    import time

    # from sympy.parsing.sympy_parser import parse_expr
    import sys

    random.seed(100)

    # kwargs=opdict
    trans = list(
        [
            "linfcns",
            "expfcns",
            "logfcns",
            "sinfcns",
            "cosfcns",
            "monomialpower",
            "multi2power",
            "multi3power",
            "ratiopower",
        ]
    )
    et = list(["multi2power", "multi3power", "ratiopower"])
    # datacc = {}
    ins = list(["cvfold"])
    for opt in ins:
        if opt in kwargs.keys():
            cvfold = kwargs["cvfold"]
    ntrans = 0
    ndata = np.shape(xdata)[0]
    ninputs = np.shape(xdata)[1]

    for t in trans:
        if t in kwargs.keys():
            if type(kwargs[t]) == list():
                nt = len(kwargs[t])
            else:
                nt = 1
            if t not in et:
                ntrans = ntrans + ninputs * nt
            else:
                ntrans = ntrans + math.factorial(ninputs) * nt
        else:
            kwargs[t] = 0
    # mseold = 0.0
    rmseold = 0.0
    startt = time.time()
    oldres = {}
    oldp = ()
    ntrans = min(ntrans, 1000)

    # split training and validationd ata before looping through cc
    tlist = list([])
    vlist = list([])
    # ndlist = range(ndata)
    if cvfold == "loo":
        for i in range(ndata):
            vlist.append(np.asarray([i]))
            temp = [x for x in range(ndata) if x != i]
            tlist.append(np.asarray(temp))
    else:
        temp = range(ndata)
        shuffle(temp)
        # tlist = np.asarray(temp)
        if cvfold == "valset":
            tlist = np.asarray(temp)[0 : int(0.7 * ndata)]
            vlist = np.asarray(temp)[int(0.7 * ndata) : -1]
        else:
            vlist = np.array_split(np.asarray(temp), int(cvfold))
            tlist = [1] * int(cvfold)
            for v in range(len(vlist)):
                tlist[v] = range(ndata)
                for this in vlist[v]:
                    tlist[v].remove(this)

    for ccon in list(range(1, ntrans + 1)):
        # reload(alamopy)
        # try:
        #     del cvalsim, almsim
        # except Exception:
        #     pass

        # res = alamopy.doalamo(xdata,zdata,kwargs.values())
        if cvfold != "valset":
            res = alamopy.doalamo(
                xdata,
                zdata,
                xval=xval,
                zval=zval,
                linfcns=kwargs["linfcns"],
                expfcns=kwargs["expfcns"],
                logfcns=kwargs["logfcns"],
                sinfcns=kwargs["sinfcns"],
                cosfcns=kwargs["cosfcns"],
                monomialpower=kwargs["monomialpower"],
                multi2power=kwargs["multi2power"],
                multi3power=kwargs["multi3power"],
                ratiopower=kwargs["ratiopower"],
                sigma=kwargs["sigma"],
                xlabels=kwargs["xlabels"],
                zlabels=kwargs["zlabels"],
                modeler=6,
                convpen=0,
                maxterms=ccon,
                almname=kwargs["almname"],
                expandoutput=kwargs["expandoutput"],
                xmax=kwargs["xmax"],
                xmin=kwargs["xmin"],
                savescratch=kwargs["savescratch"],
            )
        else:
            res = alamopy.doalamo(
                xdata[tlist, :],
                zdata[tlist],
                xval=xdata[vlist, :],
                zval=zdata[vlist],
                linfcns=kwargs["linfcns"],
                expfcns=kwargs["expfcns"],
                logfcns=kwargs["logfcns"],
                sinfcns=kwargs["sinfcns"],
                cosfcns=kwargs["cosfcns"],
                monomialpower=kwargs["monomialpower"],
                multi2power=kwargs["multi2power"],
                multi3power=kwargs["multi3power"],
                ratiopower=kwargs["ratiopower"],
                sigma=kwargs["sigma"],
                xlabels=kwargs["xlabels"],
                zlabels=kwargs["zlabels"],
                modeler=6,
                convpen=0,
                maxterms=ccon,
                almname=kwargs["almname"],
                expandoutput=kwargs["expandoutput"],
                xmax=kwargs["xmax"],
                xmin=kwargs["xmin"],
                savescratch=kwargs["savescratch"],
            )
        for fn in ["cvalsim.py", "cvalsim.pyc", "almsim.py", "almsim.pyc"]:
            try:
                os.remove(fn)
            except Exception:
                pass
        os.system("cp " + kwargs["almname"].split(".")[0] + "cv.py cvalsim.py")
        os.system("cp " + kwargs["almname"].split(".")[0] + "alm.py almsim.py")

        if ccon == 1:
            if ndata < 10:  # not enough data to do cross validation
                writethis("Not enough data to facilitate cross validation")
                endt = time.time()
                res["totaltime"] = endt - startt
                return res

        # import cvalsim
        import almsim

        xalm = alamopy.almfeatmat(xdata, ccon)
        # Lets do the cross validatione error
        mseval = 0.0
        rmse = {}
        if cvfold == "valset":
            resid = np.sum((zdata[vlist] - almsim.f(xdata[vlist, :])) ** 2)
            mseval = resid / float(len(vlist))
            params = "ALM params used for valset"
            rmse = {}
            rmse["val"] = res["rmseval"]
            rmse["train"] = res["rmse"]
        else:
            mseval = 0.0
            # track = 0
            for tl, vl in zip(tlist, vlist):
                # initb = [ 0.0 for x in range(ccon)]
                xd, xmax, xmin = alamopy.mapminmax(xalm)
                zd, zmax, zmin = alamopy.mapminmax(zdata)
                fitres = lstsq(xd[tl], zd[tl])
                params = fitres[0]
                resid = zd[vl] - np.matmul(xd[vl, :], params[:])
                resid = alamopy.remapminmax(resid, zmax, zmin)
                resid = sum(resid**2)
                if cvfold == "loo":
                    mseval = mseval + resid
                else:
                    mseval = mseval + resid / float(len(vl))
            mseval = mseval / float(len(vlist))
            rmse["val"] = np.sqrt(mseval)
            rmse["train"] = "not"
        if ccon > 1:
            if float(rmse["val"]) >= float(rmseold):
                sys.stdout.write(
                    "              Problem name   : " + kwargs["almname"] + "\n"
                )
                sys.stdout.write(
                    "  rMSEval : "
                    + str(rmse["val"])
                    + "    MSEold : "
                    + str(rmseold)
                    + "\n"
                )
                if cvfold == "valset":
                    sys.stdout.write(
                        "   Ntrain : "
                        + str(len(tlist))
                        + "    Nval : "
                        + str(len(vlist))
                        + "\n"
                    )
                else:
                    sys.stdout.write(
                        "   Ntrain : "
                        + str(len(tlist[0]))
                        + "    Nval : "
                        + str(len(vlist[0]))
                        + "\n"
                    )
                sys.stdout.write(
                    "       optimal model size is : " + str(ccon - 1) + "\n"
                )
                sys.stdout.write("    optimal coefficients are : " + str(oldp) + "\n")
                os.remove(kwargs["almname"].split(".")[0] + "alm.py")
                os.system(
                    "cp " + "oldalmsim.sv " + kwargs["almname"].split(".")[0] + "alm.py"
                )
                endt = time.time()
                oldres["totaltime"] = endt - startt
                return oldres
            elif ccon == ntrans:
                endt = time.time()
                sys.stdout.write("optimal model size is :" + str(ccon) + "\n")
                res["totaltime"] = endt - startt
                return res
            else:
                # mseold = mseval
                oldres = res
        else:
            # mseold = float(mseval)
            rmseold = float(rmse["val"])
            oldres = res
            oldp = params
        # keep track of alm model of old iteratoin
        try:
            os.remove("oldalmsim.sv")
        except Exception:
            pass
        os.system("mv almsim.py oldalmsim.sv")


def almlsq(params, *X):
    # return lsq objective
    # import numpy as np
    try:
        reload(cvalsim)
    except Exception:
        # import cvalsim
        from cvalsim import f as sim

    xdata, zdata = X[0]
    return sum((zdata - sim(xdata, params)) ** 2)


def almlsqjac(params, *X):
    try:
        reload(cvalsim)
    except Exception:
        # import cvalsim
        from cvalsim import f as sim
    import numpy as np

    xdata, zdata = X[0]
    rp = np.ones(len(params))
    for i in range(len(params)):
        dparams = [0.0 for x in range(len(params))]
        dparams[i] = 1.0
        rp[i] = sum(-2.0 * (zdata - sim(xdata, params)) * sim(xdata, dparams))
    return rp
    # return np.dot(params,np.matmul(np.transpose(xdata),xdata))+np.matmul(zdata,xdata)


def almfeatmat(X, nparams):
    import numpy as np

    try:
        reload(cvalsim)
    except Exception:
        # import cvalsim
        from cvalsim import f as sim

    dims = np.shape(X)
    fm = np.ones([dims[0], nparams])
    dparams = [0.0 for x in range(nparams)]
    for i in range(nparams):
        dparams[i] = 1.0
        fm[:, i] = sim(X, dparams)
        dparams[i] = 0.0
    return fm
