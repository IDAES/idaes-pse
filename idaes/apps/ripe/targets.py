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
# This subroutine generates target rates of genration for each species
# Currently, alamo is called explicitly in RIPE
# alamopy dependencies exist through calls to alamopy.multos in order
# to manipulate files alamo creates
# soubroutines in this file
# doalamo         - generate alamo models for dynamic targets
# dopwalamo       - wrapper for pw alamo application
# gentartgs       - generate target values from process data
# sstargets       - called in gentargets for steady-state
# dynamictargets  - called in gentargets for dynamic data

import numpy as np
import os


def doalamo(data, time, xlo, xup, rpc, rspec, sharedata):
    # This subroutine calls alamo to generate targets for dynamic problems
    # this call could be replaced by a call to alamopy.doalamo()
    # Inputs:
    # data     - process data from ripemodel()
    # time     - time from the pc dictionary in ripemodel()
    # xlo/xup  - bounds for time in dynamic problems
    # rpc      - list containing range of process conditions
    # rspec    - list of species
    # Outputs:
    # profiles - dictionary of prfiles from recursive piece-wise alamo
    # r2ret    - r2 of submodels

    # import sympy
    from idaes.apps.alamopy_depr import multos as mos
    from idaes.apps.ripe import debug

    # from sympy import symbols

    if len(np.shape(data)) == 2:
        data = np.expand_dims(data, axis=-1)

    # Convert data to python index
    xlo = xlo - 1
    xup = xup - 1

    r2ret = {}
    profiles = {}
    # st = symbols("st")
    # Build alamo models for each process condition
    for p in rpc:
        r2ret[p] = {}
        profiles[p] = {}
        with open("temp.alm", "w") as a:
            a.write("linfcns 1\n")
            a.write("funform 1\n")
            a.write("trace 1\n")
            a.write("tracefname temptrace.trc\n")
            a.write("solvemip 1\n")
            a.write("ninputs 1\n")
            a.write("noutputs " + str(len(rspec)) + "\n")
            a.write("xlabels st\n")
            a.write("modeler 1\n")
            a.write("monomialpower 0.5 1 1.5 2 \n")
            a.write("logfcns 1\n")
            a.write("expfcns 1\n")
            a.write("constant 0\n")
            a.write("xmin 0\n")
            a.write("xmax " + str(time[-1]) + "\n")
            a.write("cvxbic 1\n")
            a.write("initialpoints " + str(len(range(xlo, xup))) + "\n")
            a.write("ndata " + str(len(range(xlo, xup))) + "\n")
            a.write("begin_data\n")
            for i in range(xlo, xup):  #
                tstr = str(time[i]) + " "
                for j in rspec:
                    tstr = tstr + str(data[i, j, p]) + " "
                a.write(tstr + "\n")
            a.write("end_data\n")
        os.system(sharedata["alamo_path"] + " temp.alm > tmpscratch")

    f = open("temptrace.trc")
    lf = f.read()
    f.close()
    lf2 = lf.split("\n")
    z = 0
    for i in rpc:
        for j in rspec:
            if "#filename" in lf2[z]:
                z = z + 1
            profiles[i][j] = lf2[z].split(",")[-1].split(" = ")[-1]
            r2ret[i][j] = lf2[z].split(",")[18]
            z = z + 1
    mos.copyfile("sv.alm", "temp.alm")
    if not debug["savescratch"]:
        mos.deletefile("sv.alm", "tmpscratch", "temp.alm", "temptrace.trc", "temp.lst")

    return [profiles, r2ret, sharedata]


def dopwalamo(data, tdata, profile, r2, xlo, xup, rpc, rspec, sharedata):
    from sympy import symbols
    from sympy.parsing.sympy_parser import parse_expr

    # debug = ripe.debug

    # python syntax of data indicies
    xlm = xlo - 1
    xum = xup - 1

    st = symbols("st")
    for pc in rpc:
        for s in rspec:
            amodel = parse_expr(profile[pc][s])
            if float(r2[pc][s]) > sharedata["pwalamo"]["err"]:
                sharedata["pwalamo"]["profiles"][pc][s].append(profile[pc][s])
                sharedata["pwalamo"]["tlo"][pc][s].append(xlo)
                sharedata["pwalamo"]["tup"][pc][s].append(xup)
                sharedata["pwalamo"]["etrack"][pc][s].append(r2)
            elif (xup - xlo) > 2 * sharedata["pwalamo"]["nmin"]:
                cutoff = (1.0 - float(sharedata["pwalamo"]["err"])) * np.sum(
                    (data[xlm:xum, s, pc] - np.average(data[xlm:xum, s, pc])) ** 2
                )
                temp = 0
                i = xlm
                while temp < 0.5 * cutoff:
                    temp = temp + (data[i, s, pc] - amodel.subs(st, tdata[i])) ** 2
                    i = i + 1
                if (i - xlo) < sharedata["pwalamo"]["nmin"]:
                    i = xlo + sharedata["pwalamo"]["nmin"]
                if (xup - i) > sharedata["pwalamo"]["nmin"]:
                    i = xup - sharedata["pwalamo"]["nmin"]
                [lp1, lr1, sharedata] = doalamo(
                    data, tdata, xlo, i, [pc], [s], sharedata
                )
                sharedata = dopwalamo(
                    data, tdata, lp1, lr1, xlo, i, [pc], [s], sharedata
                )
                [lp2, lr2, sharedata] = doalamo(
                    data, tdata, i, xup, [pc], [s], sharedata
                )
                sharedata = dopwalamo(
                    data, tdata, lp2, lr2, i, xup, [pc], [s], sharedata
                )
            else:
                sharedata["pwalamo"]["profiles"][pc][s].append(profile[pc][s])
                sharedata["pwalamo"]["tlo"][pc][s].append(xlo)
                sharedata["pwalamo"]["tup"][pc][s].append(xup)
                sharedata["pwalamo"]["etrack"][pc][s].append(r2[pc][s])
    return sharedata


def gentargets(data, kwargs, fdata, pc, sd):
    # this subroutine is a wrapper for gene
    # Determine target rate values
    # Sharedata is passed to track pw-alamo progress
    if "t" in kwargs.keys():
        targets, sd = dynamictargets(data, fdata, pc, sd)
    else:
        targets = sstargets(fdata, pc)

    if sd["zscale"]:
        s_targets = np.max([np.abs(np.ndarray.flatten(targets))])
        savetargets = targets
        targets = np.divide(targets, s_targets)
    else:
        s_targets = 1.0
        savetargets = targets

    return targets, s_targets, savetargets, sd


def sstargets(fdata, pc):
    # Determine steady-state target values via mass balance
    # 0 = (flow/volume) * (Conc_out - Conc_in)
    n, s = np.shape(fdata)
    targets = np.zeros([n, s])
    for i in range(n):
        for j in range(s):
            targets[i, j] = (pc["flow"][i][j] / pc["vol"][i]) * (
                fdata[i, j] - pc["x0"][i][j]
            )
    return targets


def dynamictargets(data, fdata, pc, sharedata):
    # import sympy
    from sympy import symbols, diff

    n, spec = np.shape(fdata)
    dshape = np.shape(data)
    ntimedata = dshape[0]
    targets = np.ones([n, spec])
    tempn = n
    st = symbols("st")

    while tempn >= ntimedata:
        [prof, r2, sharedata] = doalamo(
            data[tempn - ntimedata : tempn, :, :],
            pc["t"][tempn - ntimedata : tempn],
            1,
            ntimedata,
            range(int(pc["npc"])),
            range(int(spec)),
            sharedata,
        )
        for l in list(["profiles", "tlo", "tup", "etrack"]):
            for prc in range(pc["npc"]):
                sharedata["pwalamo"][l][prc] = {}
                for s in range(spec):
                    sharedata["pwalamo"][l][prc][s] = list()
        sharedata = dopwalamo(
            data[tempn - ntimedata : tempn, :, :],
            pc["t"][tempn - ntimedata : tempn],
            prof,
            r2,
            1,
            ntimedata,
            range(pc["npc"]),
            range(spec),
            sharedata,
        )

        for p in range(spec):
            for i in range(pc["npc"]):
                pwt = 0
                for j in range(tempn - ntimedata, tempn):
                    profile = sharedata["pwalamo"]["profiles"][i][p][pwt]
                    targets[j, p] = diff(profile, st).subs(st, pc["t"][j])

                    if j == sharedata["pwalamo"]["tup"][i][p][pwt]:
                        pwt = pwt + 1

        tempn -= ntimedata
    return targets, sharedata
