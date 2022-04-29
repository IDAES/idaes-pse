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
# This file contians subroutines related to finding the bounds for big M constraints
# additionally, depedent reaction stoichiometries are detected here
from idaes.apps import ripe
import numpy as np


def stoich_cons(stoichs):
    # This subroutine identifies dependent stoichiometries
    smat = np.array(stoichs)
    nt, ns = np.shape(stoichs)
    scons_small = []
    scons_large = []
    fam_cons = []
    nfcons = []
    ncons = []
    for con in [3, 4]:
        for ins in nt.permutations(range(nt), con):  # was it
            inds = list(ins)
            tmat = smat[inds, :]
            rank = np.linalg.matrix_rank(tmat)
            if rank < con:
                msize = count_neg(tmat)
                if all(msize[0] == item for item in msize):
                    fam_cons.append([inds[i] + 1 for i in inds])
                    nfcons.append(con - 1)
                else:
                    bigones = np.argwhere(msize == np.max(msize)).flatten().tolist()
                    smallones = np.argwhere(msize != np.max(msize)).flatten().tolist()
                    scons_small.append([inds[i] + 1 for i in smallones])
                    scons_large.append([inds[i] + 1 for i in bigones])
                    ncons.append(con - 1)
    return scons_small, scons_large, ncons, fam_cons, nfcons


def count_neg(mat):
    r = []
    for i in range(len(mat)):
        vec = mat[i]
        count = 0
        for item in vec:
            # This line is commented out to enforce total reaction order (not just reactants)
            #            if item < 0:
            count += abs(item)
        r.append(count)
    return r


def get_bounds(inargs):
    sharedata = inargs[-1]
    inargs[-1]["maxmiptime"] = 60.0
    pc = inargs[-2]
    # try except commented out
    #    try:
    res = ripe.genpyomo.ripeomo(*inargs)
    #    except:
    #        res = {}
    #        res['maxk'] = sharedata['bounds']['k']['max']
    #        res['maxe'] = sharedata['bounds']['e']['max']
    # Calculate upperbound of E from closed form expression
    # c_1,2 = (-1/(R*T_min,max))
    # E_max = ln(c2/c1)/(c1-c2)
    T = [pc["T"][i][0] for i in range(pc["npc"])]
    Tmax = np.max(T)
    Tmin = np.min(T)
    c1 = -1 / (sharedata["gasconst"] * Tmin)
    c2 = -1 / (sharedata["gasconst"] * Tmax)
    arrpen = c2 - c1
    # sharedata['bounds']['e']['max'] = np.log(c2/c1)/(c1-c2)
    kmax = res["maxk"]
    # emax = res["maxe"]
    if kmax >= sharedata["bounds"]["k"]["max"] or kmax == 0.0:
        kmax = 10 * sharedata["bounds"]["k"]["max"]
    elif kmax == 0.0:
        kmax = 10 * sharedata["bounds"]["k"]["max"]
    else:
        kmax = 10 * kmax
    #    if emax >= sharedata['bounds']['e']['max'] or emax == 0.0:
    #        emax = 5*sharedata['bounds']['e']['max']
    #    else:
    #        emax = 1000*emax
    sharedata["bounds"]["k"]["max"] = kmax
    sharedata["arrpen"] = arrpen
    return sharedata
