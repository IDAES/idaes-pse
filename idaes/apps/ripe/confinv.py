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


def confinv(aterm, targets, results, ccon, ptype, pc, sharedata, sigma):
    # This subroutine calculates 95% confidence intervals for estimated ripe models
    # inputs:
    # aterm     - activity array
    # targets   - observed rates of generation
    # results   - selected ripe model
    # ccon      - selected cardinality
    # ptype     - flag for isothermal/arrhenious
    # pc        - process condition array
    # sharedata - share data array defined in shared.py
    # sigma     - estimated/defined residual variance
    # Outputs:
    # kres      - preexponential factor estimate/confidence
    # eres      - activation eenrgy " "
    # stoichres - selected stoichiometry
    # mechres   - selected mechanisms
    # r2        - R2 on training data

    # Alpha value defined 95% CI, can be made variable with ease
    alpha = 0.05
    import numpy as np
    import idaes.apps.ripe as ripe

    # from scipy.optimize import curve_fit
    from scipy.stats.distributions import t

    n, ns, nm, nh = np.shape(aterm)
    # Reshape targets and aterm
    flat_targets = np.ndarray.flatten(targets, order="F")
    flat_sigma = np.ndarray.flatten(sigma, order="F")
    # Tdata needs to be appended to aterm in the case of arrhenious
    tempaterm = np.zeros([n, ns, ccon])
    idx = 0
    for res in results["ind"]:
        tempaterm[:, :, idx] = aterm[:, :, int(res[0]) - 1, int(res[1]) - 1]
        idx += 1
    aterm = tempaterm

    # collect results
    p0 = []
    kres = []
    eres = []
    stoichres = []
    mechres = []
    # kinformres = []
    for inds in results["ind"]:
        i1 = inds[0]
        i2 = inds[1]
        mechres.append(i1)
        stoichres.append(i2)
        for k in results["k"]:
            if k[1][0] == i1 and k[1][1] == i2:
                # unwind scaling from model selection here
                #                tempk = s_target * (k[0] / scales['nfactor'][int(i1)-1][int(i2)-1])
                tempk = k[0]
                kres.append(tempk)
                p0.append(tempk)
    if "E" in results.keys():
        for inds in results["ind"]:
            i1 = inds[0]
            i2 = inds[1]
            for e in results["E"]:
                if e[1][0] == i1 and e[1][1] == i2:
                    eres.append(e[0])
                    p0.append(e[0])

    # analyze what kinetic model is used to determine models to include in conf_inv calculation
    st = np.average(flat_targets)
    nf = len(flat_targets)
    p = len(p0)
    dof = max(1, int(nf / ns) - p)
    tval = t.ppf(1.0 - alpha / 2.0, dof)

    # The file kinforms.py contains definitions for types kinetic models
    # defaults include: isothermal, arrhenious, arrhenious-reformulation
    if ptype == "arr":
        if "Tref" not in pc.keys():
            sr = ripe.kinforms.arr(
                targets, np.ndarray.flatten(pc["T"]), aterm, p0, sharedata["gasconst"]
            )  # ,sigma)
            sensmat = ripe.kinforms.arrjac(
                targets, np.ndarray.flatten(pc["T"]), aterm, p0, sharedata["gasconst"]
            )
        if "Tref" in pc.keys():
            sr = ripe.kinforms.refarr(
                targets,
                np.ndarray.flatten(pc["T"]),
                sharedata["Tref"],
                aterm,
                p0,
                sharedata["gasconst"],
            )  # ,sigma)
            sensmat = ripe.kinforms.refarrjac(
                targets,
                np.ndarray.flatten(pc["T"]),
                sharedata["Tref"],
                aterm,
                p0,
                sharedata["gasconst"],
            )
    else:
        # Calculate squared residuals and sensitivity matrix
        sr = ripe.kinforms.lin(targets, aterm, p0)
        sensmat = ripe.kinforms.linjac(targets, aterm, p0)

    sst = sum((flat_targets - st) ** 2)  # ,flat_sigma))
    r2 = 1.0 - sum(sr) / sst

    # Calculate confidence intervals
    stmat = np.matmul(np.ndarray.transpose(sensmat), sensmat)
    mmat = np.matmul(
        np.matmul(np.ndarray.transpose(sensmat), np.diag(flat_sigma)), sensmat
    )
    invmat = np.linalg.inv(stmat)
    covar = np.multiply(np.multiply(invmat, mmat), invmat)
    confinvs = [tval * (s**0.5) for s in np.diag(covar)]
    itemp = 0
    if ptype == "arr":
        for i in range(int(0.5 * len(p0))):
            kres[itemp] = [kres[itemp], confinvs[i]]
            eres[itemp] = [eres[itemp], confinvs[int(0.5 * len(p0)) + i]]
            itemp += 1
    else:
        for i in range(len(p0)):
            kres[i] = [kres[i], confinvs[i]]

    return [kres, eres, stoichres, mechres, r2]
