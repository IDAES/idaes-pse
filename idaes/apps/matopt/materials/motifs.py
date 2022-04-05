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
from copy import deepcopy
from ..util.util import myPointEq


def areMotifViaTransF(D1, D2, TransF, blnPreserveIndexing=False, blnIgnoreVoid=True):
    """

    Args:
        D1: param D2:
        TransF: param blnPreserveIndexing:  (Default value = False)
        blnIgnoreVoid: Default value = True)
        D2:
        blnPreserveIndexing:  (Default value = False)

    Returns:

    """
    return D1.isEquivalentTo(
        D2.getTransformed(TransF),
        blnPreserveIndexing=blnPreserveIndexing,
        blnIgnoreVoid=blnIgnoreVoid,
    )


def areMotifViaTransFs(D1, D2, TransFs, blnPreserveIndexing=False, blnIgnoreVoid=True):
    """

    Args:
        D1: param D2:
        TransFs: param blnPreserveIndexing:  (Default value = False)
        blnIgnoreVoid: Default value = True)
        D2:
        blnPreserveIndexing:  (Default value = False)

    Returns:

    """
    for TransF in TransFs:
        if areMotifViaTransF(
            D1,
            D2,
            TransF,
            blnPreserveIndexing=blnPreserveIndexing,
            blnIgnoreVoid=blnIgnoreVoid,
        ):
            return True
    return False


def getEnumConfs(
    Motifs, TransFs, MotifToConfMap=None, blnPreserveMotifLocs=True, DBL_TOL=1e-5
):
    """

    Args:
        DBL_TOL:
        blnPreserveMotifLocs:
        Motifs: param TransFs:
        MotifToConfMap: Default value = None)
        TransFs:

    Returns:

    """
    if MotifToConfMap is not None:
        # del MotifToConfMap[:]
        MotifToConfMap.clear()
    result = []
    for iMotif, Motif in enumerate(Motifs):
        if MotifToConfMap is not None:
            MotifToConfMap.append([len(result)])
        result.append(Motif)
        for TransF in TransFs:
            Conf = Motif.getTransformed(TransF)
            blnUniqueConf = True
            for r in result:
                if Conf.isEquivalentTo(r):
                    blnUniqueConf = False
                    break
            if blnUniqueConf:
                if MotifToConfMap is not None:
                    MotifToConfMap[iMotif].append(len(result))
                # NOTE: This section is important if we want to just assume
                #       that each index is the same location.
                #       This is important for use in models.
                if blnPreserveMotifLocs:
                    ConfCopy = Conf
                    Conf = deepcopy(Motif)  # So they have the same points<->indexes
                    for l in range(len(Motif)):
                        for l2 in range(len(ConfCopy)):
                            if myPointEq(
                                Motif.Canvas.Points[l],
                                ConfCopy.Canvas.Points[l2],
                                DBL_TOL,
                            ):
                                Conf.setContent(l, ConfCopy.Contents[l2])
                                break
                result.append(Conf)
    return result
