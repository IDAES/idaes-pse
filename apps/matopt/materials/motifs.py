
from copy import deepcopy
import os
from matopt.materials.design import Design,loadFromPDBs
from matopt.util.util import myPointEq

def areMotifViaTransF(D1,D2,TransF,
                      blnPreserveIndexing=False,
                      blnIgnoreVoid=True):
    """

    Args:
        D1: param D2:
        TransF: param blnPreserveIndexing:  (Default value = False)
        blnIgnoreVoid: Default value = True)
        D2: 
        blnPreserveIndexing:  (Default value = False)

    Returns:

    """
    return D1.isEquivalentTo(D2.getTransformed(TransF),
                             blnPreserveIndexing=blnPreserveIndexing,
                             blnIgnoreVoid=blnIgnoreVoid)

def areMotifViaTransFs(D1,D2,TransFs,
                       blnPreserveIndexing=False,
                       blnIgnoreVoid=True):
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
        if(areMotifViaTransF(D1,D2,TransF,
                             blnPreserveIndexing=blnPreserveIndexing,
                             blnIgnoreVoid=blnIgnoreVoid)):
            return True
    return False

def getEnumConfs(Motifs,TransFs,MotifToConfMap=None,blnPreserveMotifLocs=True,DBL_TOL=1e-5):
    """

    Args:
        Motifs: param TransFs:
        MotifToConfMap: Default value = None)
        TransFs: 

    Returns:

    """
    if(MotifToConfMap is not None):
        #del MotifToConfMap[:]
        MotifToConfMap.clear()
    result = []
    for iMotif,Motif in enumerate(Motifs):
        if(MotifToConfMap is not None): 
            MotifToConfMap.append([len(result)])
        result.append(Motif)
        for TransF in TransFs:
            Conf = Motif.getTransformed(TransF)
            blnUniqueConf = True
            for r in result:
                if(Conf.isEquivalentTo(r)):
                    blnUniqueConf = False
                    break
            if(blnUniqueConf):
                if(MotifToConfMap is not None):
                    MotifToConfMap[iMotif].append(len(result))
                # NOTE: This section is important if we want to just assume 
                #       that each index is the same location.
                #       This is important for use in models. 
                if(blnPreserveMotifLocs):
                    ConfCopy = Conf
                    Conf = deepcopy(Motif) # So they have the same points<->indexes
                    for l in range(len(Motif)):
                        for l2 in range(len(ConfCopy)):
                            if(myPointEq(Motif.Canvas.Points[l],ConfCopy.Canvas.Points[l2],DBL_TOL)):
                                Conf.setContent(l,ConfCopy.Contents[l2])
                                break
                result.append(Conf)
    return result

