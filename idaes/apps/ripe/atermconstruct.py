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
import numpy as np
from idaes.apps.ripe import mechs
from .shared import sharedata as sd

# from ripe.shared import debug as debug

consivars = sd["ivars"]


def makeaterm(data, stoich, rxn_mechs, kwargs, ncons, mechlist, fixarray, sharedata):
    # This subroutine constructs an activity matrix that is used to encode
    # considered stoichimetry and mechanisms
    # Input :
    # data      - process data from ripemodel()
    # stoichs   - considered reactions toichiometries
    # mechs     - considered reaction mechanisms
    # kwargs    - keyword arguments from ripemodel()
    # mechlist  - list containing which mechs correspond to which stoich
    # fixarray  - [#mech by #stoich] matrix that encodes valid combinations
    # Outputs :
    # aterm     - activity matrix
    # fdata     - formatted data for use in writing pyomo models
    # pc        - formatted process conditions
    # data      - data is reshaped and returned
    # scales    -  dictionary containing relevant information for scaling

    # Call subroutine to analyze kwargs and format data
    data, kwargs, fdata, pc, alld = formatinputs(data, kwargs)

    # Find sizes and initialize aterm
    ndata, ns = np.shape(fdata)
    nm = len(mechlist)
    nh = len(stoich)
    aterm = np.zeros([ndata, ns, nm, nh])

    # This section of code loops thorugh the aterm array
    # Callling s_mech in order to determine numerical values
    # mechanisms can be specified with variable stoichiometry
    for i1 in range(ndata):
        # index i1 over observations
        for i2 in range(ns):
            # index i2 over species
            i3 = 0
            # index i3 over mechanisms
            for tempi in range(len(rxn_mechs)):
                mechline = rxn_mechs[tempi]
                for mspec in mechline[1]:
                    for h_ind in list(mechline[0]):
                        # final index over stoichiometries
                        if mspec == "massact" or mechline[2]:
                            s_mech = mechs.mechperstoich(mspec, stoich[h_ind])
                        else:
                            s_mech = mspec
                        if "T" in kwargs.keys():
                            inv = np.hstack((fdata[i1, :], pc["T"][i1][0]))
                        else:
                            inv = fdata[i1, :]
                        atemp = stoich[h_ind][i2] * s_mech(*inv)
                        aterm[i1, i2, i3, h_ind] = atemp
                    i3 += 1

    # Scale data if specified
    if sharedata["ascale"]:
        aterm, scales = normalizefeatures(aterm, fixarray)
    else:
        scales = []

    return [aterm, fdata, pc, data, scales]


def formatinputs(data, kwargs):
    # This subroutine formats inputs supplied to ripemodel()
    # Inputs:
    # data     - Input data that contains variable size and type
    # kwargs   - key word arguments supplied to ripemodel()
    # Outputs:
    # data     - Format and shape of input data is made consistent
    # kwargs   - kwargs is checked and returned
    # fdata    - formatted input data
    # pc       - process conditions are puleld from kwargs and formatted
    # alldata  - another reformatting of existing data required by ripeems()

    # input data may be 1,2, or 3 dimensional at this point

    # Define what kwargs are used
    inkeys = list(kwargs)  # kwargs.keys() python 2
    # Define what keys are not used, dynamic data is handled differently so t is removed
    inkeys_not = list(set(inkeys) - set(["t", "other"]))
    # This line ensures additional kwargs are not analyzed/can be defined
    notinkeys = list(set(consivars) - set(inkeys))

    # initialize process condition dictionary and find shape of data provided
    pc = {}
    dshape = np.shape(data)

    if "t" in inkeys:
        # dictionaries do not perserve order - ensure 't' is at end
        inkeys.remove("t")
        inkeys.append("t")
        # if data is 1 dimensional and dynamic then no process conditions should be specified
        if len(dshape) == 1:
            data = np.expand_dims(np.expand_dims(data, axis=-1), axis=-1)
            ns = 1
            nobs = dshape[0]
            npc = 1
        elif len(dshape) == 2:
            data = np.expand_dims(data, axis=-1)
            # first test to see if other process conditions exist
            if inkeys == inkeys_not:
                # ripe problem is dynamic with 1 process condition
                ns = dshape[1]
                nobs = dshape[0]
                npc = 1
            else:
                # data must have only 1 species
                ns = 1
                npc = dshape[1]
                nobs = dshape[0]
        else:
            nobs = dshape[0]
            ns = dshape[1]
            npc = dshape[2]
    else:
        # nobs is used to track observations per process condition without time this is 1
        nobs = 1
        if len(dshape) == 1:
            # data has one species
            data = np.expand_dims(np.expand_dims(data, axis=-1), axis=-1)
            npc = dshape[0]
            ns = 1
        else:
            # third dimension is redudant be needed for future code
            data = np.expand_dims(np.expand_dims(data, axis=-1), axis=-1)
            ns = dshape[1]
            npc = dshape[0]

    # ensure consistency in process condition formats
    kwargs = checkargs(kwargs, inkeys, npc, nobs, ns)

    # ndata is always npc*nobs
    ndata = npc * nobs

    # initialize formatted data array
    fdata = np.zeros([ndata, ns])

    # Construct process condition array
    for key in inkeys_not:
        pc[key] = kwargs[key]

    # construction of fdata makes construction of aterm much easier
    ind = 0
    if "t" in inkeys:
        pc["t"] = np.ndarray.flatten(kwargs["t"])
        if npc == 1 and ns == 1:
            for i in range(nobs):
                fdata[i, :] = data[i, 0, 0]
        elif npc == 1:
            for i in range(nobs):
                for j in range(ns):
                    fdata[i, j] = data[i, j, 0]
        elif ns == 1:
            while ind < ndata:
                for i in range(nobs):
                    for j in range(npc):
                        fdata[ind, 0] = data[i, j, 0]
                        ind += 1
        else:
            while ind < ndata:
                for k in range(npc):
                    for i in range(nobs):
                        for j in range(ns):
                            fdata[ind, j] = data[i, j, k]
                        ind += 1
    else:
        if ns == 1:
            for i in range(npc):
                fdata[i, 0] = data[i, 0, 0]
        else:
            for i in range(npc):
                for j in range(ns):
                    fdata[i, j] = data[i, j, 0]

    # Initialize notinkeys for easy code in pyomo models
    for key in notinkeys:
        if key == "vol":
            pc["vol"] = [1.0] * ndata
        elif key == "flow":
            pc["flow"] = [[1.0] * ns] * ndata
        elif key == "x0":
            pc["x0"] = [[0.0] * ns] * ndata
        elif key == "T":
            pc["T"] = [[-1.0]] * ndata

    # define npc and nobs for easy access later
    #  additional process-specific information could
    #  be define in the dictionary pc without issue
    pc["npc"] = npc
    pc["nobs"] = nobs

    if "Tr" in kwargs.keys():
        pc["Tref"] = kwargs["Tr"][0][0]

    # return concatenated array for use in returned models in ems
    num_other = 0
    alldata = np.zeros([ndata, 3 * ns + 2 + num_other])
    for i in range(ndata):
        alldata[i, :ns] = fdata[i, :]
        alldata[i, ns : 2 * ns] = pc["x0"][i]
        alldata[i, 2 * ns] = pc["T"][i][0]
        alldata[i, 2 * ns + 1 : 3 * ns + 1] = pc["flow"][:][i]
        alldata[i, 3 * ns + 1] = pc["vol"][i]

    # append other process conditions
    if "other" in pc.keys():
        alldata = np.concatenate((alldata, pc["other"]), axis=0)

    return [data, kwargs, fdata, pc, alldata]


def checkargs(kwargs, inkeys, npc, nobs, ns):
    # This subroutine checks and formats kwargs in order to construct a process condition array
    # Inputs:
    # kwargs  - kwargs from ripemodel()
    # inkeys  - kwargs specified
    # npc     - number of process conditions
    # nobs    - number of observations
    # ns      - number of species
    # Outputs:
    # kwargs - reformatted kwargs

    for key in list(set(inkeys) - set(["t", "Tref"])):
        # reformat kwargs so it is easier to process
        if isinstance(kwargs[key], type(0.0)):
            kwargs[key] = [[kwargs[key]]] * npc
        elif isinstance(kwargs[key], type(0)):
            kwargs[key] = [[float(kwargs[key])]] * nobs
        elif key in ["x0", "flow"]:
            if len(np.shape(kwargs[key])) == 1 and ns != 1:
                kwargs[key] = [kwargs[key]] * npc
            elif len(np.shape(kwargs[key])) == 1 and ns == 1:
                kwargs[key] = np.expand_dims(np.asarray(kwargs[key]), axis=-1)
        else:
            kwargs[key] = np.expand_dims(np.asarray(kwargs[key]), axis=-1)
    if "t" in kwargs:
        if len(np.shape(kwargs["t"])) == 1 and npc != 1:
            kwargs["t"] = [kwargs["t"]] * npc

    return kwargs


def normalizefeatures(aterm, fixarray):
    # This subroutine calculate scaling values for an activity array
    # Inputs:
    # aterm      - activity array
    # fixarray   - array containing mech/stoich pairs
    # Outputs:
    # aterm      - scaled aterm
    # scale_dict - information for unwinding scales

    nd, ns, nm, nh = np.shape(aterm)

    scale_dict = {}
    scale_dict["aterm_u"] = aterm
    n_factor = np.zeros([nm, nh])
    for k in range(nm):
        for l in range(nh):  # L
            if fixarray[l, k] == 1:
                maxval = np.max(
                    [x for x in np.ndarray.flatten(aterm[:, :, k, l])]
                )  # if x > 0.0])
                minval = np.min(
                    [x for x in np.ndarray.flatten(aterm[:, :, k, l])]
                )  # if x > 0.0])
                n_factor[k, l] = np.max([np.abs(maxval), np.abs(minval)])
                for i in range(nd):
                    for j in range(ns):
                        aterm[i, j, k, l] = aterm[i, j, k, l] / n_factor[k, l]
    scale_dict["nfactor"] = n_factor
    return aterm, scale_dict
