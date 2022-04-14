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
"""
Run ALAMO .
"""
import sys
import collections
import numpy as np
import os

from idaes.apps import alamopy_depr as alamopy
from idaes.apps.alamopy_depr import almerror
from idaes.apps.alamopy_depr.multos import deletefile, has_alamo


def doalamo(xdata, zdata, **kwargs):
    """
    Warning: doalamo is deprecated.
    please use ``alamopy.doalamo.alamo(xdata, zdata, **kwargs)``
    """
    return alamo(xdata, zdata, **kwargs)


def alamo(xdata, zdata, **kwargs):
    """
    [almmodel] = doalamo(xdata,zdata, xvaldata, zvaldata,addopt=vals)

    Args:
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real)
        kwargs: Additional options may be specified and will be applied
                to the .alm. Example -  monomialpower=(1,2,3,4)

          -  xlabels       : labels given to input variables
          -  zlabels       : labels given to outputs
          -  xval          : validaiton data for alamo
          -  zval          : response validation data for alamo
          -  modeler       : modeler value used in alamo
          -  solvemip      : force alamo to solve mip if gams is availible
          -  linfcns       : 0-1 option to include linear transformations
          -  expfcns       : 0-1 option to include exponential transformations
          -  logfcns       : 0-1 option to include logarithmic transformations
          -  sinfcns       : 0-1 option to include sine transformations
          -  cosfcns       : 0-1 option to include cosine transformations
          -  monomialpower : list of monomial powers
          -  multi2power   : list of binomial powers
          -  multi3power   : list of trinomials
          -  ratiopower    : list of ratio powers
          -  screener      : screening method
          -  almname       : specify a name for the .alm file
          -  savescratch   : saves .alm and .lst
          -  savetrace     : saves trace file
          -  expandoutput  : add a key to the output dictionary for the output
                             (must be on for inputs(outputs?#Engle)>1)
          -  almopt        : direct text appending
                             the option almopt=<file> will append a file to the
                             end of the .alm and can be used to facilitate
                             direct access to the .alm (no current checks)
          -  loo           : leave one out evaluation
          -  lmo           : leave many out evaluation

    Returns:
        dict: An ALAMO model with the following keys

          -  'model'    : algebraic form of model
          -  'f(model)' : a callable lambda function
          -   Syntax is depended on expandout
               syntax => almmodel['f(model)']['out'](inputs,sep,by,comma)
                         almmodel['f(model)'](inputs,sep,by,comma)
          -  'ssr'      : SSE on training set provided
          -  'R2'       : R2 on training set provided
          -  'ssrval'   : SSE on testing set if provided
          -  'R2val'    : R2 on testing set if provided

    """
    alamopy.initialize()
    data, debug = alamopy.data, alamopy.debug

    # patched together validation data check
    if kwargs.get("xval", None) is not None and kwargs.get("zval", None) is not None:
        vargs = (kwargs["xval"], kwargs["zval"])
    else:
        vargs = ()

    xdata, zdata, xvaldata, zvaldata = setupData(
        data, debug, xdata, zdata, vargs, kwargs
    )
    manageArguments(xdata, zdata, data, debug, kwargs)

    data["results"] = {}

    writeCustomALAMOOptions(kwargs)  # New Custom Options MENGLE

    # Cross Validation
    if debug["loo"]:
        q2 = []
        if debug["outkeys"] or debug["expandoutput"]:
            q2 = {}

        data["opts"]["ndata"] = data["opts"]["ndata"] - 1
        kwargValidation = debug["validation"]
        kwargSaveTrace = debug["savetrace"]
        if kwargValidation:
            kwargNvaldata = data["opts"]["nvaldata"]

        data["opts"]["nvaldata"] = 1
        debug["validation"] = True
        debug["savetrace"] = False
        for i in range(0, len(xdata)):
            cvxdata = [x for y, x in enumerate(xdata) if y != i]
            cvzdata = [x for y, x in enumerate(zdata) if y != i]
            alamopy.almwriter(
                data, debug, (cvxdata, cvzdata, [xdata[i][:]], [zdata[i][:]]), kwargs
            )

            # Calling ALAMO
            if not debug["mock"]:
                os.system(
                    debug["almloc"]
                    + " "
                    + str(data["stropts"]["almname"])
                    + " > logscratch"
                )

            data["results"] = {}
            readTraceFile([xdata[i][:], zdata[i][:]], data, debug)

            if debug["outkeys"] or debug["expandoutput"]:
                for k in data["results"]["R2"].keys():
                    if k not in q2.keys():
                        q2[k] = [float(data["results"]["R2val"][k])]
                    else:
                        q2sub = q2[k]
                        q2sub.append(float(data["results"]["R2val"][k]))
                        q2[k] = q2sub
            else:
                q2.append(float(data["results"]["R2val"]))
            cleanFiles(data, debug)

        if debug["outkeys"] or debug["expandoutput"]:
            data["results"]["Q2"] = {}
            for k in q2.keys():
                Q2 = np.mean(q2[k])
                data["results"]["Q2"][k] = Q2
                print("%s: Running cross validation LOO, evaluated Q2:%f" % (k, Q2))
        else:
            Q2 = np.mean(q2)
            data["results"]["Q2"] = Q2
            print("Running cross validation LOO, evaluated Q2:%f" % Q2)

        del data["opts"]["nvaldata"]
        debug["validation"] = kwargValidation
        debug["savetrace"] = kwargSaveTrace
        if kwargValidation:
            data["opts"]["nvaldata"] = kwargNvaldata
        data["opts"]["ndata"] = data["opts"]["ndata"] + 1
    elif debug["lmo"] > 0:
        q2 = []
        if debug["outkeys"] or debug["expandoutput"]:
            q2 = {}

        kwargNdata = data["opts"]["ndata"]
        kwargValidation = debug["validation"]
        kwargSaveTrace = debug["savetrace"]
        if kwargValidation:
            kwargNvaldata = data["opts"]["nvaldata"]

        debug["validation"] = True
        debug["savetrace"] = False
        numOfFolds = debug["lmo"]

        if numOfFolds > len(xdata):
            raise Exception(
                "Number of Cross validation \
                            folds exceeds the number of data points %i"
                % debug["lmo"]
            )

        sizeOfFolds = int(len(xdata) / numOfFolds)
        r = len(xdata) % numOfFolds
        remS = 0
        remE = 1

        for i in range(numOfFolds):
            if i < r + 1:
                remS = i
                remE = i + 1
            cvvalxdata = xdata[remS + sizeOfFolds * i : sizeOfFolds * (i + 1) + remE]
            cvvalzdata = zdata[remS + sizeOfFolds * i : sizeOfFolds * (i + 1) + remE]
            if i == 0:
                cvxdata = xdata[sizeOfFolds * (i + 1) + remE :]
                cvzdata = zdata[sizeOfFolds * (i + 1) + remE :]
            else:
                cvxdata = np.concatenate(
                    [
                        xdata[0 : remS + sizeOfFolds * i],
                        xdata[sizeOfFolds * (i + 1) + remE :],
                    ]
                )
                cvzdata = np.concatenate(
                    [
                        zdata[0 : remS + sizeOfFolds * i],
                        zdata[sizeOfFolds * (i + 1) + remE :],
                    ]
                )

            data["opts"]["nvaldata"] = len(cvvalxdata)
            data["opts"]["ndata"] = len(cvxdata)

            alamopy.almwriter(
                data, debug, (cvxdata, cvzdata, cvvalxdata, cvvalzdata), kwargs
            )

            # Calling ALAMO
            if not debug["mock"]:
                os.system(
                    debug["almloc"]
                    + " "
                    + str(data["stropts"]["almname"])
                    + " > logscratch"
                )

            data["results"] = {}
            expandOutput(xdata, zdata, [cvvalxdata, cvvalzdata], data, debug)
            readTraceFile([cvvalxdata, cvvalzdata], data, debug)

            if debug["outkeys"] or debug["expandoutput"]:
                for k in data["results"]["R2"].keys():
                    if k not in q2.keys():
                        q2[k] = [float(data["results"]["R2val"][k])]
                else:
                    q2sub = q2[k]
                    q2sub.append(float(data["results"]["R2val"][k]))
                    q2[k] = q2sub
            else:
                q2.append(float(data["results"]["R2val"]))
            cleanFiles(data, debug)

        if debug["outkeys"] or debug["expandoutput"]:
            data["results"]["Q2"] = {}
            for k in q2.keys():
                Q2 = np.mean(q2[k])
                data["results"]["Q2"][k] = Q2
                print("%s: Running cross validation LMO, evaluated Q2:%f" % (k, Q2))
        else:
            Q2 = np.mean(q2)
            data["results"]["Q2"] = Q2
            print("Running cross validation LMO, evaluated Q2:%f" % Q2)

        del data["opts"]["nvaldata"]
        debug["validation"] = kwargValidation
        debug["savetrace"] = kwargSaveTrace
        if kwargValidation:
            data["opts"]["nvaldata"] = kwargNvaldata
        data["opts"]["ndata"] = kwargNdata

    # Write alamo file
    if debug["validation"]:
        alamopy.almwriter(data, debug, (xdata, zdata, xvaldata, zvaldata), kwargs)
    else:
        alamopy.almwriter(data, debug, (xdata, zdata), kwargs)

    # Call alamo from the terminal
    if not debug["mock"]:
        if debug["showalm"]:
            os.system(debug["almloc"] + " " + str(data["stropts"]["almname"]))
        else:
            os.system(
                debug["almloc"]
                + " "
                + str(data["stropts"]["almname"])
                + " > logscratch"
            )

    # Check to see if additional data was sampled and add it
    if kwargs.get("simulator", None) is not None:
        xdata, zdata = checkForSampledData(data, debug)

    # calculate additional statistics
    expandOutput(xdata, zdata, vargs, data, debug)

    # Open the trace file and pull appropriate results
    readTraceFile(vargs, data, debug)

    # write python file of regressed model
    alamopy.almpywriter(data, debug)
    if debug["cvfun"]:
        alamopy.almcvwriter(data)

    # add <>alm.py to results dict
    data["results"]["pymodel"] = data["stropts"]["almname"].split(".")[0] + "alm"

    if debug["loo"] or debug["lmo"] > 0:
        Q2, R2, diff = 0, 0, 0
        if debug["outkeys"] or debug["expandoutputs"]:
            for k in data["results"]["R2"].keys():
                R2 = data["results"]["R2"][k]
                if k in data["results"]["Q2"].keys():
                    Q2 = data["results"]["Q2"][k]
                    diff = float(R2) - float(Q2)
                    if Q2 < 0.5:
                        print(
                            "%s: Q2 suggests this is not a predictive model, \
                              Q2: %f, R2: %s"
                            % (k, Q2, R2)
                        )
                    elif diff < 0.3:
                        print(
                            "%s: The difference of R2-Q2 is  %f. This is an acceptable \
                              difference for predictability, Q2: %f, R2: %s"
                            % (k, diff, Q2, R2)
                        )
                    else:
                        print(
                            "%s: The difference of R2-Q2 is %f. The surrogate model is \
                              not able to predict the data reliably, Q2: %f, R2: %s"
                            % (k, diff, Q2, R2)
                        )

        else:
            R2 = data["results"]["R2"]
            Q2 = data["results"]["Q2"]
            diff = float(R2) - float(Q2)

            if Q2 < 0.5:
                print(
                    "Q2 suggests this is not a predictive model, Q2: %f, R2: %s"
                    % (Q2, R2)
                )
            elif diff < 0.3:
                print(
                    "The difference of R2-Q2 is  %f. This is an acceptable difference for \
                      predictability, Q2: %f, R2: %s"
                    % (diff, Q2, R2)
                )
            else:
                print(
                    "The difference of R2-Q2 is %f. The surrogate model is not able to \
                      predict the data reliably, Q2: %f, R2: %s"
                    % (diff, Q2, R2)
                )

    cleanFiles(data, debug, pywrite=True, **kwargs)

    return data["results"]


# Data Management


def setupData(data, debug, xdata, zdata, vargs, kwargs):
    """
    [xdata, zdata, xvaldata, zvaldata] =
    setupData(data, debug, xdata,zdata, vargs, kwargs)

    Checks inputted data and resturctures the data for the .alm file

    Args:
        data/debug: shared default options for .alm file
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real)
        vargs:  Validation data
        kwargs: Additional options may be specified and will be applied
                  to the .alm
    """

    xvaldata, zvaldata = list(), list()
    checkinput(data, debug, xdata, zdata, vargs, kwargs)
    xdata, zdata = getTrainingData(xdata, zdata, data, debug)
    if len(vargs) > 0:
        xvaldata, zvaldata = getValidationData(vargs, data, debug)
    else:
        debug["validation"] = False
    getlabels(data, debug, kwargs)

    return xdata, zdata, xvaldata, zvaldata


def checkinput(data, debug, xdata, zdata, vargs, kwargs):
    """
    Check the input data into doalamo for errors.

    Args:
        data/debug: shared default options for .alm file
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real)
        vargs:  Validation data
        kwargs: Additional options may be specified and will be applied
                to the .alm
    """

    # t = np.zeros([1, 1])
    kk = kwargs.keys()
    # First check vargs for validation set viability
    if len(vargs) > 0 and len(vargs) != 2:
        raise almerror.AlamoInputError(
            "Validation Xdata and Zdata must be"
            "specified in the following form: "
            "alamopy.doalamo(x,z,xval,zval)"
        )

    for arg in kk:
        if arg not in data["pargs"]["opts"] + data["pargs"]["stropts"] + data["pargs"][
            "lstopts"
        ] + data["pargs"]["set4"] + debug["pargs"] + ["xlabels", "zlabels"] + [
            "xval",
            "zval",
        ]:
            raise almerror.AlamoInputError(
                "The following argument given to"
                "alamo() is not understood: {}".format(arg)
            )

    # read keys for debug
    for key in debug["pargs"]:
        if key in kk:
            debug[key] = kwargs[key]

    if kwargs.get("almname", None) is not None:
        data["stropts"]["almname"] += ".alm"


def getTrainingData(xdata, zdata, data, debug):
    """
    Structure data for training the model. Modifies data['opts']

    Args:
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real)
        data: shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
    """

    dshape = np.shape(xdata)
    if len(dshape) == 0:
        debug["traindata"] = False
    elif len(dshape) == 1:
        data["opts"]["ndata"] = np.shape(xdata)[0]
        data["opts"]["ninputs"] = 1
    else:
        data["opts"]["ndata"] = np.shape(xdata)[0]
        data["opts"]["ninputs"] = np.shape(xdata)[1]
    xdata = np.asarray(xdata)

    # Check training data
    if len(np.shape(zdata)) == 1:
        zdata = np.reshape(zdata, (data["opts"]["ndata"], 1))
        data["opts"]["noutputs"] = 1
    else:
        data["opts"]["noutputs"] = np.shape(zdata)[1]
    # PYLINT-TODO-FIX use almerror correctly instead of calling it directly
    if np.shape(zdata)[0] != data["opts"]["ndata"]:
        almerror("p1")  # pylint: disable=not-callable
    elif np.shape(xdata)[0] != data["opts"]["ndata"]:
        almerror("p1")  # pylint: disable=not-callable
    zdata = np.asarray(zdata)
    return xdata, zdata


def getValidationData(vargs, data, debug):
    """
    Structure data for validating the model. Modifies data['opts']

    Args:
        vargs: validation data valxdata, valzdata
        data: shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
    """

    if vargs != ():
        debug["validation"] = True
        xvaldata = vargs[0]
        zvaldata = vargs[1]
        temp = np.shape(xvaldata)
        data["opts"]["nvaldata"] = temp[0]
        if len(np.shape(zvaldata)) == 1:
            zvaldata = np.reshape(zvaldata, (data["opts"]["nvaldata"], 1))
        if temp[1] != data["opts"]["ninputs"]:
            raise almerror.AlamoInputError(
                "Number of input variables inconsistent between x and xval"
            )
        temp = np.shape(zvaldata)
        if temp[0] != data["opts"]["nvaldata"] or temp[1] != data["opts"]["noutputs"]:
            raise almerror.AlamoInputError("Problem with zval")
        return xvaldata, zvaldata


# Function/variable labels


def getlabels(data, debug, kwargs):
    """
    Creates labels for data and output. Modifies data['labs'].
    Makes labels if no labels are given.

    Args:
        data: shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
        vargs: validation data valxdata, valzdata
        kwargs: keyword arguments
    """

    # This function generates labels if they are not provided
    # Check to see if labels have been specified before we generate them
    if kwargs.get("xlabels", None) is not None:
        data["labs"]["savexlabels"] = kwargs["xlabels"]
    else:
        # Make savexlabels
        temp = list()
        for i in range(data["opts"]["ninputs"]):
            temp.append("x" + str(i + 1))
        data["labs"]["savexlabels"] = temp

    if kwargs.get("zlabels", None) is not None:
        data["labs"]["savezlabels"] = kwargs["zlabels"]
    else:
        # Make savezlabels
        temp = list()
        for i in range(data["opts"]["noutputs"]):
            temp.append("z" + str(i + 1))
        data["labs"]["savezlabels"] = temp

    # create temp labels for alm file
    # This avoids having to check user labels for alm rules
    data["labs"]["xlinks"] = [0 for i in range(data["opts"]["ninputs"])]
    data["labs"]["zlinks"] = [0 for i in range(data["opts"]["noutputs"])]
    makelabs(data, debug, "ninputs")
    makelabs(data, debug, "noutputs")


def makelabs(data, debug, param):
    """
    Constructs labels for alamo

    Args:
        data: shared alamo data options
        debug: Additional options may be specified and will be applied to the .alm
        param: 'ninputs' or 'noutputs'
    """

    temp = list([])
    for i in range(data["opts"][param]):
        if param == "ninputs":
            data["labs"]["xlinks"][i] = ["a", "b"]
            temp.append("almx" + str(i + 1) + "d")
            data["labs"]["xlinks"][i][0] = "almx" + str(i + 1) + "d"
            data["labs"]["xlinks"][i][1] = data["labs"]["savexlabels"][i]
            key = "xlabels"
        elif param == "noutputs":
            data["labs"]["zlinks"][i] = ["a", "b"]
            temp.append("almz" + str(i + 1) + "d")
            data["labs"]["zlinks"][i][0] = "almz" + str(i + 1) + "d"
            data["labs"]["zlinks"][i][1] = data["labs"]["savezlabels"][i]
            key = "zlabels"
    data["lstopts"][key] = temp


# Extra capabilities
group_list = []


def addBasisGroups(groups):
    """
    format: index_num  Type_of_function Member_indices(indice of the input) <pow>
    types: LIN, LOG, EXP, SIN, COS, MONO, MULTI2, MULTI3, RATIO, RBf, CUST, CONST
    Only one CONST
    One input: MONO, EXP, LOG, SIN, COS - integer from 1- NINPUTs, (-1) - all inputs
    multiple input: Multi2, Multi3, ratio
    Need <pow> for MONO, Multi2, Multi3, Ratio, (-1111) - all powers
    """

    print("Adding groups", groups)
    for g in groups:
        if len(g) == 1:
            addBasisGroup(g[0])
        elif len(g) == 2:
            addBasisGroup(g[0], g[1])
        elif len(g) == 3:
            addBasisGroup(g[0], g[1], g[2])


def addBasisGroup(type_of_function, input_indices="", powers=""):
    """
    Include the checks
    """

    global group_list
    group_index = len(group_list) + 1
    if type_of_function == "CONST":
        group_list.append("%i %s" % (group_index, type_of_function))
    elif type_of_function in ["LIN", "LOG", "EXP", "SIN", "COS"]:
        if input_indices == "":
            input_indices = "-1"
        group_list.append("%i %s %s" % (group_index, type_of_function, input_indices))
    elif type_of_function in ["MONO", "MULTI2", "MULTI3", "RATIO"]:
        if input_indices == "":
            input_indices = "-1"
        if powers == "":
            powers = "-1111"
        group_list.append(
            "%i %s %s %s" % (group_index, type_of_function, input_indices, powers)
        )
    elif type_of_function in ["GRP", "RBF", "CUST"]:
        group_list.append("%i %s %s" % (group_index, type_of_function, input_indices))


def clearBasisGroups():
    global group_list
    group_list = []


basis_constraint_list = []


def addBasisConstraints(groups_constraints):
    """
    NMT (no more than), ATL (at least), REQ (required if main group), XCL (exclude)
    format: group-id output_id constraint_type integer_parameter
    """

    # PYLINT-TODO-FIX "group_constraints" and "basis_constraints" are undefined
    # pylint: disable=undefined-variable
    global group_constraints, basis_constraints
    print("Adding basis selection constraint", group_constraints)

    for b in basis_constraints:
        addBasisConstraint(b[0], b[1], b[2], b[3])


def addBasisConstraint(group_id, output_id, constraint_type, intParam):
    global basis_constraint_list
    basis_constraint_list.append(
        "%i %i %s %i" % (group_id, output_id, constraint_type, intParam)
    )


def clearBasisConstraints():
    global basis_constraint_list
    basis_constraint_list = []


extrapxmin = None
extrapxmax = None
surface_constraint_list = []


def addCustomConstraints(constraint_list, **kwargs):
    """
    constraint args = CRTOL, CRNINITIAL, CRMAXITER, CRNVIOL, CRNTRIALS
    """

    print("Adding constraints", constraint_list)

    global extrapxmin, extrapxmax, surface_constraint_list
    if "extrapxmin" in kwargs.keys():
        extrapxmin = kwargs["extrapxmin"]
    if "extrapxmax" in kwargs.keys():
        extrapxmax = kwargs["extrapxmax"]
    for c in constraint_list:
        surface_constraint_list.append(c)


custom_fxn_list = []


def addCustomFunctions(fxn_list):
    print("Adding Functions", fxn_list)

    global custom_fxn_list

    for f in fxn_list:
        custom_fxn_list.append(f)
    print(custom_fxn_list)


def clearCustomFunctions():
    global custom_fxn_list
    custom_fxn_list = []


def writeCustomALAMOOptions(kwargs):
    global custom_fxn_list, surface_constraint_list, extrapxmax, extrapxmin, group_list, basis_constraint_list
    name = "almopt.txt"

    if kwargs.get("almopt", None) is not None:
        name = kwargs["almopt"]

    with open(name, "w") as r:
        if len(custom_fxn_list) > 0:
            r.write("\nNCUSTOMBAS %d\n" % len(custom_fxn_list))
            r.write("BEGIN_CUSTOMBAS\n")
            for c_fxn in custom_fxn_list:
                r.write("%s\n" % c_fxn)
            r.write("END_CUSTOMBAS\n\n")

        if len(surface_constraint_list) > 0:
            if extrapxmin is not None:
                r.write("extrapxmin %d\n" % extrapxmin)
            if extrapxmax is not None:
                r.write("extrapxmax %d\n" % extrapxmax)

        r.write("crncustom %i\n" % len(surface_constraint_list))
        i = 1
        r.write("BEGIN_CUSTOMCON\n")
        for c in surface_constraint_list:
            r.write("%s\n" % (c))
            i += 1
        r.write("END_CUSTOMCON\n")

        if len(group_list) > 0:
            r.write("\nNGROUPS %d\n" % len(group_list))
            r.write("BEGIN_GROUPS\n")
            for group in group_list:
                r.write("%s\n" % group)
            r.write("END_GROUPS\n\n")

        if len(basis_constraint_list) > 0:
            r.write("BEGIN_GROUPCON\n")
            for b in basis_constraint_list:
                r.write("%s\n" % b)
            r.write("END_GROUPCON\n\n")


# Argument management


def manageArguments(xdata, zdata, data, debug, kwargs):
    """
    Parse additional input options

    The 'pargs' library is used to keep track of options a user has availible
    descriptions of the dictionaries data, and debug are given in shared.py
    Multiple keys used to make writing the .alm file easier

    Args:
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real])
        data:  shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
    """

    parseKwargs(data, debug, kwargs)

    # Check to see if a simwrapper should be built
    if debug.get("simwrap", None) or kwargs.get("simulator", None) is not None:
        buildSimWrapper(data, debug)

    # Specific check to see if the labels of the response variables
    # should be used in the output dictionary
    # This is important for systematic testing vs. single model input
    if data["opts"]["noutputs"] > 1:
        # 'Must use outkeys for multiple outputs'
        debug["outkeys"] = True

    # Construct xmin and xmax vector based on training data if not provided
    if kwargs.get("xmin", None) is None:
        constructXBounds(xdata, zdata, data, debug)


def parseKwargs(data, debug, kwargs):
    """
    Parse keyword arguments

    Args:
        data: shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
        kwargs: keyword arguments
    """

    for arg in kwargs.keys():
        if kwargs[arg] is not None:
            if arg in data["pargs"]["opts"]:
                data["opts"][arg] = kwargs[arg]
            elif arg in data["pargs"]["lstopts"]:
                data["lstopts"][arg] = list()
                for term in list([kwargs[arg]]):
                    data["lstopts"][arg].append(term)
            elif arg in data["pargs"]["stropts"]:
                data["stropts"][arg] = kwargs[arg]
            elif arg in data["pargs"]["set4"]:
                data["set4"][arg] = kwargs[arg]
            elif arg in debug["pargs"]:
                debug[arg] = kwargs[arg]
            else:
                if arg not in (["xlabels", "zlabels", "xval", "zval"]):
                    sys.stdout.write("Problem with option : " + arg)
                    # PYLINT-TODO-FIX use almerror correctly insted of calling it directly
                    # pylint: disable=not-callable
                    almerror("p3")


def buildSimWrapper(data, debug):
    """
    Builds an executable simulator to sample for data

    Args:
        data: shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
    """

    if not isinstance(data["stropts"]["simulator"], type("string")):
        try:
            data["stropts"]["simulator"] = data["stropts"]["simulator"].__name__
        except Exception:
            raise almerror.AlamoInputError(
                "Simulator must be provided as a string"
                "and obey ALAMOs simulator conventions"
                "OR must be a python function whose name"
                "can be obtained via .__name__"
            )
        data["stropts"]["simulator"] = alamopy.wrapwriter(data["stropts"]["simulator"])
        debug["simwrap"] = True


def constructXBounds(xdata, zdata, data, debug):
    """
    Construct xmin,xmax and zmin, zmax for alamo if none are given

    Args:
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real])
        data: shared alamo data options
        debug: Additional options may be specified and will be applied
                to the .alm
    """

    xmin = ""
    xmax = ""
    if data["opts"]["ninputs"] > 1:
        for i in range(data["opts"]["ninputs"]):
            tn = debug["bignum"]
            tx = -1 * debug["bignum"]
            for j in range(data["opts"]["ndata"]):
                if float(xdata[j][i]) < float(tn):
                    tn = xdata[j][i]
                if float(xdata[j][i]) > float(tx):
                    tx = xdata[j][i]
            xmin = xmin + str(tn) + " "
            xmax = xmax + str(tx) + " "
    else:
        tn = debug["bignum"]
        tx = -1 * debug["bignum"]
        for j in range(data["opts"]["ndata"]):
            if float(xdata[j]) < float(tn):
                tn = xdata[j]
            if float(xdata[j]) > float(tx):
                tx = xdata[j]
        xmin = "%s %f " % (xmin, tn)
        xmax = xmax + str(tx) + " "
    data["set4"]["xmax"] = xmax
    data["set4"]["xmin"] = xmin


def checkForSampledData(data, debug):
    """
    Check to see if data has been sampled and update ndata

    Args:
        data: shared alamo data options
        debug: Additional options may be specified
               and will be applied to the .alm
    """

    lst_name = data["stropts"]["almname"].split(".")[0] + ".lst"
    with open(lst_name) as infile, open("awkres", "w") as outfile:
        copy = False
        for line in infile:
            if "Errors on observed data points" in line.strip():
                copy = True
            elif "Maximum absolute errors" in line.strip():
                copy = False
            elif copy:
                outfile.write(line)
    f = open("awkres")
    lf = f.read()
    f.close()
    lf2 = lf.split("\n")
    lf2 = lf2[1:-1]
    sys.stdout.write(
        "Updating number of training points from "
        + str(data["opts"]["ndata"])
        + " to "
        + str(len(lf2))
        + "\n"
    )
    data["opts"]["ndata"] = len(lf2)
    xdata = np.zeros([data["opts"]["ndata"], data["opts"]["ninputs"]])
    zdata = np.zeros([data["opts"]["ndata"], data["opts"]["noutputs"]])
    for i in range(len(lf2)):
        lf3 = lf2[i].split(" ")
        while "" in lf3:
            lf3.remove("")
        for j in range(data["opts"]["ninputs"]):
            xdata[i][j] = float(lf3[j])
        for j in range(data["opts"]["noutputs"]):
            zdata[i][j] = float(lf3[data["opts"]["ninputs"] + j])
    deletefile("awkres")
    return xdata, zdata


def expandOutput(xdata, zdata, vargs, data, debug):
    """
    Expand output to validation metrics and labels

    Args:
        data/debug: shared default options for .alm file
        xdata: (numpy.array or list[real])
        zdata: (numpy.array or list[real])
        vargs: Validation data
    """

    if debug["expandoutput"]:
        data["results"]["xdata"] = xdata
        data["results"]["zdata"] = zdata
        data["results"]["xlabels"] = data["labs"]["savexlabels"]
        data["results"]["zlabels"] = data["labs"]["savezlabels"]

    if debug["expandoutput"]:
        data["results"]["model"] = {}
        data["results"]["f(model)"] = {}
    else:
        for key in list(
            [
                "model",
                "f(model)",
                "ssr",
                "R2",
                "size",
                "rmse",
                "nbas",
                "totaltime",
                "olrtime",
                "miptime",
                "clrtime",
                "othertime",
                "version",
                "status",
                "madp",
                "numolr",
                "nummip",
                "numclr",
                "ninputs",
            ]
        ):
            data["results"][key] = collections.OrderedDict()
        if len(vargs) > 0:
            data["results"]["ssrval"] = collections.OrderedDict()
            data["results"]["R2val"] = collections.OrderedDict()
            data["results"]["rmseval"] = collections.OrderedDict()
            data["results"]["madpval"] = collections.OrderedDict()
    if debug["loo"]:
        data["results"]["Q2"] = data["results"]["Q2"]


# External File management


def readTraceFile(vargs, data, debug):
    """
    Read the alamo trace file to read in the model and metrics

    Args:
        data/debug: shared default options for .alm file
        vargs: Validation data
    """

    trace_file = data["stropts"]["tracefname"]
    trace_str = trace_file  # currentDirectory + "/" + trace_file
    try:
        lf = open(trace_str).read()
    except (IOError, FileNotFoundError) as err:
        if debug["mock"]:
            _construct_mock(data)
            return
        else:
            error_message = _diagnose_alamo_failure(trace_str, err)
            raise almerror.AlamoError(error_message)

    try:
        from sympy.parsing.sympy_parser import parse_expr
        from sympy import symbols, lambdify
    except Exception:
        # PYLINT-TODO-FIX "writethis" is not defined
        # pylint: disable=undefined-variable
        writethis("Cannot import sympy")

    lf2 = lf.split("\n")
    lf2_ind = 0  # Allows for multiple writings to trace.trc file

    dict_out_str = (
        "#filename, NINPUTS, NOUTPUTS, INITIALPOINTS, OUTPUT, SET, "
        "INITIALIZER, SAMPLER, MODELER, BUILDER, GREEDYBUILD, "
        "BACKSTEPPER, GREEDYBACK, REGULARIZER, SOLVEMIP, SSEOLR, SSE, "
        "RMSE, R2, ModelSize, BIC, RIC, Cp, AICc, HQC, MSE, SSEp, MADp, "
        "OLRTime, numOLRs, OLRoneCalls, OLRoneFails, OLRgsiCalls, OLRgsiFails, "
        "OLRdgelCalls, OLRdgelFails, OLRclrCalls, OLRclrFails, OLRgmsCalls, "
        "OLRgmsFails, CLRTime, numCLRs, MIPTime, NumMIPs, LassoTime, "
        "Metric1Lasso, Metric2Lasso, LassoSuccess, LassoRed, nBasInitAct, "
        "nBas, SimTime, SimData, TotData, NdataConv, OtherTime, NumIters, "
        "IterConv, TimeConv, Step0Time, Step1Time, Step2Time, TotalTime, "
        "AlamoStatus, AlamoVersion, Model"
    )

    lf2_ind = len(lf2) - 1 - lf2[::-1].index(dict_out_str)
    tkeys = lf2[lf2_ind].split(",")
    kl1 = list(
        [
            "ssr",
            "rmse",
            "R2",
            "size",
            "nbas",
            "totaltime",
            "olrtime",
            "miptime",
            "clrtime",
            "othertime",
            "version",
            "status",
            "madp",
            "numolr",
            "nummip",
            "numclr",
            "ninputs",
        ]
    )
    kl2 = list(
        [
            " SSE",
            " RMSE",
            " R2",
            " ModelSize",
            " nBasInitAct",
            " TotalTime",
            " OLRTime",
            " MIPTime",
            " CLRTime",
            " OtherTime",
            " AlamoVersion",
            " AlamoStatus",
            " MADp",
            " numOLRs",
            " NumMIPs",
            " numCLRs",
            " NINPUTS",
        ]
    )
    # Construct results for training data (&val if provided)
    ln = 1
    wlparam = data["opts"]["noutputs"]

    # initialize multiple output expanded dictionary
    if debug["expandoutput"]:
        data["results"]["f(model)"] = {}
        data["results"]["model"] = {}
        for i in range(len(kl1)):
            data["results"][kl1[i]] = {}
        if len(vargs) > 0:
            for i in ["ssrval", "R2val", "rmseval", "madpval"]:
                data["results"][i] = {}

    if len(vargs) > 0:
        wlparam = 2 * wlparam
    else:
        wlparam = wlparam + 1
    while ln < wlparam:
        lf3 = lf2[lf2_ind + ln].split(",")
        # Reapply the saved labels for the output
        model = lf3[tkeys.index(" Model")]
        #    for label in data['labs']['savexlabels']:
        for i in range(data["opts"]["ninputs"]):
            label = data["labs"]["xlinks"][i][0]
            # Now is a convenient time to collect information that will be used in the
            # confidence interval analysis
            model = model.replace(str(label), str(data["labs"]["xlinks"][i][1]))
        for i in range(data["opts"]["noutputs"]):
            label = data["labs"]["zlinks"][i][0]
            model = model.replace(str(label), str(data["labs"]["zlinks"][i][1]))

        # determine which output label to write
        # if debug['outkeys'] == True use olab as a key if not dont
        if debug["expandoutput"] or debug["outkeys"]:
            olab = model.split("=")[0]
            olab = olab.replace(" ", "")

            data["results"]["model"][olab] = model

            # Record tokenized model for each output
            data["results"]["f(model)"][olab] = lambdify(
                [symbols(data["labs"]["savexlabels"])],
                parse_expr(model.split("=")[1].replace("^", "**")),
                "numpy",
            )
        else:
            data["results"]["model"] = model
            data["results"]["f(model)"] = lambdify(
                [symbols(data["labs"]["savexlabels"])],
                parse_expr(model.split("=")[1].replace("^", "**")),
                "numpy",
            )

        if debug["expandoutput"]:
            if debug["outkeys"]:
                for i in range(len(kl1)):
                    data["results"][kl1[i]][olab] = lf3[tkeys.index(kl2[i])]
                # Check for validation set
                if len(vargs) > 0:
                    lf3 = lf2[lf2_ind + 2].split(",")
                    data["results"]["ssrval"][olab] = lf3[tkeys.index(" SSE")]
                    data["results"]["R2val"][olab] = lf3[tkeys.index(" R2")]
                    data["results"]["rmseval"][olab] = lf3[tkeys.index(" RMSE")]
                    data["results"]["madpval"][olab] = lf3[tkeys.index(" MADp")]
            else:
                for i in range(len(kl1)):
                    data["results"][kl1[i]] = lf3[tkeys.index(kl2[i])]
                # Check for validation set
                if len(vargs) > 0:
                    lf3 = lf2[lf2_ind + 2].split(",")
                    data["results"]["ssrval"] = lf3[tkeys.index(" SSE")]
                    data["results"]["R2val"] = lf3[tkeys.index(" R2")]
                    data["results"]["rmseval"] = lf3[tkeys.index(" RMSE")]
                    data["results"]["madpval"] = lf3[tkeys.index(" MADp")]
        else:
            if debug["outkeys"]:
                data["results"]["ssr"][olab] = lf3[tkeys.index(kl2[0])]
            else:
                data["results"]["ssr"] = lf3[tkeys.index(kl2[0])]
        ln = ln + 1


def cleanFiles(data, debug, pywrite=False, **kwargs):
    """
    Removes intermediate files

    Args:
        data/debug: shared default options for .alm file
        vargs: Validation data
    """

    # Delete files
    if debug["mock"]:
        try:
            deletefile("temp.alm")
            if pywrite:
                for z in alamopy.data["zlabels"]:
                    deletefile("%s.py" % z)
            return
        except Exception:
            pass

    if not debug["savepyfcn"]:
        for z in alamopy.data["results"]["zlabels"]:
            deletefile("%s.py" % z)

    if not debug["savescratch"]:
        deletefile(
            str(data["stropts"]["almname"])
            + " "
            + str(data["stropts"]["almname"]).split(".")[0]
            + ".lst"
        )
    if not debug["savetrace"]:
        deletefile(data["stropts"]["tracefname"])
    if not debug["saveopt"]:
        global custom_fxn_list, surface_constraint_list, extrapxmin, extrapxmax
        custom_fxn_list = []
        surface_constraint_list = []
        extrapxmax = None
        extrapxmin = None
        if kwargs.get("almopt", None) is not None:
            deletefile(kwargs["almopt"])

    if debug["simwrap"]:
        deletefile("simwrapper.py")


def _construct_mock(data):
    data["results"]["clrtime"] = "0"
    data["results"]["size"] = "6"
    data["results"]["numolr"] = "16960"
    data["results"]["othertime"] = "0.8799995E-01"
    data["results"]["olrtime"] = "0.10800002"
    data["results"]["miptime"] = "0"
    data["results"]["version"] = "2018.4.3"
    data["results"]["status"] = "0"
    data["results"]["R2"] = "1"
    data["results"]["numclr"] = "0"
    data["results"]["nummip"] = "0"
    data["results"]["ssr"] = "0.169E-21"
    data["results"]["pymodel"] = "cam6alm"
    data["results"]["totaltime"] = "0.1760001"
    data["results"]["rmse"] = "0.255E-11"
    data["results"]["madp"] = "0.814E-09"
    data["results"][
        "model"
    ] = "  z1 = 3.9999999999884194856747 * x1^2 \
         - 3.9999999999873385725380 * x2^2 - 2.0999999999876837186719 \
         * x1^4 + 3.9999999999879496392907 * x2^4 + 0.33333333333014281141260 \
         * x1^6 + 1.0000000000008837375276 * x1*x2"
    data["results"]["nbas"] = "15"

    # PYLINT-TODO-FIX "debug" is undefined
    # pylint: disable=undefined-variable
    if debug["expandoutput"]:
        data["results"]["ssrval"] = 0
        data["results"]["R2val"] = 0
        data["results"]["rmseval"] = 0
        data["results"]["madpval"] = 0


def _diagnose_alamo_failure(trace_output, error):
    error_message = None
    b_alamo = has_alamo()
    lf_logscratch = open("logscratch").read()
    if not b_alamo:
        error_message = "Alamo cannot be found. Please check Alamo is installed."
    elif "termination code" in lf_logscratch:
        error_message = "{}".format(lf_logscratch)

    else:
        error_message = 'Cannot read from trace file "{}": {}'.format(
            trace_output, error
        )
    return error_message


def get_alamo_version():
    x = [0, 1]
    z = [0, 1]

    res = alamo(x, z, saveopt=False)  # xval=xival, zval=zival, mock=True)

    try:
        deletefile("logscratch")
        deletefile("almopt.txt")
        deletefile("z1.py")
    except Exception:
        pass
    return res["version"]
