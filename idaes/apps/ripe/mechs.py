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

import pyomo.environ as pyo

_all__ = [
    "powerlawp5",
    "powerlaw2",
    "powerlaw3",
    "powerlaw4",
    "avrami2",
    "avrami3",
    "avrami4",
    "avrami5",
    "randomnuc",
    "ptompkins",
    "jander",
    "antijander",
    "valensi",
    "parabolic",
    "gb3d",
    "zlt",
    "grain",
    "massact",
    "massactm",
]

# file containing possible reaction mechanisms
# Due to contraints placed in other subroutines any mechanisms defined here can only take
# process data of an anticipated form [indep variables, conc data]

# Mechanisms defined in ripe have a standard input form due
# To restraints placed on the models by rbfopt in ems
# RIPE standard input form is :
# data = [x,T,flow,vol,t] where x is length ns
# data is not a list of list but a continuous list
# Additional options can be specified for user-specified mechanisms


def powerlawp5(*data):
    pd = data[0]
    return 0.5 * (pd ** (0.5 - (1.0 / 0.5)))


def powerlaw2(*data):
    pd = data[0]
    return 2.0 * (pd ** (2.0 - (1.0 / 2.0)))


def powerlaw3(*data):
    pd = data[0]
    return 3.0 * (pd ** (3.0 - (1.0 / 3.0)))


def powerlaw4(*data):
    pd = data[0]
    return 4.0 * (pd ** (4.0 - (1.0 / 4.0)))


def avrami2(*data):
    pd = data[0]
    return 2.0 * (1.0 - pd) * (-1 * pyo.log(1 - pd)) ** (2.0 - (1.0 / 2.0))


def avrami3(*data):
    pd = data[0]
    return 3.0 * (1.0 - pd) * (-1 * pyo.log(1 - pd)) ** (3.0 - (1.0 / 3.0))


def avrami4(*data):
    pd = data[0]
    return 4.0 * (1.0 - pd) * (-1 * pyo.log(1 - pd)) ** (4.0 - (1.0 / 4.0))


def avrami5(*data):
    pd = data[0]
    return 5.0 * (1.0 - pd) * (-1 * pyo.log(1 - pd)) ** (5.0 - (1.0 / 5.0))


def randomnuc(*data):
    pd = data[0]
    return 1.0 - pd


def ptompkins(*data):
    pd = data[0]
    return pd * (1.0 - pd)


def jander(*data):
    pd = data[0]
    return (
        3.0 * ((1.0 - pd) ** (1.0 / 3.0)) * (1.0 / (1.0 + pd) ** ((-1.0 / 3.0) - 1.0))
    )


def antijander(*data):
    pd = data[0]
    return ((3.0 / 2.0) * (1.0 - pd) ** (2.0 / 3.0)) * (
        1.0 / (1.0 + pd) ** ((-1.0 / 3.0) - 1.0)
    )


def valensi(*data):
    pd = data[0]
    return 1.0 / (-1.0 * pyo.log(1.0 - pd))


def parabolic(*data):
    pd = data[0]
    return 1.0 / (2.0 * pd)


def gb3d(*data):
    pd = data[0]
    return (
        (3.0 / 2.0)
        * (1.0 - pd) ** (4.0 / 3.0)
        * 1.0
        / ((1.0 - pd) ** ((-1.0 / 3.0) - 1.0))
    )


def zlt(*data):
    pd = data[0]
    # zhuralev lesohkin tempelman model
    return (3.0 / 2.0) / ((1.0 - pd) ** ((-1.0 / 3.0) - 1.0))


def grain(*data):
    pd = data[0]
    # grain model
    return (1.0 - pd) ** (2.0 / 3.0)


# Mechanisms can be specified such that they require a mechanism
# EMS requires these models do not take the mechanism as an argument


def massactm(data, def_stoich):
    import numpy as np

    ns = len(def_stoich)
    pwr = []
    x = []
    for i in range(ns):
        if def_stoich[i] < 0:
            pwr.append(abs(def_stoich[i]))
            x.append(data[i])
    return np.product(np.power(x, pwr))


def mechperstoich(mech, stoich):
    import numpy as np

    ns = len(stoich)
    # out_mechs = []

    def massact(*data):
        pwr = []
        x = []
        for i in range(ns):
            if stoich[i] < 0:
                pwr.append(abs(stoich[i]))
                x.append(data[i])
        return np.product(np.power(x, pwr))

    def usr_f(*data):
        #        stoich = stoich
        return mech(data, stoich)

    if mech == "massact":
        f = massact
    else:
        f = usr_f

    return f


# This subroutine analyzes inputs to construct kinetic mechanisms
def getmechs(kwargs):
    import numpy as np

    # Analyze mechanisms in call to ripe
    # determine which stoichiometries to apply each mechanism to
    inkeys = kwargs.keys()
    # Analyze stoichiometries for consistency and relevant information
    # Ideally this input is formatted as a list of list, but can be int, float, or list
    if "stoich" not in inkeys:
        stoich = [[1]]
    else:
        stoich = kwargs["stoich"]
        if isinstance(stoich, type(0)) or isinstance(stoich, type(0.0)):
            stoich = [kwargs["stoich"]]
        elif len(np.shape(stoich)) == 1:
            stoich = [stoich]
        else:
            stoich = kwargs["stoich"]
    dshape = np.shape(stoich)
    if dshape == 1:
        stoich = [stoich]
        dshape = np.shape(stoich)
    nstoich = dshape[0]
    # if len(dshape) > 1:
    #     nspec = dshape[1]
    # else:
    #     nspec = 1
    # Initialize list for later use
    sfixlist = []
    mechlist = []
    # Check for specified mechanisms
    # massaction kinetics internal to ripe are used by default
    # Ideally mechanisms is a list of list of the following format
    # [[[list-of-stoich-indexes],[mechanisms]],...[[next],[sets]]
    if "mech" not in kwargs.keys():
        mechanisms = [[range(nstoich), ["massact"]]]
    else:
        dshape = np.shape(kwargs["mechanisms"])
        if len(dshape) == 1:  # and len(kwargs['mechanisms']) == 1:
            mechanisms = [[range(nstoich), kwargs["mechanisms"]]]
        elif len(dshape) == 0:
            mechanisms = [[range(nstoich), [kwargs["mechanisms"]]]]
        else:
            mechanisms = kwargs["mechanisms"]
    # Determine number of mechanisms specified
    ncons = 0
    for i in range(len(mechanisms)):
        # replace token 'all' with range(stoich)
        # Additionally, either element in mechanisms (may) not be list erroniously
        if not isinstance(mechanisms[i][0], (type([]), range)):
            mechanisms[i][0] = [mechanisms[i][0]]
        if mechanisms[i][0][0] == "all":
            mechanisms[i][0] = range(nstoich)
        # Ensure the mechanism(s) specified are in a list
        if not isinstance(mechanisms[i][1], type([])):
            mechanisms[i][1] = [mechanisms[i][1]]
        for m in mechanisms[i][1]:
            mechlist.append(m)
            sfixlist.append(mechanisms[i][0])
            ncons += len(mechanisms[i][0])
        if len(mechanisms[i]) < 3:
            mechanisms[i].append(False)
    # nm is the number of mechanisms across stoichiometries
    nm = len(mechlist)
    # fix array is used in the pyomo model to fix binaries for reaction not associated with stoichs
    fixarray = np.zeros([nstoich, nm])
    for j in range(nm):
        for st in sfixlist[j]:
            fixarray[st, j] = 1
    # Now that fixarray is constructed simplify mechanisms
    return stoich, mechanisms, fixarray, ncons, mechlist
