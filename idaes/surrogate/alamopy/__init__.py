###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
#
###############################################################################
__all__ = ['alamo', 'doalamo', 'getlabels', 'makelabs', 'addCustomFunctions', 
           'addCustomConstraints', 'addBasisGroups', 'addBasisGroup', 
           'addBasisConstraints', 'addBasisConstraint',
           'get_alamo_version', 'almwriter',
           'almconfidence', 'almpickle', 'postpickle', 'almunpickle',
           'AlamoError', 'AlamoInputError',
           'almpywriter', 'almcvwriter', 'wrapwriter',
           'almplot', 'mapminmax', 'remapminmax', 
           'data', 'debug', 'deletefile', 'movefile', 'catfile', 
           'copyfile', 'has_alamo', 'allcard', 'almlsq',
           'almlsqjac', 'almfeatmat',
           'ackley', 'branin', 'sixcamel', 'col']

from .doalamo import alamo, doalamo, getlabels, makelabs, addCustomFunctions, \
    addCustomConstraints, addBasisGroups, addBasisGroup, addBasisConstraints, \
    addBasisConstraint, get_alamo_version
from .allcard import allcard, almlsq, almlsqjac, almfeatmat
from .almwriter import almwriter
from .almconfidence import almconfidence
from .almpickle import almpickle, postpickle, almunpickle
from .almerror import AlamoError, AlamoInputError
from .almpywriter import almpywriter, almcvwriter, wrapwriter
from .mapminmax import mapminmax
from .remapminmax import remapminmax
from .shared import data, debug, initialize
from .almplot import almplot
from .multos import deletefile, movefile, catfile, copyfile, has_alamo
from .examples import sixcamel, ackley, branin

# Initializes all values to be used in the .alm file
import collections as col
