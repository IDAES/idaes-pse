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
