__all__ = ['doalamo', 'writethis', 'almwriter', 'almconfidence', 'almpickle',
           'postpickle', 'almunpickle', 'almerror', 'almplot', 'almpywriter', 'allcard',
           'mapminmax', 'remapminmax', 'almcvwriter','wrapwriter']

from .doalamo import doalamo, getlabels, makelabs
from .allcard import allcard, almlsq, almlsqjac, almfeatmat
from .writethis import writethis
from .almwriter import almwriter
from .almconfidence import almconfidence
from .almpickle import almpickle, postpickle, almunpickle
from .almerror import almerror, AlamoError, AlamoInputError
from .almpywriter import almpywriter, almcvwriter, wrapwriter
from .mapminmax import mapminmax
from .remapminmax import remapminmax
from .shared import data, debug
from .almplot import almplot
from .multos import deletefile, movefile, catfile, copyfile

# Initializes all values to be used in the .alm file
import collections as col

