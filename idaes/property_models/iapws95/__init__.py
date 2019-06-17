import os

def iapws95_available():
    plib = os.path.join(os.path.dirname(__file__), "iapws95.so")
    return os.path.isfile(plib)

from idaes.property_models.iapws.iapws95 import *
