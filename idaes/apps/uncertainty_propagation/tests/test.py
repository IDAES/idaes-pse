import sys
import os
sys.path.append(os.path.abspath('..')) # current folder is ~/tests
import numpy as np
import pandas as pd
from scipy import sparse
import pytest
from pytest import approx
from mock import patch
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_uncertainty, propagate_uncertainty,clean_variable_name
from pyomo.opt import SolverFactory
from pyomo.environ import *
import pyomo.contrib.parmest.parmest as parmest
from sens import get_dsdp,get_dfds_dcds
ipopt_available = SolverFactory('ipopt').available()
kaug_available = SolverFactory('k_aug').available()
dotsens_available = SolverFactory('dot_sens').available()


variable_name = ['asymptote', 'rate_constant']
theta={'asymptote': 19.142575284617866, 'rate_constant': 0.53109137696521}
cov=np.array([[ 6.30579403, -0.4395341 ],[-0.4395341 ,  0.04193591]])
model_uncertain= ConcreteModel()
model_uncertain.asymptote = Var(initialize = 15)
model_uncertain.rate_constant = Var(initialize = 0.5)
model_uncertain.obj = Objective(expr = model_uncertain.asymptote*( 1 - exp(-model_uncertain.rate_constant*10  )  ), sense=minimize)
theta= {'asymptote': 19.142575284617866, 'rate_constant': 0.53109137696521}
for v in variable_name:
    getattr(model_uncertain, v).setlb(theta[v])
    getattr(model_uncertain, v).setub(theta[v])
gradient_f, gradient_c, col,row, line_dic=  get_dfds_dcds(model_uncertain, variable_name)
dsdp, col1 =  get_dsdp(model_uncertain, variable_name, theta, {})

print(gradient_f)
print(gradient_c)
print(col)
print(row)
print(line_dic)

print(dsdp.toarray())
print(col1)
