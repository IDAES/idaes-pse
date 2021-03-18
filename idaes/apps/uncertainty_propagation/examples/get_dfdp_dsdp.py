##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import sys
import os
sys.path.append(os.path.abspath('..')) # current folder is ~/examples
from idaes.apps.uncertainty_propagation.uncertainties import get_sensitivity, get_dsdp, clean_variable_name
import pandas as pd
from pyomo.environ import *
import pyomo.contrib.parmest.parmest as parmest
'''
min f:  p1*x1+ p2*(x2^2) + p1*p2 
 s.t  c1: x1 + x2 = p1   
      c2: x2 + x3 = p2
          0 <= x1, x2, x3 <= 10
          p1 = 10
          p2 = 5 
'''

variable_name = ['p1', 'p2']

m= ConcreteModel()
m.x1 = Var(initialize = 0)
m.x2 = Var(initialize = 0)
m.x3 = Var(initialize = 0)
m.x1.setlb(0)
m.x1.setub(10)
m.x2.setlb(0)
m.x2.setub(10)
m.x3.setlb(0)
m.x3.setub(10)

m.p1 = Var(initialize = 0)
m.p2 = Var(initialize = 0)

m.obj = Objective(expr = m.p1*m.x1+ m.p2*m.x2*m.x2 + m.p1*m.p2 , sense=minimize)
m.c1 = Constraint(expr = m.x2 + m.x3 == m.p2)
m.c2 = Constraint(expr = m.x1 + m.x2 == m.p1)
theta= {'p1': 10.0, 'p2': 5.0}
for v in variable_name:
    getattr(m, v).setlb(theta[v])
    getattr(m, v).setub(theta[v])
gradient_f,gradient_f_dic, gradient_c,gradient_c_dic, line_dic =  get_sensitivity(m, variable_name)
dsdp_dic = get_dsdp(m, variable_name, theta)


