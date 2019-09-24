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
# This file contains an example application of the RIPE software
# The goal is to identify the reaction kinetics present in a reactor
# where ethylbenzene is converted to styrene, with a number of side products
# 
# Initial training sets provided through random or space-filling sampling 
# result in incorrect model identification
# Error maximization sampling can be used to refine the model

import pyomo.environ as pyo
from idaes.surrogate import ripe
import numpy as np
import random
from . import cracsim

np.random.seed(100)


def main():
    Tr = 750.0
    # Experimental variance is known in this problem,
    # it can be estimated if not provided
    noise = 0.05
    # ndata = 10 in publication example
    ndata = 30
    ns = 9
    # Define temperature bounds
    Tlo = 500
    Tup = 1000
    # Define range of inlet concentrations
    lb = [0,0,0,0,0,1,0,0,1]
    ub = [3]*ns

    gc = .008314

    Temp = np.linspace(Tlo,Tup,ndata)
    # Initialize concentration data

    # Inlet concentrations are fixed in publication example
    # cdata0 = [[.5,0,0,0,0,4.5,0,0,4.5]]*ndata
    cdata0=np.zeros([ndata,ns])
    for i in range(ndata):
        for j in range(ns):
            cdata0[i,j] = random.uniform(lb[j],ub[j])


    # Calculate steady-state concentration values from simulator cracsim.py
    cdata = cracsim.sim(np.hstack((cdata0,np.expand_dims(Temp, axis=1))))

    # In this example, we know the true stoichiometries. Lets define them first for clarity
    t_stoich = [[-1,1,0,0,0,0,1,0,0],[-1,0,1,1,0,0,0,0,0],[0,0,0,-1,2,0,-2,0,0],[0,0,0,0,-1,-2,4,1,0]]
    # Additional considered stoichiometries are defined
    a_stoich = [[-1,0,0,0,0,-16,21,8,0],[-1,0,0,4,0,0,-3,0,0],[0,0,0,-1,0,-4,6,2,0],[-1,0,1,0,2,0,-2,0,0]]
    # Index 0-3 are the true reactions, 4-7 are considered reactions
    stoichs = t_stoich+a_stoich

    # Define kinetic mechanisms, adsorption parameters must be known a-priori
    kco = 35
    kst = 1.5

    # Mechanisms can be defined for each stoichiometry using a list-of-list
    mechs = [[[0,1,3,4,7],eb_dep],[[0],[t_st_prod,cat_st_prod_t1]],[[0,1],cat_ben_prod_t2],[[2,3],meth_prod_t3],[[2,3],ch4_to_co_t4],[[2,3,4,5,6,7],[ma_g,ma_h]]]

    # Experimental variance is known in this case
    sigma = np.multiply(noise**2,cdata)

    results = ripe.ripemodel(cdata,stoich = stoichs,mechanisms=mechs,x0=cdata0,temp=Temp,sigma=sigma,tref=Tr)


def keq(*x):
    a,b,c,d,f,g,h,i,j,T = x
    temp = 0.1 + 300/T
    return pyo.exp(temp)

# These mechanisms are present in the simulation
def cat_st_prod_t1(*x):
    a,b,c,d,f,g,h,i,j,Temp = x
    return (a - (b * h) / keq(*x)) * (1/((1+kst*b)*(1+kco*i)))
def cat_ben_prod_t2(*x):
    # Mechanism for EB > B + C2H4
    a,b,c,d,f,g,h,i,j,Temp = x
    return a / (1+kco*i)
def meth_prod_t3(*x):
    # Mechanism for C2H4+4H2O > 2CO2+6H
    a,b,c,d,f,g,h,i,j,Temp = x
    return d * h
def ch4_to_co_t4(*x):
     a,b,c,d,f,g,h,i,j,Temp = x
    #mechanism for CH4+2H2O > CO2+4H2
     return f * g

# Additional mechanisms are specified for the true stoichiometries and additional stoichs
def eb_dep(*x):
    a,b,c,d,f,g,h,i,j,Temp = x
    return a
def t_st_prod(*x):
    a,b,c,d,f,g,h,i,j,Temp = x
    return (a - (b * h) / keq(*x))
def ma_h(*x):
    a,b,c,d,f,g,h,i,j,Temp = x
    return h
def ma_g(*x):
    a,b,c,d,f,g,h,i,j,Temp = x
    return g


# The following section of code can be un-commented in order to continue sampling data until an accurate model is obtained
# note that this is considerably more computationally expensive
'''
#Append T bounds for ems variables
lb.append(Tlo)
ub.append(Tup)

# Call RIPE ems in order to identify the next best sample point
[new_points, err] = ripe.ems(results,cracsim.sim,lb,ub,ns,x=cdata,x0=cdata0,Temp=Temp,Tref=Tr)#,frac=fracfun)
new_res = cracsim.sim(new_points)
ite = 0
data = cdata
data0 = cdata0
tdata = Temp.tolist()
# print 'maximum allowable tolerances : ', [2*noise*s for s in new_res]
while any(err >  [2*noise*s for s in new_res] ):
#    print 'Which concentration : ', err > [noise*s for s in new_res]
    data = np.vstack((data,new_res))
    data0 = np.vstack((data0,new_points[:-1]))
    tdata.append(new_points[-1])
    results = {}
    ite+=1
    sigma =  np.multiply(noise**2,np.array(data))
    results = ripe.ripemodel(data,stoich = stoichs,mechanisms=mechs,x0=data0,sigma=sigma,temp=tdata,tref=Tr, hide_output=True)
    [new_points, err] = ripe.ems(results,cracsim.sim,lb,ub,10,x=data,x0=data0,temp=tdata,tref=Tr)
    new_res = cracsim.sim(new_points)
#    print 'currently at '+str(len(data))+' data points'
#    print 'proposed new conc : ', new_res
#    print 'maximum allowable tolerances : ', [noise*s for s in new_res]
'''


if __name__ == "__main__":
    main()