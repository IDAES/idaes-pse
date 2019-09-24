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
# This file contains a reactor simulator built in pyomo
# The simulator emulates behavior observed in the production of
# Styrene from ethylbenzene

import numpy as np
import pyomo.environ as pyo

# define fractional variance of noise, SNR = 1 / noise
noise = 0.05

# Reference temperature defined
Tr = 750.0
# Kinetic parameters are hard coded in the form [[k1,k2...],[E1,E2,...]]
kinetic_params = [[250,220,38,25],[115,131,55,75]]


def sim(data):
    import numpy as np
    # Enable 1/2d calls
    # Ensure that data sizes and shapes are consistent
    try:
        x0 = data[:,:9]
        Temp = data[:,9]
    except:
        x0 = data[:9]
        Temp = data[9]
    params = kinetic_params
    dshape = np.shape(x0)
    if len(dshape) == 1:
        x0 = np.expand_dims(x0, axis=-1)
        x0 = np.ndarray.transpose(x0)
        dshape = np.shape(x0)
    npc = dshape[0]
    ns = dshape[1]
    # Define the reactor simulation in pyomo
    def pyomosim(data):
        # Define rate parameters 
        kco = 35
        kst = 1.5
        sparam = 1.0
        flow = 1.0
        vol = 1.0
        gc = .008314
        # Define kinetic rate parameters
        k = params[0]
        E = params[1]
        # define UB for concentration
        #nound_ub = 20
        # pyomo solver options
        opt = pyo.SolverFactory('baron')
        cracmodel = pyo.ConcreteModel()
        ca0,cb0,cc0,cd0,cf0,cg0,ch0,ci0,cj0,T = [float(v) for v in data]
        
        bound_ub = 100.0
        # Define cracmodel variables
        # A = Eb , B = St , C = Bz , D = Et, E = Tl, F = Me, G = Water, H = H2 I = CO2, J = N2
        cracmodel.a = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = ca0)
        cracmodel.b = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = cb0)
        cracmodel.c = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = cc0)
        cracmodel.d = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = cd0)
        cracmodel.f = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = cf0)
        cracmodel.g = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = cg0)
        cracmodel.h = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = ch0)
        cracmodel.i = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = ci0)
        cracmodel.j = pyo.Var(domain = pyo.NonNegativeReals, bounds = (0,bound_ub), initialize = cj0)
        cracmodel.r1 = pyo.Var(domain = pyo.Reals)
        cracmodel.r2 = pyo.Var(domain = pyo.Reals)
        cracmodel.r4 = pyo.Var(domain = pyo.Reals)
        cracmodel.r5 = pyo.Var(domain = pyo.Reals)

        def keq(T):
            return pyo.exp(0.1 + (300 / T))
#            return pyo.exp(-1.0*(122700-126.3*T-0.002194*T**2)/(gc*T)) * sparam
        # define reaction rate variable
        
        def fr1(cracmodel):
            # A <> B + H
            return cracmodel.r1 == k[0] * pyo.exp(-(E[0]/(gc))*((1/T)-(1/Tr))) * (cracmodel.a - (cracmodel.b * cracmodel.h)/keq(T)) * (1.0/((1+kst*cracmodel.b)*(1+kco*cracmodel.i)))
        def fr2(cracmodel):
            # A > C + D
            return cracmodel.r2 == k[1] * pyo.exp(-(E[1]/(gc))*((1/T)-(1/Tr)))  * cracmodel.a * (1.0/(1+kco*cracmodel.i))

        def fr4(cracmodel):
            # D + 2H > 2F
            return cracmodel.r4 == k[2] * pyo.exp(-(E[2]/(gc))*((1/T)-(1/Tr)))  * cracmodel.d * cracmodel.h
        def fr5(cracmodel):
            # F+G > I + 4H
            return cracmodel.r5 == k[3] * pyo.exp(-(E[3]/(gc))*((1/T)-(1/Tr)))  * cracmodel.f * cracmodel.g

        cracmodel.er1 = pyo.Constraint( rule = fr1)
        cracmodel.er2 = pyo.Constraint( rule = fr2)
        cracmodel.er4 = pyo.Constraint( rule = fr4)
        cracmodel.er5 = pyo.Constraint( rule = fr5)
        cracmodel.sets = pyo.RangeSet(9)
        cracmodel.dum = pyo.Var(cracmodel.sets, domain = pyo.Reals)

        def fra(cracmodel):  # A - 1,2,3
            return cracmodel.dum[1]  == ca0-cracmodel.a - cracmodel.r1 -cracmodel.r2 
        def frb(cracmodel): # B - 1
            return cracmodel.dum[2]  == cb0-cracmodel.b+cracmodel.r1
        def frc(cracmodel): # C - 2
            return cracmodel.dum[3]  == cc0-cracmodel.c+cracmodel.r2
        def frd(cracmodel): # D - 2,4
            return cracmodel.dum[4]  == cd0-cracmodel.d+cracmodel.r2-cracmodel.r4
        def frf(cracmodel): # F - 3,4,5
            return cracmodel.dum[5]  == cf0-cracmodel.f+2*cracmodel.r4-cracmodel.r5
        def frg(cracmodel): # G - 5
            return cracmodel.dum[6]  == cg0-cracmodel.g-2*cracmodel.r5
        def frh(cracmodel): # H - 1,3,4,5
            return cracmodel.dum[7]  == ch0-cracmodel.h+cracmodel.r1-2*cracmodel.r4+4*cracmodel.r5
        def fri(cracmodel): # I - 5
            return cracmodel.dum[8]  == ci0-cracmodel.i + cracmodel.r5
        def frj(cracmodel): # J - N2 is inert
            return cracmodel.dum[9] == cj0-cracmodel.j

        cracmodel.era = pyo.Constraint( rule = fra)
        cracmodel.erb = pyo.Constraint( rule = frb)
        cracmodel.erc = pyo.Constraint( rule = frc)
        cracmodel.erd = pyo.Constraint( rule = frd)
        cracmodel.erf = pyo.Constraint( rule = frf)
        cracmodel.erg = pyo.Constraint( rule = frg)
        cracmodel.erh = pyo.Constraint( rule = frh)
        cracmodel.eri = pyo.Constraint( rule = fri)
        cracmodel.erj = pyo.Constraint( rule = frj)

        # minimize square of dummy variables to find steady-state concentrations
        def objf(cracmodel):
            return sum(cracmodel.dum[s]**2 for s in cracmodel.sets)
        
        cracmodel.OBJ = pyo.Objective(rule = objf)
        results = opt.solve(cracmodel)
        cracmodel.solutions.store_to(results)
        klist = ['a','b','c','d','f','g','h','i','j']
        # Add noise of the specifiec SNR, noise has variance eps ~ N(0,noise*conc)
        vn = [results.Solution.Variable[key]['Value']+np.random.normal(0,noise*results.Solution.Variable[key]['Value']) for key in klist]
        return vn
    
    # Simulate over requested datapoints
    # requested data may have 1 or more points
    concentrations = []
    if npc != 1:
        for i in range(npc):
            try:
                t2 = Temp[0][i]
            except:
                t2 = Temp[i]            
            conres = pyomosim(np.ndarray.tolist(x0[i,:])+[t2])

            concentrations.append(conres)
    else:
        conres = pyomosim(data)
        concentrations = conres
    return concentrations
