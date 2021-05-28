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
import numpy as np
noise = 0.1


def sim(data):
    x0 = data
    params = [1.5,2.1,0.9]
    #pyomo simulator for cracking example
    import numpy as np
    import pyomo.environ as pyo

    dshape = np.shape(x0)
    
    if len(dshape) == 1:
        x0 = np.expand_dims(x0, axis=-1)
        x0 = np.ndarray.transpose(x0)
        dshape = np.shape(x0)
    npc = dshape[0]
    ns = dshape[1]
    if dshape[0] == 1:
        xo = [x0]*npc
    
    def pyomosim(x0):
        # Define rate parameters 
        flow = 1.0
        vol = 1.0
        gc = 8.314
        # Define kinetic rate parameters
        k = params        
        a0,b0,c0,d0,e0 = np.array(x0).T.tolist()
        opt = pyo.SolverFactory('baron')
        model = pyo.ConcreteModel()
#        bound_ub = 100
        # Define model variables
        model.a = pyo.Var(domain = pyo.NonNegativeReals)#, bounds = (0.0,bound_ub))
        model.b = pyo.Var(domain = pyo.NonNegativeReals)#, bounds = (0.0,bound_ub))
        model.c = pyo.Var(domain = pyo.NonNegativeReals)#, bounds = (0.0,bound_ub))
        model.d = pyo.Var(domain = pyo.NonNegativeReals)#, bounds = (0.0,bound_ub))
        model.e = pyo.Var(domain = pyo.NonNegativeReals)#, bounds = (0.0,bound_ub))
        model.dset = pyo.RangeSet(5)
        model.dum = pyo.Var(model.dset)

        model.r1 = pyo.Var(domain = pyo.Reals)
        model.r2 = pyo.Var(domain = pyo.Reals)
        model.r3 = pyo.Var(domain = pyo.Reals)
        
        def fr1(model):
            return model.r1 == k[0] * model.a * model.b
        def fr2(model):
            return model.r2 == k[1] * model.b * model.c
        def fr3(model):
            return model.r3 == k[2] * model.a * model.d

        model.er1 = pyo.Constraint( rule = fr1)
        model.er2 = pyo.Constraint( rule = fr2)
        model.er3 = pyo.Constraint( rule = fr3)
        
        num = 1.0
        def fra(model):
            return num * model.dum[1] == (flow/vol)*(a0-model.a ) - model.r1 -model.r3
        def frb(model):
            return num * model.dum[2] == (flow/vol)*(b0-model.b ) - model.r1 - model.r2
        def frc(model):
            return num * model.dum[3] == (flow/vol)*(c0-model.c) + model.r1 - model.r2
        def frd(model):
            return num * model.dum[4] == (flow/vol)*(d0-model.d) + model.r2 - model.r3
        def fre(model):
            return num * model.dum[5] == (flow/vol)*(e0-model.e) + model.r3

        model.era = pyo.Constraint(rule=fra)
        model.erb = pyo.Constraint( rule = frb)
        model.erc = pyo.Constraint( rule = frc)
        model.erd = pyo.Constraint( rule = frd)
        model.ere = pyo.Constraint( rule = fre)

        def objf(model):
            return sum([ model.dum[i]**2 for i in model.dset])
        
        model.OBJ = pyo.Objective(rule = objf)

        results = opt.solve(model, tee=False)
        model.solutions.store_to(results)
        # Note: noise is not truly normally distributed as concentration values cannot be negative
#        print 'debug this : ', results.Solution.Variable['a']['Value'],np.random.normal(0,noise*results.Solution.Variable['a']['Value']),np.argmax([0.0,results.Solution.Variable['a']['Value']+np.random.normal(0,noise*results.Solution.Variable['a']['Value'])])
        v = [results.Solution.Variable[key]['Value'] for key in ['a','b','c','d','e']]
        vn = [results.Solution.Variable[key]['Value']+np.random.normal(0,noise*results.Solution.Variable[key]['Value']) for key in ['a','b','c','d','e']]
        tsum = 0
        for i in range(5):
            if vn[i] < 0:
                vn[i] = v[i]
                tsum+= 1
#        print 'total number of zeros : ', tsum
        return vn
    
    # Simulate over requested datapoints
    concentrations = []
    for i in range(npc):
        conres = pyomosim(x0[i])
        concentrations.append(conres)
    return concentrations
