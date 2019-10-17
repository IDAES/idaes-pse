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
import pyomo.environ as pyo
from idaes.surrogate import ripe
import numpy as np
import random
from . import isotsim

np.random.seed(20)


def main():

    #ndata = 100
    noise = 0.1
    ns = 5
    lb_conc = [0,0,0,0,0]
    ub_conc = [10,10,0,0,0]


    # Initialize concentration arrays

    # initial concentrations - only 2 data points at bounds
    cdata0 = [[1,1,0,0,0],[10,10,0,0,0]]
    cdata = isotsim.sim(cdata0)
    nd = len(cdata0)

    # Considered reaction stoichiometries
    stoich = [[-1,-1,1,0,0] ,[0,-1,-1,1,0],[-1,0,0,-1,1],[-1,-2,0,1,0] ,[-2,-2,0,0,1],[-1,-1,-1,0,1],[-2,-1,1,-1,1]]

    # IRIPE internal mass action kinetics are specified
    mechs = [['all','massact']]

    # Use expected variance - estimated from data if not provided
    sigma = np.multiply(noise**2,np.array(cdata))

    # Call to RIPE
    results = ripe.ripemodel(cdata,stoich = stoich,mechanisms=mechs,x0=cdata0,hide_output=False,sigma=sigma,deltaterm=0,expand_output=True)

    # Adaptive experimental design using error maximization sampling
    [new_points, err] = ripe.ems(results,isotsim.sim,lb_conc,ub_conc,5,x=cdata,x0=cdata0)

    # Implement EMS as described in the RIPE publication
    new_res = isotsim.sim(new_points)[0]
    ite = 0
    # print 'maximum allowable tolerances : ', [noise*s for s in new_res]
    while any(err >  [2*noise*s for s in new_res] ):
    #    print 'Which concentrations violate error (True=violation) : ', err > [noise*s for s in new_res]
        results = {}
        ite+=1
        # Data updated explicitly so RBFopt subroutines produce consistent results
        new_cdata0 = np.zeros([nd+ite,ns])
        new_cdata  = np.zeros([nd+ite,ns])
        new_cdata0[:-1][:] = cdata0[:][:]
        new_cdata[:-1][:] = cdata[:][:]
        new_cdata0[-1][:] = new_points[:]
        res = isotsim.sim(new_points)[0]
        for j in range(len(res)):
            new_cdata[-1][j] = res[j]

        #Update weight parameters
        sigma =  np.multiply(noise**2,np.array(new_cdata))

        # Build updated RIPE model
        results = ripe.ripemodel(new_cdata,stoich = stoich,mechanisms=mechs,x0=new_cdata0,sigma=sigma,expand_output=True)

        # Another call to EMS
        [new_points, err] = ripe.ems(results,isotsim.sim,lb_conc,ub_conc,5,x=cdata,x0=cdata0)

        # Update results
        new_res = isotsim.sim(new_points)[0]
        cdata0 = new_cdata0
        cdata = new_cdata

    # Final call to RIPE to get concise output
    results = ripe.ripemodel(cdata,stoich = stoich,mechanisms=mechs,x0=cdata0,sigma=sigma,expand_output=False)
    #print results


if __name__ == "__main__":
    main()



