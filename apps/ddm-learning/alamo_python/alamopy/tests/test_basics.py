##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Smoke tests, to make sure things are working at all.
"""
from alamopy import doalamo, almconfidence, almplot, wrapwriter
from alamopy.multos import deletefile
import numpy as np
import examples

def test_basic(**kwargs):

	x = [[ 0.17361977, -0.44326123], [-0.30192964,  0.68955226],[-1.98112458, -0.75686176],[0.68299634,  0.65170551],[-1.45317364,  0.15018666],[ 1.56528782, -0.58159576] ,[-1.25868712, -0.78324622],[-1.12121003,  0.95724757] ,[ 1.2467326,  -0.65611797],[ 1.26489899, -0.45185251]] 
	z = [-0.58978634828943055, -0.85834512885363479, 4.0241154669754113, 0.91057814668811488, 1.9147616212616931, 0.29103827202206878, 2.4290896722960778, 0.99199475534877579, 0.59688699266830847, 1.167850366995701]
	xival = [[5,5], [2,2]]
	zival = [5,2]

	doalamo(x, z) #, xval=xival, zval=zival, mock=True)    
	doalamo(x, z, mock=True, xlabels=["T", "V"], zlabels= ["P"])
	doalamo(x, z, xval=xival, zval=zival, mock=True, lmo=5)


	ndata=10
	x = np.random.uniform([-2,-1],[2,1],(ndata,2))
	z = [0]*ndata
	# specify simulator as examples.sixcamel
	sim = examples.sixcamel
	for i in range(ndata):
	    z[i]=sim(x[i][0],x[i][1])

	# Use alamopy's python function wrapper to avoid using ALAMO's I/O format
	almsim = wrapwriter(sim)
	# Call alamo through the alamopy wrapper
	res = doalamo(x,z,almname='cam6',monomialpower=(1,2,3,4,5,6), multi2power=(1,2), simulator=almsim, expandoutput=True,maxiter=20, mock=True)#,cvfun=True)
	conf_inv = almconfidence(res)

	res = doalamo(x, z, xval = xival, zval = zival) #, expandoutput=True) #, mock=True, loo=True) #, xval=xival, zval=zival, mock=True, loo=True)    

	conf_inv = almconfidence(res)
	print('Model: {}'.format(res['model']))
	print('Confidence Intervals : {}'.format(conf_inv['conf_inv']))
	almplot(res)

	try:
		deletefile("logscratch")
		deletfile("../cam6")
		deletefile("../cam6alm.py")
	except:
		pass

test_basic()