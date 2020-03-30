
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
# This file applies the ALAMOpython module to the six hump camel problem
# More information on this function can be found at :
#            https://www.sfu.ca/~ssurjano/camel6.html
# This problem utilizes ALAMO's sampling features

from idaes.surrogate.main import Pysmo_rbf, Pysmo_kriging, Pysmo_polyregression, Alamopy
from pyomo.environ import Var, ConcreteModel, Objective

import math
import examples
import numpy as np

def main():
    # if not alamopy.has_alamo():
    #     return False
    _main()
    return True

def _main():
    def branin_function(x1, x2):
        pi = 3.1417
        t1 = (x2 - (5.1 * x1 * x1 / (4 * pi * pi)) + (5 * x1 / pi) - 6) ** 2
        t2 = (10 * (1 - 1/(8 * pi))*np.cos(x1))
        y = t1 + t2 + 10
        return y

    # Create 25 random samples for training, 100 for testing
    np.random.seed(100)
    ndata = 25
    nval = 100
    x = np.random.uniform([-5, 0], [10, 15], (ndata, 2))
    y = branin_function(x[:, 0], x[:, 1])
    # xval = np.random.uniform([-5, 0], [10, 15], (nval, 2))
    # yval = branin_function(xval[:, 0], xval[:, 1])

    m = ConcreteModel()
    m.x = Var([1, 2])

    pysmo_rbf_settings = {'basis_function': 'gaussian',
                          'regularization': True,
                          'pyomo_vars': [m.x[1], m.x[2]],
                          'overwrite': True}
    modeler = Pysmo_rbf(**pysmo_rbf_settings)

    pysmo_krg_settings = {'numerical_gradients': True,
                          'regularization': True,
                          'pyomo_vars': [m.x[1], m.x[2]],
                          'overwrite': True}
    modeler = Pysmo_kriging(**pysmo_krg_settings)

    pysmo_pr_settings = {'maximum_polynomial_order':4,
                         'multinomials':1,
                         'pyomo_vars': [m.x[1], m.x[2]],
                         'training_split':0.9,
                         'number_of_crossvalidations': 5,
                         'overwrite': True}
    modeler = Pysmo_polyregression(**pysmo_pr_settings)

    alamo_settings = {'monomialpower':(1, 2, 3, 4, 5, 6),
                      'multi2power':(1, 2),
                      'expandoutput':True}
    modeler = Alamopy(**alamo_settings)


    modeler.regressed_data(x,y)
    print(modeler.get_regressed_data())

    modeler.build_model()

    print(modeler.get_results())
    m.obj = Objective(expr=modeler._model)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    modeler2 = Pysmo_polyregression()
    modeler2.load_results('results.pickle')
    print(modeler2._model)



if __name__ == "__main__":
    np.random.seed(100)
    main()