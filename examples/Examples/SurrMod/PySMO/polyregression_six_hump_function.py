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
from idaes.surrogate.pysmo.polynomial_regression import *
from idaes.surrogate.pysmo import sampling as sp
import pandas as pd
import numpy as np
import os


def main():
    # Load XY data from high fidelity model from tab file using Pandas. Y-data must be in the last column.
    current_path = os.path.dirname(os.path.realpath(__file__))
    data = pd.read_csv(os.path.join(current_path, 'data_files', 'six_hump_function_data.tab'), sep='\t', header=0, index_col=0)

    b = sp.LatinHypercubeSampling(data, 30, 'selection')
    tr_data = b.sample_points()

    # Carry out polynomial regression
    d = PolynomialRegression(tr_data, tr_data, maximum_polynomial_order=8, multinomials=1, overwrite=True)
    p = d.get_feature_vector()
    d.training()

    # Print pyomo expression
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1,2])
    print("")
    print(d.generate_expression([m.x[1], m.x[2]]))

if __name__ == "__main__":
    main()
