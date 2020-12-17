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
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_unucertainty
import pandas as pd
from NRTL_model_scripts import NRTL_model, NRTL_model_opt


variable_name = ["fs.properties.tau[benzene,toluene]", "fs.properties.tau[toluene,benzene]"]

#  https://github.com/IDAES/examples-pse/blob/main/src/Tutorials/Advanced/ParamEst/BT_NRTL_dataset.csv
data = pd.read_csv('BT_NRTL_dataset.csv')

def SSE(model, data):
    expr = ((float(data["vap_benzene"]) -
             model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"])**2 +
            (float(data["liq_benzene"]) -
             model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"])**2)
    return expr*1E4

obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(NRTL_model,NRTL_model_opt, data, variable_name, SSE)

print("\nEstiamted theta")
print(theta)
print("\nCovarinace of theta")
print(cov)
print("\nError propagation in the objective function, var(f(theta))")
print(propagation_f)
print("\nError propagation in the constraints, var(c(theta))")
print(propagation_c)
