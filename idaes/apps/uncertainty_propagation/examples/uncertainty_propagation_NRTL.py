#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
import sys
import os

sys.path.append(os.path.abspath(".."))  # current folder is ~/examples
from idaes.apps.uncertainty_propagation.uncertainties import (
    quantify_propagate_uncertainty,
)
import pandas as pd
from idaes.apps.uncertainty_propagation.examples.NRTL_model_scripts import (
    NRTL_model,
    NRTL_model_opt,
)

if __name__ == "__main__":
    variable_name = [
        "fs.properties.tau['benzene', 'toluene']",
        "fs.properties.tau['toluene','benzene']",
    ]
    current_path = os.path.dirname(os.path.realpath(__file__))
    data = pd.read_csv(os.path.join(current_path, "BT_NRTL_dataset.csv"))

    def SSE(model, data):
        expr = (
            float(data["vap_benzene"])
            - model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"]
        ) ** 2 + (
            float(data["liq_benzene"])
            - model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"]
        ) ** 2
        return expr * 1e4

    results = quantify_propagate_uncertainty(
        NRTL_model, NRTL_model_opt, data, variable_name, SSE
    )
