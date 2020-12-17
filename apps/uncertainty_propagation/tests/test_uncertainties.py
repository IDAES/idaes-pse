##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
sys.path.append(os.path.abspath('..')) # current folder is ~/tests
import numpy as np
import pandas as pd
import pytest
from mock import patch
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_unucertainty
from pyomo.opt import SolverFactory
ipopt_available = SolverFactory('ipopt').available()
kaug_available = SolverFactory('k_aug').available()
dotsens_available = SolverFactory('dot_sens').available()

class TestUncertaintyPropagation:

    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_theta_cov_propagation1(self):
        from idaes.apps.uncertainty_propagation.tests.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
    
        obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE)

        np.testing.assert_almost_equal(obj, 4.331711213656886)
        np.testing.assert_almost_equal(theta[variable_name[0]], 19.142575284617866)
        np.testing.assert_almost_equal(theta[variable_name[1]], 0.53109137696521)
        np.testing.assert_almost_equal(cov, np.array([[6.30579403, -0.4395341], [-0.4395341, 0.04193591]]))
        np.testing.assert_almost_equal(propagation_f['objective'], 5.45439337747349)
        assert propagation_c == {}
    
    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_theta_cov_propagation2(self):
        from idaes.apps.uncertainty_propagation.tests.NRTL_model_scripts import NRTL_model, NRTL_model_opt
        variable_name = ["fs.properties.tau[benzene,toluene]", "fs.properties.tau[toluene,benzene]"]
        data = pd.read_csv('BT_NRTL_dataset.csv')
        def SSE(model, data):
            expr = ((float(data["vap_benzene"]) -
                     model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"])**2 +
                    (float(data["liq_benzene"]) -
                     model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"])**2)
            return expr*1E4
        obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(NRTL_model,NRTL_model_opt, data, variable_name, SSE)

        np.testing.assert_almost_equal(obj, 0.004663348837044143)
        np.testing.assert_almost_equal(theta[variable_name[0]],  0.4781086784101746)
        np.testing.assert_almost_equal(theta[variable_name[1]],  -0.40924465377598657)
        np.testing.assert_almost_equal(cov, np.array([[0.00213426, -0.00163064], [-0.00163064, 0.00124591]]))
        np.testing.assert_almost_equal(propagation_f['objective'], 0.00014333989649382864)
        np.testing.assert_almost_equal(propagation_c['constraints 4'], 0.000084167318885)
        np.testing.assert_almost_equal(propagation_c['constraints 5'], 0.0002455439364710466)
        np.testing.assert_almost_equal(propagation_c['constraints 6'], 0.0008215307761164211)
        np.testing.assert_almost_equal(propagation_c['constraints 7'], 0.0001749469866777253)
        np.testing.assert_almost_equal(propagation_c['constraints 8'], 0.0005214685823573828)
        np.testing.assert_almost_equal(propagation_c['constraints 9'], 0.00024951071782077465)

    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_Exception1(self):
        from idaes.apps.uncertainty_propagation.tests.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        tee = 1
        with pytest.raises(Exception):
            obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE,tee)


    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_Exception2(self):
        from idaes.apps.uncertainty_propagation.tests.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        tee = False
        diagnostic_mode = 1
        with pytest.raises(Exception):
            obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE,tee,diagnostic_mode)


    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_Exception3(self):
        from idaes.apps.uncertainty_propagation.tests.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        tee = False
        diagnostic_mode = False
        solver_options = [1e-8]
        with pytest.raises(Exception):
            obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE,tee,diagnostic_mode,solver_options)

