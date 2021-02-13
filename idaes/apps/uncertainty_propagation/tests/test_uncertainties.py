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
from pytest import approx
from mock import patch
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_uncertainty, propagate_uncertainty, get_sensitivity, clean_variable_name
from pyomo.opt import SolverFactory
from pyomo.environ import *
import pyomo.contrib.parmest.parmest as parmest

ipopt_available = SolverFactory('ipopt').available()
kaug_available = SolverFactory('k_aug').available()
dotsens_available = SolverFactory('dot_sens').available()

class TestUncertaintyPropagation:

    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_quantify_propagate_uncertainty1(self):
        '''
        It tests the function quantify_propagate_uncertainty with rooney & biegler's model.
        '''
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
    
        results =  quantify_propagate_uncertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE)

        assert results.obj == approx(4.331711213656886)
        assert results.theta_out[variable_name[0]] == approx(19.142575284617866)
        assert results.theta_out[variable_name[1]] == approx(0.53109137696521)
        assert results.propagation_f['objective'] == approx(5.45439337747349)
        assert results.propagation_c == {}
        np.testing.assert_array_almost_equal(results.cov, np.array([[6.30579403, -0.4395341], [-0.4395341, 0.04193591]])) 
        
    def test_quantify_propagate_uncertainty2(self):
        '''
        This is the same test as test_quantify_propagate_uncertainty1,
        but with the second argument of quantify_propagate_uncertainty as Pyomo Concrete Model.
        '''
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        model_uncertain= ConcreteModel()
        model_uncertain.asymptote = Var(initialize = 15)
        model_uncertain.rate_constant = Var(initialize = 0.5)
        model_uncertain.obj = Objective(expr = model_uncertain.asymptote*( 1 - exp(-model_uncertain.rate_constant*10  )  ), sense=minimize)

        results =  quantify_propagate_uncertainty(rooney_biegler_model,model_uncertain, data, variable_name, SSE)

        assert results.obj == approx(4.331711213656886)
        assert results.theta_out[variable_name[0]] == approx(19.142575284617866)
        assert results.theta_out[variable_name[1]] == approx(0.53109137696521)
        assert results.propagation_f['objective'] == approx(5.45439337747349)
        assert results.propagation_c == {}
        np.testing.assert_array_almost_equal(results.cov, np.array([[6.30579403, -0.4395341], [-0.4395341, 0.04193591]]))   
 
    def test_propagate_uncertainty(self):
        '''
        It tests the function propagate_uncertainty with rooney & biegler's model.
        '''
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        parmest_class = parmest.Estimator(rooney_biegler_model, data,variable_name,SSE)
        obj, theta, cov = parmest_class.theta_est(calc_cov=True)
        model_uncertain= ConcreteModel()
        model_uncertain.asymptote = Var(initialize = 15)
        model_uncertain.rate_constant = Var(initialize = 0.5)
        model_uncertain.obj = Objective(expr = model_uncertain.asymptote*( 1 - exp(-model_uncertain.rate_constant*10  )  ), sense=minimize)

        propagation_f, propagation_c =  propagate_uncertainty(model_uncertain, theta, cov, variable_name)

        assert propagation_f['objective'] == approx(5.45439337747349)
        assert propagation_c == {}
        
    def test_get_sensitivity(self):
        '''
        It tests the function get_sensitivity with rooney & biegler's model.
        '''        
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        parmest_class = parmest.Estimator(rooney_biegler_model, data,variable_name,SSE)
        obj, theta, cov = parmest_class.theta_est(calc_cov=True)
        model_uncertain= ConcreteModel()
        model_uncertain.asymptote = Var(initialize = 15)
        model_uncertain.rate_constant = Var(initialize = 0.5)
        model_uncertain.obj = Objective(expr = model_uncertain.asymptote*( 1 - exp(-model_uncertain.rate_constant*10  )  ), sense=minimize)
        theta= {'asymptote': 19.142575284617866, 'rate_constant': 0.53109137696521}
        for v in variable_name:
            getattr(model_uncertain, v).setlb(theta[v])
            getattr(model_uncertain, v).setub(theta[v])
        gradient_f, gradient_c, line_dic =  get_sensitivity(model_uncertain, variable_name)
        
        assert gradient_f == approx(np.array([0.99506259, 0.945148]))
        assert gradient_c == approx(np.array([]))
        assert line_dic['asymptote'] == approx(1)
        assert line_dic['rate_constant'] == approx(2)
        
    
    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_quantify_propagate_uncertainty_NRTL(self):
        '''
        It tests the function quantify_propagate_uncertainty with IDAES NRTL model.
        '''
        from idaes.apps.uncertainty_propagation.examples.NRTL_model_scripts import NRTL_model, NRTL_model_opt
        variable_name = ["fs.properties.tau['benzene','toluene']", "fs.properties.tau['toluene','benzene']"]
        current_path = os.path.dirname(os.path.realpath(__file__))
        data = pd.read_csv(os.path.join(current_path, 'BT_NRTL_dataset.csv'))
        def SSE(model, data):
            expr = ((float(data["vap_benzene"]) -
                     model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"])**2 +
                    (float(data["liq_benzene"]) -
                     model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"])**2)
            return expr*1E4
        results =  quantify_propagate_uncertainty(NRTL_model,NRTL_model_opt, data, variable_name, SSE)
        assert results.obj == approx(5.074968578304798)
        assert results.theta_out[variable_name[0]] == approx(-0.8987624039723433)
        assert results.theta_out[variable_name[1]] == approx(1.4104861106603803)
        np.testing.assert_almost_equal(results.cov, np.array([[0.01194738, -0.02557055], [-0.02557055, 0.05490639]]))
        assert results.propagation_f['objective'] == approx(0.0021199499778127204)
        assert results.propagation_c['constraints 4'] == approx(0.008473482674782129)
        assert results.propagation_c['constraints 5'] == approx(0.0004612914088240769)
        assert results.propagation_c['constraints 6'] == approx(0.003094158491940104)
        assert results.propagation_c['constraints 7'] == approx(0.011260188179282066)
        assert results.propagation_c['constraints 8'] == approx(0.014194289370065335)
        assert results.propagation_c['constraints 9'] == approx( 0.005670725779499914)

    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_quantify_propagate_uncertainty_NRTL_exception(self):
        '''
        It tests an exception error when the ipopt fails for the function quantify_propagate_uncertainty with IDAES NRTL model.
        '''
        from idaes.apps.uncertainty_propagation.examples.NRTL_model_scripts import NRTL_model, NRTL_model_opt_infeasible
        variable_name = ["fs.properties.tau['benzene','toluene']", "fs.properties.tau['toluene','benzene']"]
        current_path = os.path.dirname(os.path.realpath(__file__))
        data = pd.read_csv(os.path.join(current_path, 'BT_NRTL_dataset.csv'))
        def SSE(model, data):
            expr = ((float(data["vap_benzene"]) -
                     model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"])**2 +
                    (float(data["liq_benzene"]) -
                     model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"])**2)
            return expr*1E4
        with pytest.raises(Exception):
            results =  quantify_propagate_uncertainty(NRTL_model,NRTL_model_opt_infeasible, data, variable_name, SSE)


    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_Exception1(self):
        '''
        It tests an ValueError when the tee is not bool for the function quantify_propagate_uncertainty with rooney & biegler's model.
        '''
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        tee = 1
        with pytest.raises(ValueError):
            results =  quantify_propagate_uncertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE,tee)


    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_Exception2(self):
        '''
        It tests an ValueError when the diagnostic_mode is not bool for the function quantify_propagate_uncertainty with rooney & biegler's model.
        '''
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
        variable_name = ['asymptote', 'rate_constant']
        data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],
                                  [4,16.0],[5,15.6],[7,19.8]],
                            columns=['hour', 'y'])
        def SSE(model, data):
            expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
            return expr
        tee = False
        diagnostic_mode = 1
        with pytest.raises(ValueError):
            results =  quantify_propagate_uncertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE,tee,diagnostic_mode)


    @pytest.mark.unit
    @pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'k_aug' command is not available")
    @pytest.mark.skipif(not ipopt_available, reason="The 'dot_sens' command is not available")
    def test_Exception3(self):
        '''
        It tests an ValeError when solver_options is not a dictionary for the function quantify_propagate_uncertainty with rooney & biegler's model.
        '''
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
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
        with pytest.raises(ValueError):
            results =  quantify_propagate_uncertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE,tee,diagnostic_mode,solver_options)

    def test_clean_variable_name1(self):
        '''
        It tests the function clean_variable_name when variable names contain ' and spaces.
        '''
        theta_names = ["fs.properties.tau['benzene', 'toluene']", "fs.properties.tau['toluene', 'benzene' ]"] 
        theta_names_new, var_dic = clean_variable_name(theta_names)
        theta_names_expected = ["fs.properties.tau[benzene,toluene]", "fs.properties.tau[toluene,benzene]"]
        assert len(theta_names_expected) == len(theta_names_new)
        assert all([a == b for a, b in zip(theta_names_expected, theta_names_new)]) 
        
        assert len(theta_names_expected) == len(var_dic.keys())
        assert all([a == b for a, b in zip(sorted(theta_names_expected), sorted(var_dic.keys()))])

        assert len(theta_names) == len(var_dic.values())
        assert all([a == b for a, b in zip(sorted(theta_names), sorted(var_dic.values()))])

    def test_clean_variable_name2(self):
        '''
        It tests the function clean_variable_name when variable names do not contain any ' and spaces.
        '''
        theta_names = ["fs.properties.tau[benzene,toluene]", "fs.properties.tau[toluene,benzene]"]
        theta_names_new, var_dic = clean_variable_name(theta_names)
        assert len(theta_names) == len(theta_names_new)
        assert all([a == b for a, b in zip(theta_names, theta_names_new)])

        assert len(theta_names) == len(var_dic.keys())
        assert all([a == b for a, b in zip(sorted(theta_names), sorted(var_dic.keys()))])

        assert len(theta_names) == len(var_dic.values())
        assert all([a == b for a, b in zip(sorted(theta_names), sorted(var_dic.values()))])
