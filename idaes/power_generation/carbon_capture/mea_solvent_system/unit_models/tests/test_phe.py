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
"""
Test for plate heat exchanger model

Author: Paul Akula
"""
# Import Python libraries
import pytest

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value, TerminationCondition

# Import IDAES Libraries
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver
from idaes.power_generation.carbon_capture.mea_solvent_system.properties.liquid_prop \
    import LiquidParameterBlock
from idaes.power_generation.carbon_capture.mea_solvent_system.unit_models.phe \
    import PHE


solver = get_solver()


class TestPlateHeatExchanger:
    '''
    Tests for the plate heat exchange (PHE) model

    Inputs for states variables  are in SI units
    -Flowrate: mol/s
    -Temperature: K
    -Pressure: Pa
    '''

    @pytest.fixture(scope="module")
    def measurement(self):
        nccc_data = {'dataset1':
                     {'input':
                      {'hotside':
                       {'flowrate': 60.54879,
                        'temperature': 392.23,
                        'pressure': 202650,
                        'mole_fraction':
                        {'CO2': 0.0158,
                         'MEA': 0.1095,
                         'H2O': 0.8747}},
                       'coldside':
                       {'flowrate': 63.01910,
                        'temperature': 326.36,
                        'pressure': 202650,
                        'mole_fraction':
                        {'CO2': 0.0414,
                         'MEA': 0.1077,
                         'H2O': 0.8509}}},
                      'output':
                      {'hotside':
                       {'temperature': 330.42},
                       'coldside':
                       {'temperature': 384.91}}},
                     'dataset2':
                     {"input":
                      {'hotside':
                       {'flowrate': 102.07830,
                        'temperature': 389.57,
                        'pressure': 202650,
                        'mole_fraction':
                        {'CO2': 0.0284,
                         'MEA': 0.1148,
                         'H2O': 0.8568}},
                       'coldside':
                       {'flowrate': 104.99350,
                        'temperature': 332.26,
                        'pressure': 202650,
                        'mole_fraction':
                        {'CO2': 0.0438,
                         'MEA': 0.1137,
                         'H2O': 0.8425}}},
                      'output':
                      {'hotside':
                       {'temperature': 336.70},
                       'coldside':
                       {'temperature': 383.21}}}}

        return nccc_data

    @pytest.fixture(scope="class")
    def phe_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        # Set up property package
        m.fs.hotside_properties = LiquidParameterBlock()
        m.fs.coldside_properties = LiquidParameterBlock()

        # create instance of plate heat exchanger  on flowsheet
        m.fs.unit = PHE(default={'passes': 4,
                                 'channel_list': [12, 12, 12, 12],
                                 'divider_plate_number': 2,
                                 "hot_side": {
                                     "property_package": m.fs.hotside_properties
                                 },
                                 "cold_side": {
                                     "property_package": m.fs.coldside_properties
                                 }})
        return m

    @pytest.fixture(scope="class", params=['dataset1', 'dataset2'])
    def set_input_output(self, phe_model, measurement, request):
        for t in phe_model.fs.time:
            # hot fluid
            Th_out = measurement[request.param]['output']['hotside']['temperature']
            Th = measurement[request.param]['input']['hotside']['temperature']
            Fh = measurement[request.param]['input']['hotside']['flowrate']
            Ph = measurement[request.param]['input']['hotside']['pressure']
            xh_CO2 = measurement[request.param]['input']['hotside']['mole_fraction']['CO2']
            xh_MEA = measurement[request.param]['input']['hotside']['mole_fraction']['MEA']
            xh_H2O = measurement[request.param]['input']['hotside']['mole_fraction']['H2O']

            phe_model.fs.unit.hot_inlet.flow_mol[t].fix(Fh)
            phe_model.fs.unit.hot_inlet.temperature[t].fix(Th)
            phe_model.fs.unit.hot_inlet.pressure[t].fix(Ph)
            phe_model.fs.unit.hot_inlet.mole_frac_comp[t, "CO2"].fix(xh_CO2)
            phe_model.fs.unit.hot_inlet.mole_frac_comp[t, "H2O"].fix(xh_H2O)
            phe_model.fs.unit.hot_inlet.mole_frac_comp[t, "MEA"].fix(xh_MEA)

            # cold fluid
            Tc_out = measurement[request.param]['output']['coldside']['temperature']
            Tc = measurement[request.param]['input']['coldside']['temperature']
            Fc = measurement[request.param]['input']['coldside']['flowrate']
            Pc = measurement[request.param]['input']['coldside']['pressure']
            xc_CO2 = measurement[request.param]['input']['coldside']['mole_fraction']['CO2']
            xc_MEA = measurement[request.param]['input']['coldside']['mole_fraction']['MEA']
            xc_H2O = measurement[request.param]['input']['coldside']['mole_fraction']['H2O']
            phe_model.fs.unit.cold_inlet.flow_mol[t].fix(Fc)
            phe_model.fs.unit.cold_inlet.temperature[t].fix(Tc)
            phe_model.fs.unit.cold_inlet.pressure[t].fix(Pc)
            phe_model.fs.unit.cold_inlet.mole_frac_comp[t, "CO2"].fix(xc_CO2)
            phe_model.fs.unit.cold_inlet.mole_frac_comp[t, "H2O"].fix(xc_H2O)
            phe_model.fs.unit.cold_inlet.mole_frac_comp[t, "MEA"].fix(xc_MEA)

            output = {'hotside_temperature_expectation': Th_out,
                      'coldside_temperature_expectation': Tc_out}
            return output

    @pytest.mark.unit
    def test_phe_build(self, phe_model):
        assert phe_model.fs.unit.config.dynamic is False
        assert degrees_of_freedom(phe_model) == 12

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, phe_model):
        for t in phe_model.fs.time:
            # hot fluid
            phe_model.fs.unit.hot_inlet.flow_mol[t].fix(60.54879)
            phe_model.fs.unit.hot_inlet.temperature[t].fix(392.23)
            phe_model.fs.unit.hot_inlet.pressure[t].fix(202650)
            phe_model.fs.unit.hot_inlet.mole_frac_comp[t, "CO2"].fix(0.0158)
            phe_model.fs.unit.hot_inlet.mole_frac_comp[t, "H2O"].fix(0.8747)
            phe_model.fs.unit.hot_inlet.mole_frac_comp[t, "MEA"].fix(0.1095)

            # cold fluid
            phe_model.fs.unit.cold_inlet.flow_mol[t].fix(63.01910)
            phe_model.fs.unit.cold_inlet.temperature[t].fix(326.36)
            phe_model.fs.unit.cold_inlet.pressure[t].fix(202650)
            phe_model.fs.unit.cold_inlet.mole_frac_comp[t, "CO2"].fix(0.0414)
            phe_model.fs.unit.cold_inlet.mole_frac_comp[t, "H2O"].fix(0.1077)
            phe_model.fs.unit.cold_inlet.mole_frac_comp[t, "MEA"].fix(0.8509)
        initialization_tester(phe_model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_phe_validation(self, phe_model, set_input_output):

        phe_model.fs.unit.initialize()
        res = solver.solve(phe_model)

        assert res.solver.termination_condition == TerminationCondition.optimal

        # Mass conservation test
        assert abs(value(phe_model.fs.unit.hot_inlet.flow_mol[0] -
                         phe_model.fs.unit.hot_outlet.flow_mol[0])) <= 1e-6

        assert abs(value(phe_model.fs.unit.cold_inlet.flow_mol[0] -
                         phe_model.fs.unit.cold_outlet.flow_mol[0])) <= 1e-6

        # Energy conservation test
        assert value(phe_model.fs.unit.QH[0]) == value(phe_model.fs.unit.QC[0])

        # Exit tempreture prediction test
        assert value(phe_model.fs.unit.hot_outlet.temperature[0]) == \
            pytest.approx(
            set_input_output['hotside_temperature_expectation'], abs=1.5)

        assert value(phe_model.fs.unit.cold_outlet.temperature[0]) == \
            pytest.approx(
            set_input_output['coldside_temperature_expectation'], abs=1.5)
