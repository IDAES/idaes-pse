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
""" Tests for data_manager.
"""

import pytest
import pyomo.environ as pyo
from collections import OrderedDict
from pyomo.core.base.componentuid import ComponentUID
from pyomo.common.collections import ComponentMap
import pandas as pd
from idaes.apps.caprese.dynamic_block import DynamicBlock
from idaes.apps.caprese.controller import ControllerBlock
from idaes.apps.caprese.estimator import EstimatorBlock
from idaes.apps.caprese.dynamic_var import DiffVar
from idaes.apps.caprese.tests.test_simple_model import (
        make_model,
        initialize_t0,
        )
from idaes.apps.caprese.data_manager import (
                    empty_dataframe_from_variables,
                    add_variable_values_to_dataframe,
                    add_variable_setpoints_to_dataframe,
                    PlantDataManager,
                    ControllerDataManager,
                    EstimatorDataManager,
                    )

__author__ = "Kuan-Han Lin"

class TestDataManager(object):
    @pytest.mark.unit
    def make_small_ContreteModel(self):
        m = pyo.ConcreteModel()
        m.set1 = pyo.Set(initialize = [0.0, 1.0, 2.0, 3.0, 4.0])
        m.set2 = pyo.Set(initialize = ["A", "B"])
        init_d1 = {(i, "A"):0.1+1.*i for i in m.set1}
        init_d2 = {(i, "B"):0.2+1.*i for i in m.set1}
        init_dict = {**init_d1, **init_d2}
        m.var1 = pyo.Var(m.set1, m.set2, initialize = init_dict)
        return m

    @pytest.mark.unit
    def test_empty_dataframe_from_variables(self):
        m = self.make_small_ContreteModel()

        # w/o rename_map
        variables = [pyo.Reference(m.var1[:, "A"]), pyo.Reference(m.var1[:, "B"])]
        df1 = empty_dataframe_from_variables(variables)
        assert type(df1) is pd.core.frame.DataFrame
        assert df1.empty
        assert len(df1.index) == 0
        pred_column_strings1 = ["iteration"]+\
                                [str(ComponentUID(var.referent)) for var in variables]
        assert all(i1 == i2 for i1, i2 in zip(pred_column_strings1, df1.columns))
        assert df1["iteration"].dtype == "int64"
        assert all(dtype == "float64" for dtype in df1.dtypes[1:])


        # w/ rename_map
        rename_map = {var: str(ComponentUID(var.referent))+"_random"
                      for var in variables}
        df2 = empty_dataframe_from_variables(variables, rename_map)
        assert type(df2) is pd.core.frame.DataFrame
        assert df2.empty
        assert len(df2.index) == 0
        pred_column_strings2 = ["iteration"]+\
                                [str(ComponentUID(var.referent))+"_random" for var in variables]
        assert all(i1 == i2 for i1, i2 in zip(pred_column_strings2, df2.columns))
        assert df2["iteration"].dtype == "int64"
        assert all(dtype == "float64" for dtype in df2.dtypes[1:])

    @pytest.mark.unit
    def test_add_variable_values_to_dataframe(self):
        m = self.make_small_ContreteModel()
        variables = [pyo.Reference(m.var1[:, "A"]), pyo.Reference(m.var1[:, "B"])]
        str_cuids = [str(ComponentUID(var.referent)) for var in variables]

        # Case1: Neither optional argument is given.
        df1 = empty_dataframe_from_variables(variables)
        df1 = add_variable_values_to_dataframe(df1,
                                               variables,
                                               iteration = 10,
                                               )
        assert len(df1.index) == 5
        assert all(i1 == i2 for i1, i2 in zip(df1.index, [0.0, 1.0, 2.0, 3.0, 4.0]))
        assert len(df1.columns) == 3
        assert all(ele == 10 for ele in df1["iteration"])
        assert all(ele == 0.1+ind*1. for ind, ele in enumerate(df1[str_cuids[0]]))
        assert all(ele == 0.2+ind*1. for ind, ele in enumerate(df1[str_cuids[1]]))

        # Case2: time_subset is given
        df2 = empty_dataframe_from_variables(variables)
        df2 = add_variable_values_to_dataframe(df2,
                                               variables,
                                               iteration = 100,
                                               time_subset = [0,1]
                                               )
        assert len(df2.index) == 2
        assert all(i1 == i2 for i1, i2 in zip(df2.index, [0.0,1.0]))
        assert len(df2.columns) == 3
        assert all(ele == 100 for ele in df2["iteration"])
        assert all(ele == 0.1+ind*1. for ind, ele in enumerate(df2[str_cuids[0]]))
        assert all(ele == 0.2+ind*1. for ind, ele in enumerate(df2[str_cuids[1]]))

        df2 = add_variable_values_to_dataframe(df2,
                                               variables,
                                               iteration = 200,
                                               time_subset = [2,3,4]
                                               )
        assert len(df2.index) == 5
        assert all(i1 == i2 for i1, i2 in zip(df2.index, [0.0, 1.0, 3.0, 4.0, 5.0]))
        assert len(df2.columns) == 3
        assert all(ele == 100 for ele in df2["iteration"][0:2])
        assert all(ele == 200 for ele in df2["iteration"][2:])
        assert all(ele == 0.1+ind*1. for ind, ele in enumerate(df2[str_cuids[0]]))
        assert all(ele == 0.2+ind*1. for ind, ele in enumerate(df2[str_cuids[1]]))

        #Case3: rename_map is given
        rename_map = {var: str(ComponentUID(var.referent))+"_random"
                      for var in variables}
        df3 = empty_dataframe_from_variables(variables, rename_map)
        df3 = add_variable_values_to_dataframe(df3,
                                               variables,
                                               iteration = 300,
                                               time_subset = [0,1],
                                               rename_map = rename_map)
        assert len(df3.index) == 2
        assert all(i1 == i2 for i1, i2 in zip(df3.index, [0.0,1.0]))
        assert len(df3.columns) == 3
        assert all(ele == 300 for ele in df3["iteration"])
        assert all(ele == 0.1+ind*1. for ind, ele in enumerate(df3[str_cuids[0]+"_random"]))
        assert all(ele == 0.2+ind*1. for ind, ele in enumerate(df3[str_cuids[1]+"_random"]))

        #Case4: time_map is given
        df5 = empty_dataframe_from_variables(variables)
        time_map = {t: t + 40. for t in m.set1}
        df5 = add_variable_values_to_dataframe(df5,
                                               variables,
                                               iteration = 500,
                                               time_map = time_map)
        assert len(df5.index) == 5
        assert all(i1 == i2 for i1, i2 in zip(df5.index, [40.0, 41.0, 42.0, 43.0, 44.0]))
        assert len(df1.columns) == 3
        assert all(ele == 500 for ele in df5["iteration"])
        assert all(ele == 0.1+ind*1. for ind, ele in enumerate(df5[str_cuids[0]]))
        assert all(ele == 0.2+ind*1. for ind, ele in enumerate(df5[str_cuids[1]]))

    @pytest.mark.unit
    def test_add_variable_setpoints_to_dataframe(self):
        m = self.make_small_ContreteModel()
        variables = [pyo.Reference(m.var1[:, "A"]), pyo.Reference(m.var1[:, "B"])]
        str_cuids = [str(ComponentUID(var.referent)) for var in variables]

        # Manually create differential_vars, save differential varialbes and,
        # assign their setpoints.
        m.differential_vars = []
        for ind, var in enumerate(variables):
            new_ref = pyo.Reference(var.referent, ctype = DiffVar)
            new_ref.setpoint = 10.*(ind+1)
            m.differential_vars.append(new_ref)

        # This mapping is important, because user given variables are not nmpc_var.
        # Thus, they don't have setpoint attribute.
        # This mapping maps the user given variables to corresponding nmpc_vars.
        m.var_mapping = ComponentMap((original_ref, new_ref)
                                     for original_ref, new_ref in zip(variables, m.differential_vars))


        df1 = empty_dataframe_from_variables(variables)
        df1 = add_variable_setpoints_to_dataframe(df1,
                                                  variables,
                                                  time_subset = [0.0, 1.0, 2.0],
                                                  map_for_user_given_vars = m.var_mapping,
                                                  iteration = 10,
                                                  )

        assert all(val == 10 for val in df1["iteration"])
        assert all(val == 10.0 for val in df1["var1[*,A]"])
        assert all(val == 20.0 for val in df1["var1[*,B]"])

        # Set new setpoint value
        for ind, var in enumerate(m.differential_vars):
            var.setpoint = 100.*(ind+1)

        df1 = add_variable_setpoints_to_dataframe(df1,
                                                  variables,
                                                  time_subset = [3.0, 4.0],
                                                  map_for_user_given_vars = m.var_mapping,
                                                  iteration = 20,
                                                  )

        assert all(val == 20 for val in df1["iteration"][3:])
        assert all(val == 100.0 for val in df1["var1[*,A]"][3:])
        assert all(val == 200.0 for val in df1["var1[*,B]"][3:])

    @pytest.mark.unit
    def make_plant(self):
        model = make_model(horizon = 0.5, nfe = 2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements=[model.conc[0,'A'], model.conc[0,'B']]
        plant = DynamicBlock(
                model=model,
                time=time,
                inputs=inputs,
                measurements=measurements,
                )
        plant.construct()
        plant.set_sample_time(sample_time = 0.5)
        return plant

    @pytest.mark.unit
    def make_controller(self):
        model = make_model(horizon = 1., nfe = 4)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements=[model.conc[0,'A'], model.conc[0,'B']]
        controller = ControllerBlock(
                    model=model,
                    time=time,
                    inputs=inputs,
                    measurements=measurements,
                    )
        controller.construct()
        controller.set_sample_time(sample_time = 0.5)
        return controller

    @pytest.mark.unit
    def make_estimator(self):
        model = make_model(horizon = 1., nfe = 4)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements=[model.conc[0,'A'], model.conc[0,'B']]
        estimator = EstimatorBlock(
                    model=model,
                    time=time,
                    inputs=inputs,
                    measurements=measurements,
                    sample_time = 0.5)
        estimator.construct()
        estimator.set_sample_time(sample_time = 0.5)
        return estimator

    @pytest.mark.unit
    def test_PlantDataManager(self):
        model = make_model(horizon = 0.5, nfe = 2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements=[model.conc[0,'A'], model.conc[0,'B']]
        plant = DynamicBlock(
                model=model,
                time=time,
                inputs=inputs,
                measurements=measurements,
                )
        plant.construct()

        states_of_interest = [pyo.Reference(model.conc[:, "A"]),
                              pyo.Reference(model.rate[:, "A"])]
        plant_datamanager = PlantDataManager(plant, states_of_interest,)

        # Make sure all methods are there.
        assert hasattr(plant_datamanager, "__init__")
        assert hasattr(plant_datamanager, "get_plant_dataframe")
        assert hasattr(plant_datamanager, "save_initial_plant_data")
        assert hasattr(plant_datamanager, "save_plant_data")

    @pytest.mark.unit
    def test_PlantDataManager__init__(self):
        plant = self.make_plant()
        model = plant.mod
        time = plant.time
        t0 = time.first()

        states_of_interest = [pyo.Reference(model.conc[:, "A"]),
                              pyo.Reference(model.rate[:, "A"])]
        plant_datamanager = PlantDataManager(plant, states_of_interest)

        assert hasattr(plant_datamanager, "plantblock")
        assert plant == plant_datamanager.plantblock

        # Note that model.conc[:, "A"] has alrealy existed in differential_vars
        assert hasattr(plant_datamanager, "plant_states_of_interest")
        assert all(i1[t0] == i2[t0] for i1, i2 in
                   zip(plant_datamanager.plant_states_of_interest,
                       plant.differential_vars + [states_of_interest[1]]))

        assert hasattr(plant_datamanager, "plant_vars_of_interest")
        pred_plant_vars_of_interest = plant.differential_vars + \
                                        [states_of_interest[1]] + \
                                            plant.input_vars
        assert all(i1[t0] == i2[t0] for i1, i2 in
                   zip(plant_datamanager.plant_vars_of_interest,
                       pred_plant_vars_of_interest))

        assert hasattr(plant_datamanager, "plant_df")
        # empty_dataframe_from_variables has been tested.

    @pytest.mark.unit
    def test_get_plant_dataframe(self):
        plant = self.make_plant()
        model = plant.mod
        time = plant.time
        t0 = time.first()

        states_of_interest = [pyo.Reference(model.conc[:, "A"]),
                              pyo.Reference(model.rate[:, "A"])]
        plant_datamanager = PlantDataManager(plant, states_of_interest)


        assert id(plant_datamanager.get_plant_dataframe()) == id(plant_datamanager.plant_df)

    @pytest.mark.unit
    def test_save_initial_plant_data(self):
        plant = self.make_plant()
        model = plant.mod
        time = plant.time
        t0 = time.first()

        states_of_interest = [pyo.Reference(model.conc[:, "A"]),
                              pyo.Reference(model.rate[:, "A"])]
        plant_datamanager = PlantDataManager(plant, states_of_interest)

        # Set some values at t0 for differential variables
        plant.vectors.differential[0,:].set_value(35.)
        plant.vectors.differential[1,:].set_value(45.)

        # Set value at t0 for the user-interested variable.
        plant.vectors.algebraic[1,:].set_value(0.2)

        plant_datamanager.save_initial_plant_data()

        df = plant_datamanager.plant_df
        assert len(df.index) == 1
        assert df.index[0] == 0.0

        pred_columns = [
            "iteration", "mod.conc[*,A]", "mod.conc[*,B]",
            "mod.rate[*,A]", "mod.flow_in[*]"
        ]
        assert all(i1 == i2 for i1, i2 in zip(pred_columns, df.columns))
        assert df["iteration"][0] == 0
        assert df["mod.conc[*,A]"][0] == 35.
        assert df["mod.conc[*,B]"][0] == 45.
        assert df["mod.rate[*,A]"][0] == 0.2
        import numpy as np
        assert np.isnan(df["mod.flow_in[*]"][0])

    @pytest.mark.unit
    def test_save_plant_data(self):
        plant = self.make_plant()
        model = plant.mod
        time = plant.time
        t0 = time.first()

        states_of_interest = [pyo.Reference(model.conc[:, "A"]),
                              pyo.Reference(model.rate[:, "A"])]
        plant_datamanager = PlantDataManager(plant, states_of_interest)

        # Set some initial values for differential variables
        plant.vectors.differential[0,0].set_value(35.)
        plant.vectors.differential[1,0].set_value(45.)

        # Set values for the user-interested variable
        plant.vectors.algebraic[1,0].set_value(0.2)

        plant_datamanager.save_initial_plant_data()

        # Set values at other time points
        plant.vectors.differential[0,:].set_value(55.)
        plant.vectors.differential[1,:].set_value(65.)
        plant.vectors.algebraic[1,:].set_value(0.4)
        plant.vectors.input[...].set_value(123.)

        plant_datamanager.save_plant_data(iteration = 1)

        df = plant_datamanager.plant_df
        assert len(df.index) == 5
        pred_index = [0.0, 0.083333, 0.25, 0.333333, 0.5]
        assert all(i1 == i2 for i1, i2 in zip(pred_index, df.index))

        pred_columns = [
            "iteration", "mod.conc[*,A]", "mod.conc[*,B]",
            "mod.rate[*,A]", "mod.flow_in[*]"
        ]
        assert all(i1 == i2 for i1, i2 in zip(pred_columns, df.columns))
        assert all(val == 1 for val in df["iteration"][1:])
        assert all(val == 55. for val in df["mod.conc[*,A]"][1:])
        assert all(val == 65. for val in df["mod.conc[*,B]"][1:])
        assert all(val == 0.4 for val in df["mod.rate[*,A]"][1:])
        assert all(val == 123. for val in df["mod.flow_in[*]"][1:])

    @pytest.mark.unit
    def test_ControllerDataManager(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        controller = self.make_controller()
        c_model = controller.mod
        c_time = controller.time
        c_t0 = c_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        nmpc_data = ControllerDataManager(controller,
                                          states_of_interest,)

        # Make sure all methods are there.
        assert hasattr(nmpc_data, "get_controller_dataframe")
        assert hasattr(nmpc_data, "get_setpoint_dataframe")
        assert hasattr(nmpc_data, "save_controller_data")

    @pytest.mark.unit
    def test_ControllerDataManager__init__(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        controller = self.make_controller()
        c_model = controller.mod
        c_time = controller.time
        c_t0 = c_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        nmpc_data = ControllerDataManager(controller,
                                          states_of_interest,)

        assert hasattr(nmpc_data, "controllerblock")
        assert nmpc_data.controllerblock == controller

        assert hasattr(nmpc_data, "controller_states_of_interest")
        pred_controller_states_of_interest = controller.differential_vars + \
                                                [pyo.Reference(c_model.rate[:, "A"])]
        assert all(i1[c_t0] == i2[c_t0] for i1, i2 in
                       zip(pred_controller_states_of_interest,
                           nmpc_data.controller_states_of_interest))

        assert hasattr(nmpc_data, "controller_df")

        assert hasattr(nmpc_data, "setpoint_df")

        assert hasattr(nmpc_data, "user_given_vars_map_nmpcvar")
        vardata_map = controller.vardata_map
        for var in nmpc_data.extra_vars_user_interested:
            assert nmpc_data.user_given_vars_map_nmpcvar[var] == vardata_map[var[c_t0]]

    @pytest.mark.unit
    def test_get_controller_dataframe(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        controller = self.make_controller()
        c_model = controller.mod
        c_time = controller.time
        c_t0 = c_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        nmpc_data = ControllerDataManager(controller,
                                          states_of_interest,)

        assert id(nmpc_data.get_controller_dataframe()) == id(nmpc_data.controller_df)

    @pytest.mark.unit
    def test_get_setpoint_dataframe(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        controller = self.make_controller()
        c_model = controller.mod
        c_time = controller.time
        c_t0 = c_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        nmpc_data = ControllerDataManager(controller,
                                          states_of_interest,)

        assert id(nmpc_data.get_setpoint_dataframe()) == id(nmpc_data.setpoint_df)

    @pytest.mark.unit
    def test_save_controller_data(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        controller = self.make_controller()
        c_model = controller.mod
        c_time = controller.time
        c_t0 = c_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        nmpc_data = ControllerDataManager(controller, states_of_interest)

        # Set values of control inputs and setpoints in the controller
        controller.vectors.input[...].set_value(5.0)
        controller.vectors.differential.set_setpoint([0.2, 0.3])
        controller.vectors.algebraic.set_setpoint([4.5, 5.5, 6.5])

        nmpc_data.save_controller_data(iteration = 10)

        pred_controller_result = OrderedDict(
            [("iteration", 10), ("mod.flow_in[*]", 5)]
        )
        pred_setpoint_result = OrderedDict([("iteration", 10),
                                            ("mod.conc[*,A]", 0.2),
                                            ("mod.conc[*,B]", 0.3),
                                            ("mod.rate[*,A]", 5.5)])

        pred_index = [0.083333, 0.25, 0.333333, 0.5]

        df = nmpc_data.controller_df
        assert all(i1 == i2 for i1, i2 in
                zip(pred_controller_result.keys(), df.columns))
        assert all(i1 == i2 for i1, i2 in zip(pred_index, df.index))
        for key, value in pred_controller_result.items():
            assert all(val == value for val in df[key])

        df_sp = nmpc_data.setpoint_df
        assert all(i1 == i2 for i1, i2 in zip(pred_setpoint_result.keys(), df_sp.columns))
        assert all(i1 == i2 for i1, i2 in zip(pred_index, df_sp.index))
        for key, value in pred_setpoint_result.items():
            assert all(val == value for val in df_sp[key])

        # Check for the other iteration
        # re-set values
        controller.vectors.input[...].set_value(15.0)
        controller.vectors.differential.set_setpoint([5.2, 5.3])
        controller.vectors.algebraic.set_setpoint([9.5, 10.5, 11.5])

        nmpc_data.save_controller_data(iteration = 20)

        pred_columns_result = OrderedDict(
            [("iteration", 20), ("mod.flow_in[*]", 15)]
        )
        pred_setpoint_result = OrderedDict([("iteration", 20),
                                            ("mod.conc[*,A]", 5.2),
                                            ("mod.conc[*,B]", 5.3),
                                            ("mod.rate[*,A]", 10.5)])
        pred_index = [0.583333, 0.75, 0.833333, 1.0]

        df = nmpc_data.controller_df
        assert all(i1 == i2 for i1, i2 in
                zip(pred_columns_result.keys(), df.columns))
        assert all(i1 == i2 for i1, i2 in zip(pred_index, df.index[4:]))
        for key, value in pred_columns_result.items():
            assert all(val == value for val in df[key][4:])

        df_sp = nmpc_data.setpoint_df
        assert all(i1 == i2 for i1, i2 in
                zip(pred_setpoint_result.keys(), df_sp.columns))
        assert all(i1 == i2 for i1, i2 in zip(pred_index, df_sp.index[4:]))
        for key, value in pred_setpoint_result.items():
            assert all(val == value for val in df_sp[key][4:])

    @pytest.mark.unit
    def test_EstimatorDataManager(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        estimator = self.make_estimator()
        e_model = estimator.mod
        e_time = estimator.time
        e_t0 = e_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"])]
        mhe_data = EstimatorDataManager(estimator,
                                        states_of_interest,)

        # Make sure all methods are there.
        assert hasattr(mhe_data, "get_estimator_dataframe")
        assert hasattr(mhe_data, "save_estimator_data")

    @pytest.mark.unit
    def test_EstimatorDataManager__init__(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        estimator = self.make_estimator()
        e_model = estimator.mod
        e_time = estimator.time
        e_t0 = e_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        mhe_data = EstimatorDataManager(estimator, states_of_interest,)

        assert hasattr(mhe_data, "estimatorblock")
        assert mhe_data.estimatorblock == estimator

        assert hasattr(mhe_data, "estimator_vars_of_interest")
        pred_estimator_vars_of_interest = estimator.differential_vars + \
                                            [pyo.Reference(e_model.rate[:, "A"])]
        assert all(i1[e_t0] == i2[e_t0] for i1, i2
                    in zip(pred_estimator_vars_of_interest,
                            mhe_data.estimator_vars_of_interest))

        assert hasattr(mhe_data, "estimator_df")

    @pytest.mark.unit
    def test_get_estimator_dataframe(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        estimator = self.make_estimator()
        e_model = estimator.mod
        e_time = estimator.time
        e_t0 = e_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        mhe_data = EstimatorDataManager(estimator, states_of_interest)

        df = mhe_data.get_estimator_dataframe()
        assert id(df) ==  id(mhe_data.estimator_df)

    @pytest.mark.unit
    def test_save_estimator_data(self):
        plant = self.make_plant()
        p_model = plant.mod
        p_time = plant.time
        p_t0 = p_time.first()

        estimator = self.make_estimator()
        e_model = estimator.mod
        e_time = estimator.time
        e_t0 = e_time.first()

        states_of_interest = [pyo.Reference(p_model.conc[:, "A"]),
                              pyo.Reference(p_model.rate[:, "A"])]
        mhe_data = EstimatorDataManager(estimator, states_of_interest,)

        # Set values for differential variables
        estimator.vectors.differential[0,:].set_value(100.)
        estimator.vectors.differential[1,:].set_value(200.)
        estimator.vectors.algebraic[1,:].set_value(30.)

        mhe_data.save_estimator_data(iteration = 10)

        pred_column_results = OrderedDict([("iteration", 10),
                                            ("mod.conc[*,A]", 100.),
                                            ("mod.conc[*,B]", 200.),
                                            ("mod.rate[*,A]", 30.)])
        df = mhe_data.estimator_df
        assert all(i1 == i2 for i1, i2 in zip(pred_column_results.keys(), df.columns))
        assert df.index[0] == 5.5
        for key, value in pred_column_results.items():
            assert df[key][5.5] == value

        # check it for the other iteration
        # Re-set values for differential variables
        estimator.vectors.differential[0,:].set_value(300.)
        estimator.vectors.differential[1,:].set_value(400.)
        estimator.vectors.algebraic[1,:].set_value(60.)

        mhe_data.save_estimator_data(iteration = 20)

        pred_column_results = OrderedDict([("iteration", 20),
                                            ("mod.conc[*,A]", 300.),
                                            ("mod.conc[*,B]", 400.),
                                            ("mod.rate[*,A]", 60.)])
        df = mhe_data.estimator_df
        assert all(i1 == i2 for i1, i2 in zip(pred_column_results.keys(), df.columns))
        assert df.index[1] == 10.5
        for key, value in pred_column_results.items():
            assert df[key][10.5] == value
