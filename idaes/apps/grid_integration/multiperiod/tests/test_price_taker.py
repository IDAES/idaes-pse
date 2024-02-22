#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

from pathlib import Path
import pytest
import numpy as np
import pandas as pd

from pyomo.environ import (
    Var,
    Binary,
    NonNegativeReals,
    Constraint,
    Expression,
)

from idaes.apps.grid_integration import DesignModel, OperationModel, deepgetattr

import idaes.logger as idaeslog

from idaes.apps.grid_integration.multiperiod.price_taker_model import PriceTakerModel

import matplotlib.pyplot as plt


@pytest.fixture
def excel_data():
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS.xlsx"
    data = pd.read_excel(file_path, sheet_name=1)
    return data


@pytest.mark.unit
def test_daily_data_size(excel_data):
    m = PriceTakerModel()

    # Generate price data for each hour of every day in the data
    daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])

    # Check that there is a row for each horizon length in a representative day
    assert len(daily_data) == m.horizon_length


@pytest.mark.unit
def test_determine_optimal_num_clusters(excel_data):
    # Test fails because the optimal number of clusters is different
    # based on which version of scikit-learn and kneed you have.
    # Older versions get n_clusters = 10, Newer versions n_clusters = 11
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
    n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data)

    assert n_clusters == 10


@pytest.mark.unit
def test_elbow_plot(excel_data):
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
    m.get_optimal_n_clusters(daily_data, plot=True)

    assert plt.gcf() is not None


@pytest.mark.unit
def test_cluster_lmp_data(excel_data):
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
    n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data)
    lmp_data, weights = m.cluster_lmp_data(excel_data['BaseCaseTax'], n_clusters=n_clusters)

    sum_of_weights = 0
    for i in range(1, n_clusters + 1):
        sum_of_weights = sum_of_weights + weights[0][i]
    assert sum_of_weights == 365

    assert len(lmp_data) == n_clusters


@pytest.mark.unit
def test_init_logger_messages(excel_data, caplog):
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
        m.get_optimal_n_clusters(daily_data)

        assert f"kmax was not set - using a default value of 30." in caplog.text

    # Testing horizon_length input value errors
    value = 0
    with pytest.raises(
        ValueError,
        match=(f"horizon_length must be > 0, but {value} is provided."),
    ):
        m = PriceTakerModel()
        m.horizon_length = value
    value = 12.34
    with pytest.raises(
        ValueError,
        match=(f"horizon_length must be an integer, but {value} is not an integer"),
    ):
        m = PriceTakerModel()
        m.horizon_length = value
    
    # Testing seed errors
    value = 12.34
    with pytest.raises(
        ValueError,
        match=(f"seed must be an integer, but {value} is not an integer"),
    ):
        m = PriceTakerModel()
        m.seed = value


@pytest.mark.unit
def test_min_up_down_time_logger_messages(excel_data):
    des = DesignModel()
    oper = OperationModel()
    build_bin_var = 'build'
    use_min_time=True
    up_time = [10, -5, 2.2]
    down_time = [10, -5, 2.2]
    with pytest.raises(
        ValueError,
        match=(f"up_time must be an integer, but {up_time[2]} is not an integer"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, use_min_time, up_time[2], down_time[0])
    
    with pytest.raises(
        ValueError,
        match=(f"down_time must be an integer, but {down_time[2]} is not an integer"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, use_min_time, up_time[0], down_time[2])
    
    with pytest.raises(
        ValueError,
        match=(f"up_time must be >= 1, but {up_time[1]} is not"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, use_min_time, up_time[1], down_time[0])
    
    with pytest.raises(
        ValueError,
        match=(f"down_time must be >= 1, but {down_time[1]} is not"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, use_min_time, up_time[0], down_time[1])


@pytest.mark.unit
def test_optimal_clusters_logger_messages(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    with pytest.raises(
        ValueError,
        match=(f"kmin must be an integer, but {kmin[1]} is not an integer"),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
        n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data, kmin=kmin[1], kmax=kmax[2])
    
    with pytest.raises(
        ValueError,
        match=(f"kmax must be an integer, but {kmax[1]} is not an integer"),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
        n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data, kmin=kmin[2], kmax=kmax[1])

    with pytest.raises(
        ValueError,
        match=(f"kmin must be > 0, but {kmin[0]} is provided."),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
        n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data, kmin=kmin[0], kmax=kmax[2])

    
    with pytest.raises(
        ValueError,
        match=(f"kmax must be > 0, but {kmax[0]} is provided."),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
        n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data, kmin=kmin[2], kmax=kmax[0])
        m.cluster_lmp_data(excel_data, n_clusters)
    
    with pytest.raises(
        ValueError,
        match=(f"kmin must be less than kmax, but {kmin[2]} >= {kmax[2]}"),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data['BaseCaseTax'])
        n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data, kmin=kmin[2], kmax=kmax[2])
    
    # TODO: The below test is not working because our data doesn't ever seem to arrive at an n_clusters close to kmax
    # caplog.clear()
    # with caplog.at_level(idaeslog.WARNING):
    #     m = PriceTakerModel()
    #     kmin = 1
    #     kmax = 7
    #
    #     daily_data, scenarios = m.reconfigure_raw_data(excel_data)
    #     m.get_optimal_n_clusters(daily_data, kmin=kmin, kmax=kmax)
    #
    #     assert f"Optimal number of clusters is close to kmax: {kmax}. Consider increasing kmax." in caplog.text


@pytest.mark.unit
def test_ramping_constraint_logger_messages(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = 'power'
    capac_var = 'power_max'
    constraint_type = 'linear'
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    with pytest.raises(
        ValueError,
        match=(f"startup_rate fraction must be between 0 and 1, but {su_rate[0]} is not."),
    ):
        m = PriceTakerModel()
        m.add_ramping_constraints(des, oper, capac_var, ramping_var, constraint_type, linearization, 
                                  op_range_lb[1], su_rate[0], sd_rate[1], op_ru_rate[1], op_rd_rate[1])
    
    with pytest.raises(
        ValueError,
        match=(f"shutdown_rate fraction must be between 0 and 1, but {sd_rate[0]} is not."),
    ):
        m = PriceTakerModel()
        m.add_ramping_constraints(des, oper, capac_var, ramping_var, constraint_type, linearization, 
                                  op_range_lb[1], su_rate[1], sd_rate[0], op_ru_rate[1], op_rd_rate[1])
    
    with pytest.raises(
        ValueError,
        match=(f"ramp_up_rate fraction must be between 0 and 1, but {op_ru_rate[0]} is not."),
    ):
        m = PriceTakerModel()
        m.add_ramping_constraints(des, oper, capac_var, ramping_var, constraint_type, linearization, 
                                  op_range_lb[1], su_rate[1], sd_rate[1], op_ru_rate[0], op_rd_rate[1])
    
    with pytest.raises(
        ValueError,
        match=(f"ramp_down_rate fraction must be between 0 and 1, but {op_rd_rate[0]} is not."),
    ):
        m = PriceTakerModel()
        m.add_ramping_constraints(des, oper, capac_var, ramping_var, constraint_type, linearization, 
                                  op_range_lb[1], su_rate[1], sd_rate[1], op_ru_rate[1], op_rd_rate[0])
    
    with pytest.raises(
        ValueError,
        match=(f"op_range_lb fraction must be between 0 and 1, but {op_range_lb[0]} is not."),
    ):
        m = PriceTakerModel()
        m.add_ramping_constraints(des, oper, capac_var, ramping_var, constraint_type, linearization, 
                                  op_range_lb[0], su_rate[1], sd_rate[1], op_ru_rate[1], op_rd_rate[1])
    
    with pytest.raises(
        ValueError,
        match=(f"op_range_lb fraction must be <= shut_down_rate, otherwise the system cannot reach the off state."),
    ):
        m = PriceTakerModel()
        m.add_ramping_constraints(des, oper, capac_var, ramping_var, constraint_type, linearization, 
                                  op_range_lb[2], su_rate[1], sd_rate[1], op_ru_rate[1], op_rd_rate[1])


@pytest.mark.unit
def test_append_lmp_data_logger_messages(excel_data, caplog):
    file_path = "FLECCS.xlsx"
    n_clusters = [-5, 1.7, 10]
    with pytest.raises(
        ValueError,
        match=(f"n_clusters must be an integer, but {n_clusters[1]} is not an integer"),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(file_path, sheet='2035 - NREL', column_name='MiNg_$100_CAISO', n_clusters=n_clusters[1], horizon_length=24,)
    
    with pytest.raises(
        ValueError,
        match=(f"n_clusters must be > 0, but {n_clusters[0]} is provided."),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(file_path, sheet='2035 - NREL', column_name='MiNg_$100_CAISO', n_clusters=n_clusters[0], horizon_length=24,)
    
    file_path = 'trash'
    with pytest.raises(
        ValueError,
        match=(f"The file path {file_path} does not exist. Please check your file path."),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(file_path, sheet='2035 - NREL', column_name='MiNg_$100_CAISO', n_clusters=n_clusters[2], horizon_length=24,)

    file_path = "FLECCS_no_sheet.xlsx"
    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()
    
        m.append_lmp_data(file_path, column_name=1, n_clusters=n_clusters[2], horizon_length=24,)
    
        assert f"Excel file was provided but no sheet was specified. Using the first sheet of the excel file." in caplog.text

    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()
    
        m.append_lmp_data(file_path, sheet=0, n_clusters=n_clusters[2], horizon_length=24,)
    
        assert f"Data was provided but no column name was provided. Using the first column of the data." in caplog.text

@pytest.mark.unit
def test_cluster_lmp_data_logger_messages(excel_data):
    n_clusters = [-5, 1.7, 10]
    with pytest.raises(
        ValueError,
        match=(f"n_clusters must be an integer, but {n_clusters[1]} is not an integer"),
    ):
        m = PriceTakerModel()
        _, _ = m.cluster_lmp_data(excel_data, n_clusters[1])
    
    with pytest.raises(
        ValueError,
        match=(f"n_clusters must be > 0, but {n_clusters[0]} is provided."),
    ):
        m = PriceTakerModel()
        _, _ = m.cluster_lmp_data(excel_data, n_clusters[0])


@pytest.mark.unit
def test_deepgetattr_logger_messages(excel_data):
    class Data:
        def __init__(self, value):
            self.value = value
    
    A = Data(5)
    attr = 'valu'
    with pytest.raises(
        ValueError,
        match=(f"The attribute `{attr}` does not exist on the {A.__class__.__name__} instance passed. Make sure that `{attr}` is defined for that object."),
    ):
        deepgetattr(A, attr)


@pytest.mark.unit
def test_generate_daily_data_logger_messages(excel_data):
    raw_data = excel_data['BaseCaseTax']
    with pytest.raises(
        ValueError,
        match=(f"tried to generate daily data, but horizon length of {9000} exceeds raw_data length of {len(raw_data)}"),
    ):
        m = PriceTakerModel(horizon_length=9000)
        daily_data = m.generate_daily_data(raw_data)


def dfc_design(m, params, capacity_range=(650, 900)):
    _dfc_capacity = params["dfc_capacity"]
    _ng_flow = params["ng_flow"]
    _capex = params["capex"]
    _fom_factor = params["fom_factor"]

    m.capacity = Var(
        within=NonNegativeReals,
        initialize=_dfc_capacity,
        bounds=(0, capacity_range[1]),
        doc="Capacity of the power plant [in MW]",
    )

    # Define a variable that informs whether the plant needs to built or not.
    m.build_unit = Var(
        within=Binary,
        doc="1: Plant is built, 0: Plant is not built",
    )

    # Bound the capcity of the plant in the specified range
    m.capacity_lb_con = Constraint(expr=m.capacity >= m.build_unit * capacity_range[0])
    m.capacity_ub_con = Constraint(expr=m.capacity <= m.build_unit * capacity_range[1])

    # Compute the natural gas flowrate required at maximum capacity.
    m.ng_flow = Expression(
        expr=m.capacity * (_ng_flow / _dfc_capacity),
        doc="Computes the natural flowrate required [in kg/s] at full load",
    )

    # Define an expression for CAPEX and FOM
    m.capex = Expression(
        expr=_capex * (m.capacity / _dfc_capacity),
        doc="CAPEX of the power cycle [in 1000$]",
    )
    m.fom = Expression(
        expr=_fom_factor * m.capex,
        doc="Fixed O&M cost [in 1000$/year]",
    )
