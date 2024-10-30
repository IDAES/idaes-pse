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
import matplotlib.pyplot as plt
import pandas as pd

from pyomo.environ import (
    Var,
    Binary,
    NonNegativeReals,
)

import pyomo.environ as aml
from pyomo.common.config import ConfigValue

from idaes.core import FlowsheetBlock

from idaes.core.base.process_base import declare_process_block_class
from idaes.models.unit_models import SkeletonUnitModelData

from idaes.apps.grid_integration import DesignModel, OperationModel

import idaes.logger as idaeslog

from idaes.apps.grid_integration.multiperiod.price_taker_model import PriceTakerModel
from idaes.core.util.exceptions import ConfigurationError


@pytest.fixture
def excel_data():
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    data = pd.read_csv(file_path)
    return data


@pytest.mark.unit
def test_seed():
    m = PriceTakerModel()
    assert isinstance(m._seed, int)


@pytest.mark.unit
def test_seed_value():
    value = "fifty"
    with pytest.raises(
        TypeError,
        match=("seed must be an integer, but fifty is not an integer"),
    ):

        m = PriceTakerModel()
        m.seed = value


@pytest.mark.unit
def test_daily_data_size(excel_data):
    m = PriceTakerModel()

    # Generate price data for each hour of every day in the data
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])

    # Check that there is a row for each horizon length in a representative day
    assert len(daily_data) == m.horizon_length


@pytest.mark.unit
def test_determine_optimal_num_clusters(excel_data):
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data, kmax=30)

    assert n_clusters == 11


@pytest.mark.unit
def test_generate_elbow_plot(excel_data):
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    m.generate_elbow_plot(daily_data)

    # Test that a figure was created
    assert plt.gcf() is not None
    # Test that axes were created
    assert plt.gca() is not None
    # Test that the plot has data
    assert plt.gca().has_data()

    plt.close("all")


@pytest.mark.unit
def test_cluster_lmp_data(excel_data):
    # This function gets within both if statement expression in the
    # cluster_lmp_data function currently.
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data)
    lmp_data, weights = m.cluster_lmp_data(
        excel_data["BaseCaseTax"], n_clusters=n_clusters
    )

    sum_of_weights = 0
    for i in range(1, n_clusters + 1):
        sum_of_weights = sum_of_weights + weights[0][i]
    assert sum_of_weights == 365

    assert len(lmp_data) == n_clusters


@pytest.mark.unit
def test_init_logger_messages_clusters(excel_data, caplog):
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        m.get_optimal_n_clusters(daily_data)

        assert "kmax was not set - using a default value of 30." in caplog.text


@pytest.mark.unit
def test_init_logger_message1(excel_data, caplog):
    # Testing horizon_length input value errors
    value = 0
    with pytest.raises(
        ValueError,
        match=("horizon_length must be > 0, but 0 is provided."),
    ):
        m = PriceTakerModel()
        m.horizon_length = value


@pytest.mark.unit
def test_init_logger_message2(excel_data, caplog):
    value = 12.34
    m = PriceTakerModel()
    with pytest.raises(
        TypeError,
        match=("horizon_length must be an integer, but 12.34 is not an integer"),
    ):
        m.horizon_length = value


@pytest.mark.unit
def test_init_logger_message3(excel_data, caplog):
    # Testing seed errors
    value = 12.34
    m = PriceTakerModel()
    with pytest.raises(
        TypeError,
        match=("seed must be an integer, but 12.34 is not an integer"),
    ):
        m.seed = value


@pytest.mark.unit
def test_min_up_time_logger_message1(excel_data):
    des = DesignModel()
    oper = OperationModel()
    build_bin_var = "build"
    up_time = [10, -5, 2.2]
    down_time = [10, -5, 2.2]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("up_time must be an integer, but 2.2 is not an integer"),
    ):
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[2], down_time[0])


@pytest.mark.unit
def test_min_up_time_logger_message2(excel_data):
    des = DesignModel()
    oper = OperationModel()
    build_bin_var = "build"
    up_time = [10, -5, 2.2]
    down_time = [10, -5, 2.2]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("up_time must be >= 1, but -5 is not"),
    ):
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[1], down_time[0])


@pytest.mark.unit
def test_min_down_time_logger_message1(excel_data):
    des = DesignModel()
    oper = OperationModel()
    build_bin_var = "build"
    up_time = [10, -5, 2.2]
    down_time = [10, -5, 2.2]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("down_time must be an integer, but 2.2 is not an integer"),
    ):
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[0], down_time[2])


@pytest.mark.unit
def test_min_down_time_logger_message2(excel_data):
    des = DesignModel()
    oper = OperationModel()
    build_bin_var = "build"
    up_time = [10, -5, 2.2]
    down_time = [10, -5, 2.2]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("down_time must be >= 1, but -5 is not"),
    ):
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[0], down_time[1])


@pytest.mark.unit
def test_init_logger_messages_clusters_min_up_down_time(excel_data, caplog):
    # Test Not Implemented Error (Rep. Days used for su/sd code)
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
        n_clusters=2,
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    with pytest.raises(
        NotImplementedError,
        match=(
            "You tried to use representative days with minimum up or minimum downtime constraints. This is not yet supported."
        ),
    ):
        # add capacity limit constraints
        m.add_startup_shutdown(
            op_blk="fs.sofc_operation",
            design_blk="sofc_design",
            build_binary_var="build_unit",
            up_time=4,
            down_time=4,
        )


@pytest.mark.unit
def test_optimal_clusters_kmin_logger_message1(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    m = PriceTakerModel()
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    with pytest.raises(
        ValueError,
        match=("kmin must be an integer, but 10.2 is not an integer"),
    ):
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[1], kmax=kmax[2]
        )


@pytest.mark.unit
def test_optimal_clusters_kmin_logger_message2(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    m = PriceTakerModel()
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    with pytest.raises(
        ValueError,
        match=("kmin must be > 0, but -5 is provided."),
    ):
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[0], kmax=kmax[2]
        )


@pytest.mark.unit
def test_optimal_clusters_kmin_logger_message3(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    m = PriceTakerModel()
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    with pytest.raises(
        ValueError,
        match=("kmin must be less than kmax, but 9 >= 8"),
    ):
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[2], kmax=kmax[2]
        )


@pytest.mark.unit
def test_optimal_clusters_kmax_logger_message1(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    m = PriceTakerModel()
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    with pytest.raises(
        ValueError,
        match=("kmax must be an integer, but 10.2 is not an integer"),
    ):
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[2], kmax=kmax[1]
        )


@pytest.mark.unit
def test_optimal_clusters_kmax_logger_message2(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    m = PriceTakerModel()
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    with pytest.raises(
        ValueError,
        match=("kmax must be > 0, but -5 is provided."),
    ):
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[2], kmax=kmax[0]
        )
        m.cluster_lmp_data(excel_data, n_clusters)


# @pytest.mark.unit
# def test_failed_imports(excel_data):
#     m = PriceTakerModel()
#     kmin = 9
#     kmax = 10
#     daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
#     with pytest.raises(
#         ImportError,
#         match=(
#             "Optimal cluster feature requires optional imports 'scikit-learn' and 'kneed'."
#         ),
#     ):
#         n_clusters, inertia_values = m.get_optimal_n_clusters(
#             daily_data, kmin=kmin, kmax=kmax
#         )


# The following test doesn't pass on all systems, so the warning for n_clusters being
# too close to kmax will be uncovered.
# @pytest.mark.unit
# def test_optimal_clusters_close_to_kmax(excel_data, caplog):
#     # Ideally the following test will work, however this will depend on the version of scikit-learn and
#     # kneed. It is possible that if these packages change, the test will fail. In that case, the test
#     # could be removed as the function works properly, but some lines of code will not be covered.
#     caplog.clear()
#     with caplog.at_level(idaeslog.WARNING):
#         m = PriceTakerModel()
#         kmin = 9
#         kmax = 14

#         daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
#         m.get_optimal_n_clusters(daily_data, kmin=kmin, kmax=kmax)

#         assert (
#             f"Optimal number of clusters is close to kmax: {kmax}. Consider increasing kmax."
#             in caplog.text
#         )


@pytest.mark.unit
def test_ramping_constraint_logger_message1(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("startup_rate fraction must be between 0 and 1, but -0.1 is not."),
    ):
        m.add_ramping_constraints(
            des,
            oper,
            capac_var,
            ramping_var,
            constraint_type,
            linearization,
            op_range_lb[1],
            su_rate[0],
            sd_rate[1],
            op_ru_rate[1],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message2(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("shutdown_rate fraction must be between 0 and 1, but -0.1 is not."),
    ):
        m.add_ramping_constraints(
            des,
            oper,
            capac_var,
            ramping_var,
            constraint_type,
            linearization,
            op_range_lb[1],
            su_rate[1],
            sd_rate[0],
            op_ru_rate[1],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message3(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("ramp_up_rate fraction must be between 0 and 1, but -0.1 is not."),
    ):
        m.add_ramping_constraints(
            des,
            oper,
            capac_var,
            ramping_var,
            constraint_type,
            linearization,
            op_range_lb[1],
            su_rate[1],
            sd_rate[1],
            op_ru_rate[0],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message4(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("ramp_down_rate fraction must be between 0 and 1, but -0.1 is not."),
    ):
        m.add_ramping_constraints(
            des,
            oper,
            capac_var,
            ramping_var,
            constraint_type,
            linearization,
            op_range_lb[1],
            su_rate[1],
            sd_rate[1],
            op_ru_rate[1],
            op_rd_rate[0],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message5(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("op_range_lb fraction must be between 0 and 1, but -0.1 is not."),
    ):
        m.add_ramping_constraints(
            des,
            oper,
            capac_var,
            ramping_var,
            constraint_type,
            linearization,
            op_range_lb[0],
            su_rate[1],
            sd_rate[1],
            op_ru_rate[1],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message6(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=(
            "op_range_lb fraction must be <= shut_down_rate, otherwise the system cannot reach the off state."
        ),
    ):
        m.add_ramping_constraints(
            des,
            oper,
            capac_var,
            ramping_var,
            constraint_type,
            linearization,
            op_range_lb[2],
            su_rate[1],
            sd_rate[1],
            op_ru_rate[1],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message1(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]

    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    # m.sofc_design = DesignModel(model_func=SOFC_design_model, model_args={"min_power": 200, "max_power": 650})
    m.sofc_design = aml.Block()
    m.sofc_design.PMAX = 650
    m.sofc_design.PMIN = 200
    m.sofc_design.build_unit = 1
    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # Test NotImplementedError (linearization = True)
    with pytest.raises(
        NotImplementedError,
        match=(
            "You tried use nonlinear capacity with linearization. This is not yet supported."
        ),
    ):

        m.add_ramping_constraints(
            "fs.sofc_operation",
            "sofc_design",
            capac_var,
            ramping_var,
            "nonlinear",
            True,
            op_range_lb[1],
            su_rate[1],
            sd_rate[1],
            op_ru_rate[1],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_ramping_constraint_logger_message1(excel_data):
    des = DesignModel()
    oper = OperationModel()
    ramping_var = "power"
    capac_var = "power_max"
    constraint_type = "linear"
    linearization = False
    op_range_lb = [-0.1, 0.5, 1.0]
    su_rate = [-0.1, 0.5]
    sd_rate = [-0.1, 0.5]
    op_ru_rate = [-0.1, 0.5]
    op_rd_rate = [-0.1, 0.5]

    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )
    # m.sofc_design = DesignModel(model_func=SOFC_design_model, model_args={"min_power": 200, "max_power": 650})
    m.sofc_design = aml.Block()
    m.sofc_design.PMAX = 650
    m.sofc_design.PMIN = 200
    m.sofc_design.build_unit = 1
    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # Test Value Error (wrong constraint type)
    with pytest.raises(
        ValueError,
        match=(
            "constraint_type must be either linear, or nonliner, but garbage is not."
        ),
    ):

        m.add_ramping_constraints(
            "fs.sofc_operation",
            "sofc_design",
            capac_var,
            ramping_var,
            "garbage",
            True,
            op_range_lb[1],
            su_rate[1],
            sd_rate[1],
            op_ru_rate[1],
            op_rd_rate[1],
        )


@pytest.mark.unit
def test_add_capacity_limits_logger_message1(excel_data, caplog):
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # Test Value Error (wrong constraint type)
    with pytest.raises(
        ValueError,
        match=(
            "constraint_type must be either linear, or nonliner, but garbage is not."
        ),
    ):
        # add capacity limit constraints
        m.add_capacity_limits(
            op_blk="fs.sofc_operation",
            design_blk="sofc_design",
            commodity_var="power",
            capacity_max="PMAX",
            capacity_min="PMIN",
            constraint_type="garbage",
            linearization=False,
        )


@pytest.mark.unit
def test_add_capacity_limits_logger_message2(excel_data, caplog):
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # Test Not Implemented Error (linearization is True when nonlinear is chosen)
    with pytest.raises(
        NotImplementedError,
        match=(
            "You tried use nonlinear capacity with linearization. This is not yet supported."
        ),
    ):
        # add capacity limit constraints
        m.add_capacity_limits(
            op_blk="fs.sofc_operation",
            design_blk="sofc_design",
            commodity_var="power",
            capacity_max="PMAX",
            capacity_min="PMIN",
            constraint_type="nonlinear",
            linearization=True,
        )


@pytest.mark.unit
def test_append_lmp_data_logger_message1(excel_data, caplog):
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    n_clusters = [-5, 1.7, 10]
    with pytest.raises(
        ValueError,
        match=("n_clusters must be an integer, but 1.7 is not an integer"),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(
            file_path,
            column_name="BaseCaseTax",
            n_clusters=n_clusters[1],
            horizon_length=24,
        )


@pytest.mark.unit
def test_append_lmp_data_logger_message1(excel_data, caplog):
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    n_clusters = [-5, 1.7, 10]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("n_clusters must be > 0, but -5 is provided."),
    ):
        m.append_lmp_data(
            file_path,
            column_name="BaseCaseTax",
            n_clusters=n_clusters[0],
            horizon_length=24,
        )


@pytest.mark.unit
def test_append_lmp_data_logger_message1(excel_data, caplog):
    file_path = "trash"
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("The file path trash does not exist. Please check your file path."),
    ):
        m.append_lmp_data(
            file_path,
        )


@pytest.mark.unit
def test_append_lmp_data_logger_message1(excel_data, caplog):
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=(
            "Data was provided but no column name was provided. Please supply a value for column_name."
        ),
    ):
        m.append_lmp_data(
            file_path,
        )


@pytest.mark.unit
def test_cluster_lmp_data_logger_message1(excel_data):
    n_clusters = [-5, 1.7, 10]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("n_clusters must be an integer, but 1.7 is not an integer"),
    ):
        _, _ = m.cluster_lmp_data(excel_data, n_clusters[1])


@pytest.mark.unit
def test_cluster_lmp_data_logger_message2(excel_data):
    n_clusters = [-5, 1.7, 10]
    m = PriceTakerModel()
    with pytest.raises(
        ValueError,
        match=("n_clusters must be > 0, but -5 is provided."),
    ):
        _, _ = m.cluster_lmp_data(excel_data, n_clusters[0])


@pytest.mark.unit
def test_generate_daily_data_logger_messages():
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    raw_data = pd.read_csv(file_path)
    m = PriceTakerModel(horizon_length=9000)
    with pytest.raises(
        ValueError,
        match=(
            f"tried to generate daily data, but horizon length of {9000} exceeds raw_data length of {len(raw_data['BaseCaseTax'])}"
        ),
    ):
        m.append_lmp_data(
            file_path,
            column_name="BaseCaseTax",
            n_clusters=10,
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_message_no_op_blks(excel_data, caplog):
    # Tests building the model with startup/shutdown then ramping rate with LMP as a single year with all time points
    caplog.clear()
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet_no_op_blk,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # ramping and startup constraints
    m.add_startup_shutdown(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        build_binary_var="build_unit",
        up_time=4,
        down_time=4,
    )

    m.add_ramping_constraints(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        capacity_var="PMAX",
        ramping_var="power",
        constraint_type="linear",
        linearization=True,
        op_range_lb=200 / 650,
        startup_rate=1.0,
        shutdown_rate=1.0,
        ramp_up_rate=1.0,
        ramp_down_rate=1.0,
    )

    with caplog.at_level(idaeslog.WARNING):
        m.build_hourly_cashflows(
            revenue_streams=None,
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="NPV",
        )

        assert (
            "build_hourly_cashflows was called but no operation blocks were found so hourly cashflow of the model was set to 0. If you have hourly costs, please manually assign them."
            in caplog.text
        )


@pytest.mark.unit
def test_build_multiperiod_model_no_LMP_logger_message():
    # Checks the exception raised if issue arises with attaching LMP data to OperationModel
    # Create an instance of the PriceTakerModel class
    m = PriceTakerModel()

    m._n_time_points = 240
    m.set_days = None
    m.set_years = None

    with pytest.raises(
        ConfigurationError,
        match=(
            "OperationModelData has been defined to automatically "
            + "populate LMP data. However, self.LMP does not exist, where self is an instance of PriceTakerModel. "
            + "Please run the append_lmp_data method from the PriceTakerModel class first or set the "
            + "declare_lmp_param configuration option to False when configuring "
            + "your OperationModelData object."
        ),
    ):
        # Build the multiperiod model
        m.build_multiperiod_model(
            process_model_func=build_sofc_flowsheet_no_LMP,
            linking_variable_func=None,
            flowsheet_options={"sofc_design": None},
        )

@pytest.mark.unit
def test_build_multiperiod_model_before_append_lmp_data_logger_messages():
    m = PriceTakerModel()

    # First exception arises if _n_time_points is not an attribute of PriceTakerMode
    with pytest.raises(
        ConfigurationError,
        match=(
            "MultiPeriodModel requires n_time_points as an argument. Before invoking the build_multiperiod_model method, call the append_lmp_data method on PriceTakerModel class first, which will assign the number of time points, n_time_points, to be used in the MultiPeriodModel."
        ),
    ):
        # Build the multiperiod model
        m.build_multiperiod_model(
            process_model_func=build_sofc_flowsheet_no_LMP,
            linking_variable_func=None,
            flowsheet_options={"sofc_design": None},
        )
    
    # Next exception that would arise relates to either set_years or set_days attributes not existing, indicating append_lmp_data method was not run first.
    m1 = PriceTakerModel()
    m1._n_time_points = 240
    m1.set_days = None

    with pytest.raises(
    ConfigurationError,
    match=(
        "Before invoking the build_multiperiod_model method, call the append_lmp_data method on PriceTakerModel class first, which will assign the number of time points, n_time_points, to be used in the MultiPeriodModel."
    ),
    ):
        # Build the multiperiod model
        m.build_multiperiod_model(
            process_model_func=build_sofc_flowsheet_no_LMP,
            linking_variable_func=None,
            flowsheet_options={"sofc_design": None},
        )

    m2 = PriceTakerModel()
    m2._n_time_points = 240
    m2.set_years = None

    with pytest.raises(
    ConfigurationError,
    match=(
        "Before invoking the build_multiperiod_model method, call the append_lmp_data method on PriceTakerModel class first, which will assign the number of time points, n_time_points, to be used in the MultiPeriodModel."
    ),
    ):
        # Build the multiperiod model
        m.build_multiperiod_model(
            process_model_func=build_sofc_flowsheet_no_LMP,
            linking_variable_func=None,
            flowsheet_options={"sofc_design": None},
        )

@pytest.mark.unit
def test_build_hourly_cashflow_logger_message_no_des_blks(excel_data, caplog):
    # Tests building the model with startup/shutdown then ramping rate with LMP as a single year with all time points
    caplog.clear()
    # Create an instance of the Pricetrackermodel class
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    # m.sofc_design = DesignModel(model_func=SOFC_design_model, model_args={"min_power": 200, "max_power": 650})
    m.sofc_design = aml.Block()
    m.sofc_design.PMAX = 650
    m.sofc_design.PMIN = 200
    m.sofc_design.build_unit = 1
    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # ramping and startup constraints
    m.add_startup_shutdown(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        build_binary_var="build_unit",
        up_time=4,
        down_time=4,
    )

    m.add_ramping_constraints(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        capacity_var="PMAX",
        ramping_var="power",
        constraint_type="linear",
        linearization=True,
        op_range_lb=200 / 650,
        startup_rate=1.0,
        shutdown_rate=1.0,
        ramp_up_rate=1.0,
        ramp_down_rate=1.0,
    )

    with caplog.at_level(idaeslog.WARNING):
        m.build_hourly_cashflows(
            revenue_streams=None,
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="NPV",
        )

        assert (
            "build_cashflows was called, but no design blocks were found so capex and FOM are 0. Please manually add your cost objective if you require one."
            in caplog.text
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_1(excel_data, caplog):
    # Tests building the model with startup/shutdown then ramping rate with LMP as a single year with all time points
    caplog.clear()
    # Create an instance of the Pricetrackermodel class
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # ramping and startup constraints
    m.add_startup_shutdown(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        build_binary_var="build_unit",
        up_time=4,
        down_time=4,
    )

    m.add_ramping_constraints(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        capacity_var="PMAX",
        ramping_var="power",
        constraint_type="linear",
        linearization=True,
        op_range_lb=200 / 650,
        startup_rate=1.0,
        shutdown_rate=1.0,
        ramp_up_rate=1.0,
        ramp_down_rate=1.0,
    )

    # add capacity limit constraints
    m.add_capacity_limits(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        commodity_var="power",
        capacity_max="PMAX",
        capacity_min="PMIN",
        constraint_type="nonlinear",
        linearization=False,
    )

    with caplog.at_level(idaeslog.WARNING):
        m.build_hourly_cashflows(
            revenue_streams=None,
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="NPV",
        )

        assert (
            "No revenues were provided while building the hourly cashflow. Revenues will be set to 0."
            in caplog.text
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_2(excel_data, caplog):
    # Tests building the model with ramping rate then startup/shutdown with LMP as a single year with representative days
    caplog.clear()
    # Create an instance of the Pricetrackermodel class
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # ramping and startup constraints
    m.add_ramping_constraints(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        capacity_var="PMAX",
        ramping_var="power",
        constraint_type="linear",
        linearization=True,
        op_range_lb=200 / 650,
        startup_rate=1.0,
        shutdown_rate=1.0,
        ramp_up_rate=1.0,
        ramp_down_rate=1.0,
    )

    m.add_startup_shutdown(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        build_binary_var="build_unit",
        up_time=4,
        down_time=4,
    )

    # add capacity limit constraints
    m.add_capacity_limits(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        commodity_var="power",
        capacity_max="PMAX",
        capacity_min="PMIN",
        constraint_type="nonlinear",
        linearization=False,
    )

    with caplog.at_level(idaeslog.WARNING):
        m.build_hourly_cashflows(
            revenue_streams=["elec_revenue", "H2_revenue"],
            costs=None,
        )

        m.build_cashflows(
            objective="Annualized NPV",
        )

        assert (
            "No costs were provided while building the hourly cashflow. Costs will be set to 0."
            in caplog.text
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_3(excel_data, caplog):
    # Tests building the model with ramping rate then startup/shutdown with LMP as a single year with representative days
    caplog.clear()
    # Create an instance of the Pricetrackermodel class
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # add capacity limit constraints
    m.add_capacity_limits(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        commodity_var="power",
        capacity_max="PMAX",
        capacity_min="PMIN",
        constraint_type="nonlinear",
        linearization=False,
    )

    # ramping and startup constraints
    m.add_ramping_constraints(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        capacity_var="PMAX",
        ramping_var="power",
        constraint_type="linear",
        linearization=True,
        op_range_lb=200 / 650,
        startup_rate=1.0,
        shutdown_rate=1.0,
        ramp_up_rate=1.0,
        ramp_down_rate=1.0,
    )

    m.add_startup_shutdown(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        build_binary_var="build_unit",
        up_time=4,
        down_time=4,
    )

    with caplog.at_level(idaeslog.WARNING):
        m.build_hourly_cashflows(
            revenue_streams=["elec_revenue", "H2_revenue"],
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="Net Profit",
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_4(excel_data, caplog):
    # Tests building the model with ramping rate then startup/shutdown with LMP as a single year with representative days
    caplog.clear()
    # Create an instance of the Pricetrackermodel class
    m = PriceTakerModel()

    # Appending the data to the model
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton_shortened.csv"
    m.append_lmp_data(
        file_path=file_path,
        column_name="BaseCaseTax",
    )

    m.sofc_design = DesignModel(
        model_func=SOFC_design_model,
        model_args={"min_power": 200, "max_power": 650},
    )

    # Build the multiperiod model
    m.build_multiperiod_model(
        process_model_func=build_sofc_flowsheet,
        linking_variable_func=None,
        flowsheet_options={"sofc_design": None},
    )

    # add capacity limit constraints
    m.add_capacity_limits(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        commodity_var="power",
        capacity_max="PMAX",
        capacity_min="PMIN",
        constraint_type="nonlinear",
        linearization=False,
    )

    # ramping and startup constraints
    m.add_ramping_constraints(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        capacity_var="PMAX",
        ramping_var="power",
        constraint_type="linear",
        linearization=True,
        op_range_lb=200 / 650,
        startup_rate=1.0,
        shutdown_rate=1.0,
        ramp_up_rate=1.0,
        ramp_down_rate=1.0,
    )

    m.add_startup_shutdown(
        op_blk="fs.sofc_operation",
        design_blk="sofc_design",
        build_binary_var="build_unit",
        up_time=4,
        down_time=4,
    )

    with caplog.at_level(idaeslog.WARNING):
        m.build_hourly_cashflows(
            revenue_streams=["elec_revenue", "H2_revenue"],
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        bad_obj = "Garbage"
        m.build_cashflows(
            objective=bad_obj,
        )

        assert (
            "build_cashflows was called, but the objective type provided, Garbage, is invalid. The objective has been set to 0. Please manually add your cost objective if you require one."
            in caplog.text
        )


# Model for testing builds with Linear capacity constraints
#############################################################
def SOFC_design_model(
    m,
    max_power=None,
    min_power=None,
):
    # Capacity parameters
    m.PMAX = aml.Param(initialize=max_power, mutable=False)
    m.PMIN = aml.Param(initialize=min_power, mutable=False)

    m.build_unit = aml.Param(initialize=1, mutable=False)

    # Capital investment cost ($)
    m.capex = 616441.20
    # Fixed operating and investment ($ / year)
    m.fom = 433882.80


def SOFC_operation_model(
    m,
    sofc_design_blk=None,
):
    # LMP value dummy
    # m.LMP = aml.Param(initialize=1, mutable=True)

    # Operation Variables
    m.power = aml.Var(domain=aml.NonNegativeReals)

    # cost expressions
    # Linear version of surrogate cost constraint (combined fuel_cost and non_fuel_vom)
    m.fuel_cost = aml.Expression(expr=(23.2934 * m.power + 49.21504))
    # Nonlinear version of surrogate cost constraint (combined fuel_cost and non_fuel_vom)
    # m.fuel_cost = aml.Expression(expr=(23.2934 * m.power + 0.287571e-3 * m.power ** 2 + 0.1155903e-5 * m.power ** 3 + 49.21504))
    m.elec_revenue = aml.Expression(expr=(m.power * m.LMP))
    m.H2_revenue = aml.Expression(expr=(0))
    m.non_fuel_vom = aml.Expression(expr=(0))
    m.carbon_price = aml.Expression(expr=(0))

    fixed_cap_hourly = 70.37 * 1e6 / 365 / 24
    fixed_op_hourly = 49.53 * 1e6 / 365 / 24

    m.hourly_fixed_cost = aml.Expression(expr=(fixed_cap_hourly + fixed_op_hourly))


def build_sofc_flowsheet(
    m,
    sofc_design,
):

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.sofc_operation = OperationModel(
        model_func=SOFC_operation_model,
        model_args={
            "sofc_design_blk": sofc_design,
        },
    )

    # Flowsheet level aml.Variables (none for this case)

    # Flowsheet level constraints (none for this case)


#############################################################
# End model for testing builds with Linear


# Extra functions for different edge-case warning messages
#############################################################
@declare_process_block_class("SOFC_op")
class SOFC_op_blk(SkeletonUnitModelData):
    CONFIG = SkeletonUnitModelData.CONFIG()
    CONFIG.declare(
        "model_func",
        ConfigValue(
            doc="Function that builds the design model",
        ),
    )
    CONFIG.declare(
        "model_args",
        ConfigValue(
            default={},
            doc="Dictionary containing arguments needed for model_func",
        ),
    )

    def build(self):
        super().build()

        # Declare variables
        self.power = Var(
            within=NonNegativeReals,
            doc="Total power generated [in MW]",
        )
        self.op_mode = Var(
            within=Binary,
            doc="1: In Operation, 0: Shutdown",
        )
        self.startup = Var(
            within=Binary,
            doc="1: Plant is startup at this hour, 0: Otherwise",
        )
        self.shutdown = Var(
            within=Binary,
            doc="1: Plant is shutdown at this hour, 0: Otherwise",
        )


def build_sofc_flowsheet_no_op_blk(
    m,
    sofc_design,
):

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.sofc_operation = SOFC_op(
        model_func=SOFC_operation_model,
        model_args={
            "sofc_design_blk": sofc_design,
        },
    )

    # Flowsheet level aml.Variables (none for this case)

    # Flowsheet level constraints (none for this case)


def SOFC_operation_model_no_LMP(
    m,
    sofc_design_blk=None,
):
    # LMP value dummy
    # m.LMP = aml.Param(initialize=1, mutable=True)

    # Operation Variables
    m.power = aml.Var(domain=aml.NonNegativeReals)

    # cost expressions
    # Linear version of surrogate cost constraint (combined fuel_cost and non_fuel_vom)
    m.fuel_cost = aml.Expression(expr=(23.2934 * m.power + 49.21504))
    # Nonlinear version of surrogate cost constraint (combined fuel_cost and non_fuel_vom)
    # m.fuel_cost = aml.Expression(expr=(23.2934 * m.power + 0.287571e-3 * m.power ** 2 + 0.1155903e-5 * m.power ** 3 + 49.21504))
    m.elec_revenue = aml.Expression(expr=(m.power * 1.0))
    m.H2_revenue = aml.Expression(expr=(0))
    m.non_fuel_vom = aml.Expression(expr=(0))
    m.carbon_price = aml.Expression(expr=(0))

    fixed_cap_hourly = 70.37 * 1e6 / 365 / 24
    fixed_op_hourly = 49.53 * 1e6 / 365 / 24

    m.hourly_fixed_cost = aml.Expression(expr=(fixed_cap_hourly + fixed_op_hourly))


def build_sofc_flowsheet_no_LMP(
    m,
    sofc_design,
):

    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.sofc_operation = OperationModel(
        model_func=SOFC_operation_model_no_LMP,
        model_args={
            "sofc_design_blk": sofc_design,
        },
    )

    # Flowsheet level aml.Variables (none for this case)

    # Flowsheet level constraints (none for this case)


# Test cases to be done (4 total):
# (LMP1, RR1, MUD1, OBJ1) x
# (LMP2, RR2, MUD2, OBJ2) x
# (LMP3, RR3, MUD2, OBJ3)
# (LMP4, MUD1, RR1, OBJ1)
