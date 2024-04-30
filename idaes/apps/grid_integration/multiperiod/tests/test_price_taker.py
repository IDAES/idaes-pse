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

import pyomo.environ as aml
from pyomo.common.config import ConfigValue, In

from idaes.core import FlowsheetBlock

from idaes.core.base.process_base import declare_process_block_class
from idaes.models.unit_models import SkeletonUnitModelData

from idaes.apps.grid_integration import DesignModel, OperationModel

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
def test_seed_value():
    m = PriceTakerModel()

    m.seed = 50

    assert m.seed == 50


@pytest.mark.unit
def test_daily_data_size(excel_data):
    m = PriceTakerModel()

    # Generate price data for each hour of every day in the data
    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])

    # Check that there is a row for each horizon length in a representative day
    assert len(daily_data) == m.horizon_length


@pytest.mark.unit
def test_determine_optimal_num_clusters(excel_data):
    # Added a range for optimal cluster values based on how the
    # plot appears visually. Test can be removed in the future if
    # failure occurs. This may depend on scikit-learn and kneed and
    # the interaction thereof.

    # Older versions get n_clusters = 10, Newer versions n_clusters = 11
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data)

    assert 9 <= n_clusters <= 15


@pytest.mark.unit
def test_elbow_plot(excel_data):
    m = PriceTakerModel()

    daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
    m.get_optimal_n_clusters(daily_data, plot=True)
    test_fig = plt.gcf()
    plt.close("all")

    assert test_fig is not None


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
def test_init_logger_messages(excel_data, caplog):
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
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
    build_bin_var = "build"
    up_time = [10, -5, 2.2]
    down_time = [10, -5, 2.2]
    with pytest.raises(
        ValueError,
        match=(f"up_time must be an integer, but {up_time[2]} is not an integer"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[2], down_time[0])

    with pytest.raises(
        ValueError,
        match=(f"down_time must be an integer, but {down_time[2]} is not an integer"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[0], down_time[2])

    with pytest.raises(
        ValueError,
        match=(f"up_time must be >= 1, but {up_time[1]} is not"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[1], down_time[0])

    with pytest.raises(
        ValueError,
        match=(f"down_time must be >= 1, but {down_time[1]} is not"),
    ):
        m = PriceTakerModel()
        m.add_startup_shutdown(des, oper, build_bin_var, up_time[0], down_time[1])

    # Test Not Implemented Error (Rep. Days used for su/sd code)
    with pytest.raises(
        NotImplementedError,
        match=(
            f"You tried to use representative days with minimum up or minimum downtime constraints. This is not yet supported."
        ),
    ):
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
            column_name="BaseCaseTax",
            n_clusters=5,
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
        m.add_startup_shutdown(
            op_blk="fs.sofc_operation",
            design_blk="sofc_design",
            build_binary_var="build_unit",
            up_time=4,
            down_time=4,
        )


@pytest.mark.unit
def test_optimal_clusters_logger_messages(excel_data):
    kmin = [-5, 10.2, 9]
    kmax = [-5, 10.2, 8]
    with pytest.raises(
        ValueError,
        match=(f"kmin must be an integer, but {kmin[1]} is not an integer"),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[1], kmax=kmax[2]
        )

    with pytest.raises(
        ValueError,
        match=(f"kmax must be an integer, but {kmax[1]} is not an integer"),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[2], kmax=kmax[1]
        )

    with pytest.raises(
        ValueError,
        match=(f"kmin must be > 0, but {kmin[0]} is provided."),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[0], kmax=kmax[2]
        )

    with pytest.raises(
        ValueError,
        match=(f"kmax must be > 0, but {kmax[0]} is provided."),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[2], kmax=kmax[0]
        )
        m.cluster_lmp_data(excel_data, n_clusters)

    with pytest.raises(
        ValueError,
        match=(f"kmin must be less than kmax, but {kmin[2]} >= {kmax[2]}"),
    ):
        m = PriceTakerModel()

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin[2], kmax=kmax[2]
        )

    with pytest.raises(
        ValueError,
        match=(
            f"Could not find elbow point for given kmin, kmax. Consider increasing the range of kmin, kmax."
        ),
    ):
        m = PriceTakerModel()

        kmin = 9
        kmax = 10

        daily_data = m.generate_daily_data(excel_data["BaseCaseTax"])
        n_clusters, inertia_values = m.get_optimal_n_clusters(
            daily_data, kmin=kmin, kmax=kmax
        )


# The following test doesn't pass on all systems, so the warning for n_clusters being close
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
def test_ramping_constraint_logger_messages(excel_data):
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
    with pytest.raises(
        ValueError,
        match=(
            f"startup_rate fraction must be between 0 and 1, but {su_rate[0]} is not."
        ),
    ):
        m = PriceTakerModel()
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

    with pytest.raises(
        ValueError,
        match=(
            f"shutdown_rate fraction must be between 0 and 1, but {sd_rate[0]} is not."
        ),
    ):
        m = PriceTakerModel()
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

    with pytest.raises(
        ValueError,
        match=(
            f"ramp_up_rate fraction must be between 0 and 1, but {op_ru_rate[0]} is not."
        ),
    ):
        m = PriceTakerModel()
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

    with pytest.raises(
        ValueError,
        match=(
            f"ramp_down_rate fraction must be between 0 and 1, but {op_rd_rate[0]} is not."
        ),
    ):
        m = PriceTakerModel()
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

    with pytest.raises(
        ValueError,
        match=(
            f"op_range_lb fraction must be between 0 and 1, but {op_range_lb[0]} is not."
        ),
    ):
        m = PriceTakerModel()
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

    with pytest.raises(
        ValueError,
        match=(
            f"op_range_lb fraction must be <= shut_down_rate, otherwise the system cannot reach the off state."
        ),
    ):
        m = PriceTakerModel()
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

    # Test NotImplementedError (linearization = True)
    with pytest.raises(
        NotImplementedError,
        match=(
            f"You tried use nonlinear capacity with linearization. This is not yet supported."
        ),
    ):
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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

    # Test Value Error (wrong constraint type)
    with pytest.raises(
        ValueError,
        match=(
            f"constraint_type must be either linear, or nonliner, but garbage is not."
        ),
    ):
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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
def test_add_capacity_limits_logger_messages(excel_data, caplog):
    # Test Value Error (wrong constraint type)
    with pytest.raises(
        ValueError,
        match=(
            f"constraint_type must be either linear, or nonliner, but garbage is not."
        ),
    ):
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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
            constraint_type="garbage",
            linearization=False,
        )

    # Test Not Implemented Error (linearization is True when nonlinear is chosen)
    with pytest.raises(
        NotImplementedError,
        match=(
            f"You tried use nonlinear capacity with linearization. This is not yet supported."
        ),
    ):
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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
            linearization=True,
        )


@pytest.mark.unit
def test_append_lmp_data_logger_messages(excel_data, caplog):
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS.xlsx"
    n_clusters = [-5, 1.7, 10]
    with pytest.raises(
        ValueError,
        match=(f"n_clusters must be an integer, but {n_clusters[1]} is not an integer"),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(
            file_path,
            sheet="2035 - NREL",
            column_name="MiNg_$100_CAISO",
            n_clusters=n_clusters[1],
            horizon_length=24,
        )

    with pytest.raises(
        ValueError,
        match=(f"n_clusters must be > 0, but {n_clusters[0]} is provided."),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(
            file_path,
            sheet="2035 - NREL",
            column_name="MiNg_$100_CAISO",
            n_clusters=n_clusters[0],
            horizon_length=24,
        )

    file_path = "trash"
    with pytest.raises(
        ValueError,
        match=(
            f"The file path {file_path} does not exist. Please check your file path."
        ),
    ):
        m = PriceTakerModel()
        m.append_lmp_data(
            file_path,
            sheet="2035 - NREL",
            column_name="MiNg_$100_CAISO",
            n_clusters=n_clusters[2],
            horizon_length=24,
        )

    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_no_sheet.xlsx"
    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()

        m.append_lmp_data(
            file_path,
            column_name=1,
            n_clusters=n_clusters[2],
            horizon_length=24,
        )

        assert (
            f"Excel file was provided but no sheet was specified. Using the first sheet of the excel file."
            in caplog.text
        )

    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel()

        m.append_lmp_data(
            file_path,
            sheet=0,
            n_clusters=n_clusters[2],
            horizon_length=24,
        )

        assert (
            f"Data was provided but no column name was provided. Using the first column of the data."
            in caplog.text
        )


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
def test_generate_daily_data_logger_messages():
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    raw_data = pd.read_csv(file_path)
    with pytest.raises(
        ValueError,
        match=(
            f"tried to generate daily data, but horizon length of {9000} exceeds raw_data length of {len(raw_data['BaseCaseTax'])}"
        ),
    ):
        m = PriceTakerModel(horizon_length=9000)
        m.append_lmp_data(
            file_path,
            column_name="BaseCaseTax",
            n_clusters=10,
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_message_no_op_blks(excel_data, caplog):
    # Tests building the model with startup/shutdown then ramping rate with LMP as a single year with all time points
    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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

        m.build_hourly_cashflows(
            revenue_streams=None,
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="NPV",
        )

        assert (
            f"build_hourly_cashflows was called but no operation blocks were found so hourly cashflow of the model was set to 0. If you have hourly costs, please manually assign them."
            in caplog.text
        )


@pytest.mark.unit
def test_build_multiperiod_model_no_LMP_logger_message(excel_data):
    # Tests building the model with startup/shutdown then ramping rate with LMP as a single year with all time points
    with pytest.raises(
        ValueError,
        match=(
            f"OperationModelData has been defined to automatically "
            + f"populate LMP data. However, m.LMP does not exist. "
            + f"Please run the append_lmp_data function first or set the "
            + f"append_lmp_data attribute to False when configuring "
            + f"your OperationModelData object."
        ),
    ):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        m._n_time_points = 240
        m.set_days = None
        m.set_years = None

        # m.sofc_design = DesignModel(model_func=SOFC_design_model, model_args={"min_power": 200, "max_power": 650})
        m.sofc_design = aml.Block()
        m.sofc_design.PMAX = 650
        m.sofc_design.PMIN = 200
        m.sofc_design.build_unit = 1
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
    with caplog.at_level(idaeslog.WARNING):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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

        m.build_hourly_cashflows(
            revenue_streams=None,
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="NPV",
        )

        assert (
            f"build_cashflows was called, but no design blocks were found so capex and FOM are 0. Please manually add your cost objective if you require one."
            in caplog.text
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_1(excel_data, caplog):
    # Tests building the model with startup/shutdown then ramping rate with LMP as a single year with all time points
    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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

        m.build_hourly_cashflows(
            revenue_streams=None,
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        m.build_cashflows(
            objective="NPV",
        )

        assert (
            f"No revenues were provided while building the hourly cashflow. Revenues will be set to 0."
            in caplog.text
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_2(excel_data, caplog):
    # Tests building the model with ramping rate then startup/shutdown with LMP as a single year with representative days
    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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

        m.build_hourly_cashflows(
            revenue_streams=["elec_revenue", "H2_revenue"],
            costs=None,
        )

        m.build_cashflows(
            objective="Annualized NPV",
        )

        assert (
            f"No costs were provided while building the hourly cashflow. Costs will be set to 0."
            in caplog.text
        )


@pytest.mark.unit
def test_build_hourly_cashflow_logger_messages_and_build_3(excel_data, caplog):
    # Tests building the model with ramping rate then startup/shutdown with LMP as a single year with representative days
    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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
    with caplog.at_level(idaeslog.WARNING):
        # Create an instance of the Pricetrackermodel class
        m = PriceTakerModel()

        # Appending the data to the model
        DATA_DIR = Path(__file__).parent
        file_path = DATA_DIR / "FLECCS_shortened.xlsx"
        m.append_lmp_data(
            file_path=file_path,
            sheet="2030 - Princeton",
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

        m.build_hourly_cashflows(
            revenue_streams=["elec_revenue", "H2_revenue"],
            costs=["fuel_cost", "hourly_fixed_costs"],
        )

        bad_obj = "Garbage"
        m.build_cashflows(
            objective=bad_obj,
        )

        assert (
            f"build_cashflows was called, but the objective type provided, {bad_obj}, is invalid. The objective has been set to 0. Please manually add your cost objective if you require one."
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
