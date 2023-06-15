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

import pytest
import os
import pandas as pd

from pyomo.environ import (
    Var,
    Binary,
    NonNegativeReals,
    Constraint,
    Expression,
)

from idaes.apps.grid_integration import DesignModel, OperationModel

import idaes.logger as idaeslog

from idaes.apps.grid_integration.multiperiod.price_taker_model import PriceTakerModel

import matplotlib.pyplot as plt


@pytest.mark.unit
def test_daily_data_size():
    m = PriceTakerModel

    file_path = os.path.join("FLECCS.xlsx")
    file_name = "FLECCS.xlsx"
    excel_data = pd.read_excel(file_path, sheet_name="2030 - Princeton")

    # Generate price data for each hour of every day in the data
    daily_data = m.reconfigure_raw_data(excel_data)

    # Check that there is a row for each hour
    assert len(daily_data) == 24

    # Check that there is a column for each day
    assert all(len(row) == 365 for row in daily_data)


@pytest.mark.unit
def test_determine_optimal_num_clusters():
    m = PriceTakerModel

    file_name = "FLECCS.xlsx"
    file_path = os.path.join(os.getcwd(), file_name)
    excel_data = pd.read_excel(file_path, sheet_name="2030 - Princeton")

    daily_data = m.reconfigure_raw_data(excel_data)
    n_clusters, inertia_values = m.get_optimal_n_clusters(daily_data)

    assert n_clusters == 14


@pytest.mark.unit
def test_elbow_plot():
    m = PriceTakerModel

    file_name = "FLECCS.xlsx"
    file_path = os.path.join(os.getcwd(), file_name)
    excel_data = pd.read_excel(file_path, sheet_name="2030 - Princeton")

    daily_data = m.reconfigure_raw_data(excel_data)
    m.get_elbow_plot(daily_data)

    assert plt.gcf() is not None


@pytest.mark.unit
def test_logger_messages(caplog, kmin=None, kmax=None):
    file_name = "FLECCS.xlsx"
    file_path = os.path.join(os.getcwd(), file_name)
    excel_data = pd.read_excel(file_path, sheet_name="2030 - Princeton")

    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel

        daily_data = m.reconfigure_raw_data(excel_data)
        m.get_optimal_n_clusters(daily_data, kmin=kmin, kmax=kmax)

        assert f"{kmax} was not set - using a default value of 14." in caplog.text

    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel
        kmin = 1
        kmax = 10

        daily_data = m.reconfigure_raw_data(excel_data)
        sample_weight = ["1", "2", "3"]
        m.get_optimal_n_clusters(
            daily_data, kmin=kmin, kmax=kmax, sample_weight=sample_weight
        )

        assert (
            f"Ensure that the dimensions of the datasets match or set sample_weight to None."
            in caplog.text
        )

    caplog.clear()
    with caplog.at_level(idaeslog.WARNING):
        m = PriceTakerModel
        kmin = 10
        kmax = 1

        daily_data = m.reconfigure_raw_data(excel_data)
        m.get_optimal_n_clusters(daily_data, kmin=kmin, kmax=kmax)

        assert f"kmin:{kmin} needs to be less than kmax:{kmax}." in caplog.text


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
