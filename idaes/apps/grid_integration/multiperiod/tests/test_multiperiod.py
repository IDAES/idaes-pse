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

"""Tests for `MutiPeriodModel` class """

__author__ = "Radhakrishna Tumbalam Gooty"

import pytest
import matplotlib.pyplot as plt
import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.core.util.model_statistics import degrees_of_freedom


def build_flowsheet(m=None):
    """This function builds a dummy flowsheet"""
    if m is None:
        m = pyo.ConcreteModel()

    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.x = pyo.Var(within=pyo.NonNegativeReals)
    m.fs.y = pyo.Var(within=pyo.NonNegativeReals)
    m.fs.con1 = pyo.Constraint(expr=m.fs.x + m.fs.y == 1)

    return m


def fix_dof_and_initialize(m):
    """This function fixes dof and initializes the dummy flowsheet"""
    m.fs.y.fix(0.5)
    m.fs.x.set_value(0.5)


def fix_dof_and_initialize_2(m):
    """This function fixes dof and initializes the dummy flowsheet"""
    m.fs.con2 = pyo.Constraint(expr=m.fs.x + m.fs.y == 2)


def unfix_dof(m):
    """This function unfixes dof for optimization"""
    m.fs.y.unfix()


def get_linking_variable_pairs(m1, m2):
    """This function returns pairs of linking variables"""
    return [(m1.fs.y, m2.fs.y)]


@pytest.fixture(scope="module")
def build_multi_period_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_period_stochastic_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_scenarios=[1, 2, 3],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_day_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_days=["d1", "d2", "d3", "d4"],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_day_stochastic_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_days=["d1", "d2", "d3", "d4"],
        set_scenarios=[1, 2, 3],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_year_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_years=["y1", "y2", "y3"],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_year_stochastic_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_years=["y1", "y2", "y3"],
        set_scenarios=[1, 2, 3],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_day_year_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_days=["d1", "d2", "d3", "d4"],
        set_years=["y1", "y2", "y3"],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_day_year_stochastic_model():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        set_days=["d1", "d2", "d3", "d4"],
        set_years=["y1", "y2", "y3"],
        set_scenarios=[1, 2, 3],
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )

    return m


@pytest.fixture(scope="module")
def build_multi_period_model_manual():
    m = MultiPeriodModel(
        n_time_points=5,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )
    m.build_multi_period_model()
    return m


@pytest.mark.unit
def test_lmp_plot(monkeypatch):
    lmp = [2, 3, 5, 6, 8, 6, 0, 4]
    lmp2 = {"Day1": lmp, "Day2": lmp}
    lmp3 = {
        "Day1": lmp,
        "Day2": lmp,
        "Day3": lmp,
        "Day4": lmp,
        "Day5": lmp,
        "Day6": lmp,
        "Day7": lmp,
    }
    time_vec = [i for i in range(len(lmp))]
    time_vec2 = {"Day1": time_vec, "Day2": time_vec}

    monkeypatch.setattr(plt, "show", lambda: None)

    MultiPeriodModel.plot_lmp_signal(lmp)
    MultiPeriodModel.plot_lmp_signal(lmp, time_vec, x_range=(0, 8), y_range=(0, 10))
    MultiPeriodModel.plot_lmp_signal(lmp2)
    MultiPeriodModel.plot_lmp_signal(lmp2, time_vec2)

    with pytest.raises(
        Exception,
        match=(
            f"Number of LMP signals provided exceeds six: the maximum "
            f"number of subplots the function supports."
        ),
    ):
        MultiPeriodModel.plot_lmp_signal(lmp3)


@pytest.mark.unit
def test_lmp_and_schedule_plot(monkeypatch):
    lmp = [2, 3, 5, 6, 8, 6, 0, 4]
    time_vec = [i for i in range(len(lmp))]

    schedule1 = {
        "var1": [2, 3, 4, 6, 5, 8, 2, 7],
        "var2": [2, 3, 4, 6, 5, 8, 2, 7],
    }

    schedule2 = {
        "var1": [2, 3, 4, 6, 5, 8, 2, 7],
        "var2": [2, 3, 4, 6, 5, 8, 2, 7],
        "var3": [2, 3, 4, 6, 5, 8, 2, 7],
        "var4": [2, 3, 4, 6, 5, 8, 2, 7],
        "var5": [2, 3, 4, 6, 5, 8, 2, 7],
    }

    monkeypatch.setattr(plt, "show", lambda: None)

    with pytest.raises(
        Exception,
        match=(
            f"Number of elements in schedule exceeds four: "
            f"the maximum number of subplots the function supports."
        ),
    ):
        MultiPeriodModel.plot_lmp_and_schedule(lmp, schedule2)

    MultiPeriodModel.plot_lmp_and_schedule(
        lmp, schedule1, color={1: "tab:blue", 2: "magenta"}
    )
    MultiPeriodModel.plot_lmp_and_schedule(
        lmp,
        schedule1,
        time_vec,
        x_range=(0, 8),
        lmp_range=(0, 9),
        y_label={"var1": "Power [MW]"},
        y_range={"var2": (0, 9)},
    )


@pytest.mark.unit
def test_multi_period_model(build_multi_period_model):
    m = build_multi_period_model

    assert hasattr(m, "period")
    assert len(m.period) == 5
    assert hasattr(m, "link_constraints")
    for i in [1, 2, 3, 4]:
        assert hasattr(m.link_constraints[i], "link_constraints")

    assert degrees_of_freedom(m) == 1


@pytest.mark.unit
def test_multi_period_stochastic_model(build_multi_period_stochastic_model):
    m = build_multi_period_stochastic_model

    assert not hasattr(m, "period")
    assert not hasattr(m, "link_constraints")
    assert hasattr(m, "scenario")
    assert len(m.scenario) == 3

    for i in [1, 2, 3]:
        assert hasattr(m.scenario[i], "period")
        assert hasattr(m.scenario[i], "link_constraints")
        assert len(m.scenario[i].period) == 5

        for j in [1, 2, 3, 4]:
            assert hasattr(m.scenario[i].link_constraints[j], "link_constraints")

    assert degrees_of_freedom(m) == 3  # num_scenarios


@pytest.mark.unit
def test_multi_day_model(build_multi_day_model):
    m = build_multi_day_model

    assert hasattr(m, "period")
    assert hasattr(m, "link_constraints")
    for d in ["d1", "d2", "d3", "d4"]:
        for t in [1, 2, 3, 4, 5]:
            assert (t, d) in m.period

            if t != 5:
                assert (t, d) in m.link_constraints
                assert hasattr(m.link_constraints[t, d], "link_constraints")

    assert degrees_of_freedom(m) == 4  # num_days


@pytest.mark.unit
def test_multi_day_stochastic_model(build_multi_day_stochastic_model):
    m = build_multi_day_stochastic_model

    assert not hasattr(m, "period")
    assert not hasattr(m, "link_constraints")
    assert hasattr(m, "scenario")

    for s in m.scenario:
        for d in ["d1", "d2", "d3", "d4"]:
            for t in [1, 2, 3, 4, 5]:
                assert (t, d) in m.scenario[s].period

                if t != 5:
                    assert (t, d) in m.scenario[s].link_constraints
                    assert hasattr(
                        m.scenario[s].link_constraints[t, d], "link_constraints"
                    )

    assert degrees_of_freedom(m) == 12  # num_days * num_scenarios


@pytest.mark.unit
def test_multi_year_model(build_multi_year_model):
    m = build_multi_year_model

    assert hasattr(m, "period")
    assert hasattr(m, "link_constraints")
    for y in ["y1", "y2", "y3"]:
        for t in [1, 2, 3, 4, 5]:
            assert (t, y) in m.period

            if t != 5:
                assert (t, y) in m.link_constraints
                assert hasattr(m.link_constraints[t, y], "link_constraints")

    assert degrees_of_freedom(m) == 3  # num_years


@pytest.mark.unit
def test_multi_year_stochastic_model(build_multi_year_stochastic_model):
    m = build_multi_year_stochastic_model

    assert not hasattr(m, "period")
    assert not hasattr(m, "link_constraints")
    assert hasattr(m, "scenario")

    for s in m.scenario:
        for y in ["y1", "y2", "y3"]:
            for t in [1, 2, 3, 4, 5]:
                assert (t, y) in m.scenario[s].period

                if t != 5:
                    assert (t, y) in m.scenario[s].link_constraints
                    assert hasattr(
                        m.scenario[s].link_constraints[t, y], "link_constraints"
                    )

    assert degrees_of_freedom(m) == 9  # num_years * num_scenarios


@pytest.mark.unit
def test_multi_day_year_model(build_multi_day_year_model):
    m = build_multi_day_year_model

    assert hasattr(m, "period")
    assert hasattr(m, "link_constraints")
    for y in ["y1", "y2", "y3"]:
        for d in ["d1", "d2", "d3", "d4"]:
            for t in [1, 2, 3, 4, 5]:
                assert (t, d, y) in m.period

                if t != 5:
                    assert (t, d, y) in m.link_constraints
                    assert hasattr(m.link_constraints[t, d, y], "link_constraints")

    assert degrees_of_freedom(m) == 12  # num_days * num_years


@pytest.mark.unit
def test_multi_day_year_stochastic_model(build_multi_day_year_stochastic_model):
    m = build_multi_day_year_stochastic_model

    assert not hasattr(m, "period")
    assert not hasattr(m, "link_constraints")
    assert hasattr(m, "scenario")

    for s in m.scenario:
        for y in ["y1", "y2", "y3"]:
            for d in ["d1", "d2", "d3", "d4"]:
                for t in [1, 2, 3, 4, 5]:
                    assert (t, d, y) in m.scenario[s].period

                    if t != 5:
                        assert (t, d, y) in m.scenario[s].link_constraints
                        assert hasattr(
                            m.scenario[s].link_constraints[t, d, y], "link_constraints"
                        )

    assert degrees_of_freedom(m) == 36  # num_days * num_years * num_scenarios


@pytest.mark.unit
def test_no_initialization():
    # Cover warning associated with unfix_dof not provided
    m = MultiPeriodModel(
        n_time_points=2,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        initialization_func=fix_dof_and_initialize,
    )

    # Cover warning associated with initialization_func not provided
    m = MultiPeriodModel(
        n_time_points=2,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
    )

    # Cover warning associated with linking_variable_func not provided
    m = MultiPeriodModel(
        n_time_points=2,
        process_model_func=build_flowsheet,
        linking_variable_func=None,
        use_stochastic_build=True,
    )

    # Cover warning associated with periodic_variable_func provided
    m = MultiPeriodModel(
        n_time_points=2,
        process_model_func=build_flowsheet,
        linking_variable_func=get_linking_variable_pairs,
        periodic_variable_func=get_linking_variable_pairs,
        use_stochastic_build=True,
        initialization_func=fix_dof_and_initialize,
        unfix_dof_func=unfix_dof,
    )


@pytest.mark.unit
def test_initialization_fail():
    with pytest.raises(
        Exception,
        match=(
            f"Flowsheet did not converge to optimality after fixing the degrees of freedom. "
            f"To create the multi-period model without initialization, do not provide "
            f"initialization_func argument."
        ),
    ):
        m = MultiPeriodModel(
            n_time_points=5,
            process_model_func=build_flowsheet,
            linking_variable_func=get_linking_variable_pairs,
            use_stochastic_build=True,
            initialization_func=fix_dof_and_initialize_2,
            unfix_dof_func=unfix_dof,
        )


@pytest.mark.unit
def test_multi_period_model_manual(build_multi_period_model_manual):
    m = build_multi_period_model_manual

    assert hasattr(m, "blocks")
    assert len(m.blocks) == 5

    assert degrees_of_freedom(m) == 1
