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

import importlib
from numbers import Number
from pathlib import Path
from typing import Dict, Union, List

import pytest


prescient_simulator = pytest.importorskip("prescient.simulator", reason="prescient (optional dependency) not available")


@pytest.fixture(scope="module")
def base_dir() -> Path:
    pkg_init_path = Path(importlib.util.find_spec("idaes.tests.prescient").origin)
    return pkg_init_path.parent


# define custom type for type hinting
PrescientOptions = Dict[str, Union[str, bool, Number, dict]]


class Test5Bus:
    "Simple test using 5bus use case"

    @pytest.fixture
    def data_path(self, base_dir: Path) -> Path:
        return base_dir / "5bus"

    @pytest.mark.unit
    def test_data_path_available(self, data_path: Path):
        assert data_path.is_dir()

    @pytest.fixture
    def output_dir(self, tmp_path: Path) -> Path:
        path = tmp_path / "5bus_output"
        path.mkdir()
        return path

    @pytest.fixture
    def prescient_options(self, data_path: Path, output_dir: Path) -> PrescientOptions:
        return {
            "data_path": str(data_path),
            "input_format":"rts-gmlc",
            "simulate_out_of_sample":True,
            "run_sced_with_persistent_forecast_errors":True,
            "output_directory": str(output_dir),
            "start_date":"07-10-2020",
            "num_days":7,
            "sced_horizon":1,
            "ruc_mipgap":0.01,
            "reserve_factor":0.1,
            "deterministic_ruc_solver":"cbc",
            "deterministic_ruc_solver_options":{"feas":"off", "DivingF":"on",},
            "sced_solver":"cbc",
            "sced_frequency_minutes":5,
            "ruc_horizon":36,
            "compute_market_settlements":True,
            "monitor_all_contingencies":True,
            "output_solver_logs":False,
            "price_threshold":1000,
            "contingency_price_threshold":100,
            "reserve_price_threshold":5,
        }

    @pytest.fixture
    def run_simulator(self, prescient_options: PrescientOptions) -> None:
        from prescient.simulator import Prescient

        sim = Prescient()
        sim.simulate(**prescient_options)

    @pytest.fixture
    def simulation_output_files(self, run_simulator, output_dir: Path) -> List[Path]:
        return list(sorted(output_dir.rglob("*")))

    # FIXME this should be changed to pytest.mark.integration, unless we manage to
    # set the options so that in runs in < 20 s
    @pytest.mark.component
    # TODO use a more specific test to validate simulation output
    def test_simulation_results(self, simulation_output_files: List[Path], expected_number=2):
        assert len(simulation_output_files) >= expected_number
