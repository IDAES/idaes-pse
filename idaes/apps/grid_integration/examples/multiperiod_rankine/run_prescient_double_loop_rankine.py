from prescient.simulator import Prescient
import os

# NOTE: `rts_gmlc_data_dir` should point to a directory containing RTS-GMLC scenarios
rts_gmlc_data_dir = "/home/xgao1/DowlingLab/RTS-GMLC/RTS_Data/SourceData"

options = {
    "data_path": rts_gmlc_data_dir,
    "input_format": "rts-gmlc",
    "simulate_out_of_sample": True,
    "run_sced_with_persistent_forecast_errors": True,
    "output_directory": "bidding_plugin_test_multiperiod_rankine",
    "start_date": "07-11-2020",
    "num_days": 1,
    "sced_horizon": 4,
    "ruc_horizon": 48,
    "compute_market_settlements": True,
    "day_ahead_pricing": "LMP",
    "ruc_mipgap": 0.05,
    "symbolic_solver_labels": True,
    "reserve_factor": 0.0,
    "deterministic_ruc_solver": "gurobi",
    "sced_solver": "gurobi",
    "plugin": {
        "doubleloop": {
            "module": "plugin_double_loop_rankine.py",
            "bidding_generator": "102_STEAM_3",
        }
    },
}

Prescient().simulate(**options)
