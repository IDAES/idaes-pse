from prescient.simulator import Prescient

options = {'data_directory': '/home/jhjalvi/git/Prescient/downloads/rts_gmlc/deterministic_with_network_scenarios',\
           'simulate_out_of_sample': True,\
           'run_sced_with_persistent_forecast_errors': True,\
           'output_directory': 'bidding_plugin_test_multiperiod_rankine',\
           'start_date':'07-11-2020',\
           'num_days': 1,\
           'sced_horizon':4,\
           'ruc_horizon':48,\
           'compute_market_settlements': True,\
           'day_ahead_pricing': 'LMP',\
           'ruc_mipgap':0.01,\
           'symbolic_solver_labels': True,\
           'reserve_factor':0.0,\
           'deterministic_ruc_solver':'cbc',\
           'sced_solver':'cbc',\
           'plugin': {'doubleloop': {'module': 'plugin_rankine_double_loop.py',\
                                      'bidding_generator': '102_STEAM_3'}}\
            }

Prescient().simulate(**options)

