from optparse import Option
import prescient.plugins

class DoubleLoopCoordinator:

    def __init__(self, bidder, tracker, projection_tracker):

        self.bidder = bidder
        self.tracker = tracker
        self.projection_tracker = projection_tracker
        self.register_callbacks()

    def register_callbacks(self):

        self._register_custom_commandline_options()
        self._register_initialization_callbacks()
        self._register_before_ruc_solve_callbacks()
        self._register_before_operations_solve_callbacks()
        self._register_after_operations_callbacks()
        self._register_update_operations_stats_callbacks()
        self._register_after_ruc_activation_callbacks()
        self._register_after_simulation_callback()

        return

    def _register_custom_commandline_options(self):

        # Add command line options
        opt = Option('--track-ruc-signal',
                     help='When tracking the market signal, RUC signals are used instead of the SCED signal.',
                     action='store_true',
                     dest='track_ruc_signal',
                     default=False)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--track-sced-signal',
                     help='When tracking the market signal, SCED signals are used instead of the RUC signal.',
                     action='store_true',
                     dest='track_sced_signal',
                     default=False)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--hybrid-tracking',
                     help='When tracking the market signal, hybrid model is used.',
                     action='store_true',
                     dest='hybrid_tracking',
                     default=False)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--track-horizon',
                     help="Specifies the number of hours in the look-ahead horizon "
                          "when each tracking process is executed.",
                     action='store',
                     dest='track_horizon',
                     type='int',
                     default=48)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--bidding-generator',
                     help="Specifies the generator we derive bidding strategis for.",
                     action='store',
                     dest='bidding_generator',
                     type='string',
                     default='102_STEAM_3')
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--bidding',
                     help="Invoke generator strategic bidding when simulate.",
                     action='store_true',
                     dest='bidding',
                     default=False)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--deviation-weight',
                     help="Set the weight for deviation term when tracking",
                     action='store',
                     dest='deviation_weight',
                     type='float',
                     default=30)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--ramping-weight',
                     help="Set the weight for ramping term when tracking",
                     action='store',
                     dest='ramping_weight',
                     type='float',
                     default=20)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--cost-weight',
                     help="Set the weight for cost term when tracking",
                     action='store',
                     dest='cost_weight',
                     type='float',
                     default=1)
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--rts_gmlc_data_dir',
                     help="the relative path to rts gmlc data set",
                     action='store',
                     dest='rts_gmlc_data_dir',
                     type='str',
                     default='./RTS-GMLC/RTS_Data/SourceData/')
        prescient.plugins.add_custom_commandline_option(opt)

        opt = Option('--price_forecast_file',
                     help="the relative path to price forecasts",
                     action='store',
                     dest='price_forecast_file',
                     type='str',
                     default='../../prescient/plugins/price_forecasts/lmp_forecasts_concat.csv')
        prescient.plugins.add_custom_commandline_option(opt)

        return

    def _register_initialization_callbacks(self):
        prescient.plugins.register_initialization_callback(self.initialize_customized_results)

    def _register_before_ruc_solve_callbacks(self):
        prescient.plugins.register_before_ruc_solve_callback(self.bid_into_DAM)

    def _register_before_operations_solve_callbacks(self):
        prescient.plugins.register_before_operations_solve_callback(self.bid_into_RTM)

    def _register_after_operations_callbacks(self):
        prescient.plugins.register_after_operations_callback(self.track_sced_signal)

    def _register_update_operations_stats_callbacks(self):
        prescient.plugins.register_update_operations_stats_callback(self.update_observed_thermal_dispatch)

    def _register_after_ruc_activation_callbacks(self):
        prescient.plugins.register_after_ruc_activation_callback(self.after_ruc_activation)

    def _register_after_simulation_callback(self):
        prescient.plugins.register_after_simulation_callback(self.write_plugin_results)

    def initialize_customized_results(self, options, simulator):

        simulator.data_manager.extensions['customized_results'] = {}
        customized_results = simulator.data_manager.extensions['customized_results']

        customized_results['Generator'] = []
        customized_results['Date'] = []
        customized_results['Hour'] = []
        customized_results['State'] = []
        customized_results['RUC Schedule'] = []
        customized_results['SCED Schedule'] = []
        customized_results['Power Output'] = []

        return

    def _pass_DA_bid_to_prescient(self, options, ruc_instance, bids):

        '''
        This method passes the bids into the RUC model for day-ahead market clearing.

        Arguments:
            options: Prescient options.
            ruc_instance: Prescient RUC object
            bids: the bids we want to pass into the day-ahead market. It is a
            dictionary that has this structure {t: {generator: {power: cost}}}.
        Returns:
            None
        '''

        gen_name = self.bidder.generator
        gen_dict = ruc_instance.data['elements']['generator'][gen_name]
        p_cost = [list(bids[t][gen_name].items()) for t in range(options.ruc_horizon)]

        gen_dict['p_cost'] = {'data_type': 'time_series',
                              'values': [{'data_type' : 'cost_curve',
                                         'cost_curve_type':'piecewise',
                                         'values':p_cost[t]} for t in range(options.ruc_horizon)]
                             }
        return

    def assemble_project_tracking_signal(self, options, simulator, hour):

        '''
        This function assembles the signals for the tracking model to estimate the
        state of the bidding model at the begining of next RUC.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            hour: the simulation hour.
        Returns:
            market_signals: the market signals to be tracked.
        '''

        gen_name = self.bidder.generator

        current_ruc_dispatch = simulator.data_manager.ruc_market_active.thermal_gen_cleared_DA

        market_signals = []

        # append corresponding RUC dispatch
        for t in range(hour, hour+options.sced_horizon):
            if t >= 23:
                dispatch = current_ruc_dispatch[(gen_name,23)]
            else:
                dispatch = current_ruc_dispatch[(gen_name,t)]
            market_signals.append(dispatch)

        return market_signals

    def project_tracking_trajectory(self, options, simulator, ruc_hour):

        '''
        This function projects the full power dispatch trajectory from the
        tracking model so we can use it to update the bidding model, i.e. advance
        the time for the bidding model.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
        Returns:
            full_projected_trajectory: the full projected power dispatch trajectory.
        '''

        self._clone_tracking_model()

        for hour in range(ruc_hour, 24):

            # assemble market_signals
            market_signals = self.assemble_project_tracking_signal(options = options, \
                                                                   simulator = simulator, \
                                                                   hour = hour)
            # solve tracking
            self.projection_tracker.track_market_dispatch(market_dispatch = market_signals, \
                                                          date = current_date,\
                                                          hour = current_hour)

        # merge the trajectory
        full_projected_trajectory = {}
        for stat in self.tracker.daily_stats:
            full_projected_trajectory[stat] = self.tracker.daily_stats.get(stat) + \
                                              self.projection_tracker.daily_stats.get(stat)

        # clear the projection stats
        self.projection_tracker.daily_stats = None

        return full_projected_trajectory

    def _clone_tracking_model(self):
        self.projection_tracker.model = self.tracker.model.clone()

    def bid_into_DAM(self, options, simulator, ruc_instance, ruc_date, ruc_hour):

        '''
        This function uses the bidding objects to bid into the day-ahead market
        (DAM).

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            ruc_instance: Prescient RUC object.
            ruc_date: the date of the day-ahead market we bid into.
            ruc_hour: the hour the RUC is being solved in the day before.

        Returns:
            None
        '''

        # check if it is first day
        is_first_day = simulator.time_manager.current_time is None

        if not is_first_day:

            # solve rolling horizon to get the trajectory
            full_projected_trajectory = self.project_tracking_trajectory(options, \
                                                                         simulator, \
                                                                         ruc_hour)
            # update the bidding model
            self.bidder.update_model(**full_projected_trajectory)

        # generate bids
        bids = self.bidder.compute_bids(date = ruc_date)

        if is_first_day:
            self.current_bids = bids
            self.next_bids = bids
        else:
            self.next_bids = bids

        # pass to prescient
        self._pass_DA_bid_to_prescient(options, ruc_instance, bids)

        return

    def _pass_RT_bid_to_prescient(self, options, simulator, sced_instance, bids, hour):

        '''
        This method passes the bids into the SCED model for real-time market
        clearing.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            sced_instance: Prescient SCED object
            bids: the bids we want to pass into the real-time market. It is a
            dictionary that has this structure {t: {generator: {power: cost}}}.
            hour: the hour of the real-time market.
        Returns:
            None
        '''

        gen_name = self.bidder.generator
        gen_dict = sced_instance.data['elements']['generator'][gen_name]

        p_cost = list(bids[hour][gen_name].items())
        gen_dict['p_cost'] = {'data_type' : 'cost_curve',
                              'cost_curve_type':'piecewise',
                              'values':p_cost
                             }
        return

    def bid_into_RTM(self, options, simulator, sced_instance):

        '''
        This function bids into the real-time market. At this moment I just copy the
        corresponding day-ahead bid here.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            sced_instance: Prescient SCED object.

        Returns:
            None
        '''

        # fetch the bids
        hour = simulator.time_manager.current_time.hour
        bids = self.current_bids

        # pass bids into sced model
        self._pass_RT_bid_to_prescient(options, simulator, sced_instance, bids, hour)

        return

    def assemble_sced_tracking_market_signals(self, options,simulator,sced_instance, hour):

        '''
        This function assembles the signals for the tracking model.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            sced_instance: Prescient SCED object
            hour: the simulation hour.
        Returns:
            market_signals: the market signals to be tracked.
        '''

        gen_name = self.bidder.generator

        sced_dispatch = sced_instance.data['elements']['generator'][gen_name]['pg']['values']
        current_ruc_dispatch = simulator.data_manager.ruc_market_active.thermal_gen_cleared_DA
        if simulator.data_manager.ruc_market_pending is not None:
            next_ruc_dispatch = simulator.data_manager.ruc_market_pending.thermal_gen_cleared_DA

        # append the sced dispatch
        market_signals = [sced_dispatch[0]]

        # append corresponding RUC dispatch
        for t in range(hour+1, hour+options.sced_horizon):
            if t > 23 and simulator.data_manager.ruc_market_pending is not None:
                t = t % 24
                dispatch = next_ruc_dispatch[(gen_name,t)]
            elif t > 23 and simulator.data_manager.ruc_market_pending is None:
                dispatch = sced_dispatch[t-hour]
            else:
                dispatch = current_ruc_dispatch[(gen_name,t)]
            market_signals.append(dispatch)

        return market_signals

    def track_sced_signal(self, options, simulator, sced_instance):

        '''
        This methods uses the tracking object to track the current real-time market
        signals.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            sced_instance: Prescient SCED object

        Returns:
            None

        '''

        current_date = simulator.time_manager.current_time.date
        current_hour = simulator.time_manager.current_time.hour

        # get market signals
        market_signals = self.assemble_sced_tracking_market_signals(options = options, \
                                                                    simulator = simulator, \
                                                                    sced_instance = sced_instance, \
                                                                    hour = current_hour)

        # actual tracking
        self.tracker.track_market_dispatch(market_dispatch = market_signals, \
                                           date = current_date,\
                                           hour = current_hour)

        return

    def update_observed_thermal_dispatch(self, options, simulator, ops_stats):

        '''
        This methods extract the actual power delivered by the tracking model and
        inform Prescient, so Prescient can use this data to calculate the settlement
        and etc.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.
            ops_stats: Prescient operation statitstic object

        Returns:
            None

        '''

        ops_stats.observed_thermal_dispatch_levels[g] = self.tracker.get_last_delivered_power()

        return

    def after_ruc_activation(self, options, simulator):

        '''
        This function puts the day-ahead bids computed in the day before into effect,
        i.e. the bids for the next day become the bids for the current day.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.

        Returns:
            None
        '''

        # change bids
        current_bids = self.next_bids
        self.current_bids = current_bids
        self.next_bids = None

        return

    def write_plugin_results(self, options, simulator):

        '''
        After the simulation is completed, the plugins can write their own customized
        results. Each plugin will have to have a function named 'write_results'.

        Arguments:
            options: Prescient options.
            simulator: Prescient simulator.

        Returns:
            None

        '''

        self.bidder.write_results(path = options.output_directory)
        self.tracker.write_results(path = options.output_directory)

        return

'''
TODO:
    - Fix the remaining methods

Questions:
    -

Assumptions:
    - Trackers should provide the necessary stats to update the bidder
''''
