import prescient.plugins as pplugins
from prescient.simulator.config import PrescientConfig
from pyomo.common.config import ConfigDict, ConfigValue

class DoubleLoopCoordinator:

    def __init__(self, bidder, tracker, projection_tracker):

        '''
        Initializes the DoubleLoopCoordinator object and registers functionalities
        in Prescient's plugin system.

        Arguments:
            bidder: an initialized bidder object
            tracker: an initialized bidder object
            projection_tracker: an initialized bidder object, this object is for
                                projecting the system states and updating bidder
                                model

        Returns:
            None
        '''

        self.bidder = bidder
        self.tracker = tracker
        self.projection_tracker = projection_tracker

    def register_plugins(self, \
                         context: pplugins.PluginRegistrationContext,\
                         options: PrescientConfig,\
                         plugin_config: ConfigDict):

        '''
        Register functionalities in Prescient's plugin system.

        Arguments:
            None

        Returns:
            None
        '''

        self.plugin_config = plugin_config

        self._register_initialization_callbacks(context, options, plugin_config)
        self._register_before_ruc_solve_callbacks(context, options, plugin_config)
        self._register_before_operations_solve_callbacks(context, options, plugin_config)
        self._register_after_operations_callbacks(context, options, plugin_config)
        self._register_update_operations_stats_callbacks(context, options, plugin_config)
        self._register_after_ruc_activation_callbacks(context, options, plugin_config)
        self._register_finalization_callbacks(context, options, plugin_config)

        return

    def get_configuration(self, key):

        '''
        Register customized commandline options.

        Arguments:
            None

        Returns:
            None
        '''

        config = ConfigDict()

        # Add command line options
        config.declare('bidding_generator',
                       ConfigValue(domain = str,
                                   description = 'Specifies the generator we derive bidding strategis for.',
                                   default = None)).declare_as_argument('--bidding-generator')
        ## How to access this option?
        # options.plugin.{key}.bidding_generator

        return config

    def _register_initialization_callbacks(self, \
                                           context: pplugins.PluginRegistrationContext,\
                                           options: PrescientConfig,\
                                           plugin_config: ConfigDict):

        '''
        Register initialization plugins, which run before Prescient simulation
        begins.

        Arguments:
            None

        Returns:
            None
        '''

        context.register_initialization_callback(self.initialize_customized_results)

    def _register_before_ruc_solve_callbacks(self, \
                                             context: pplugins.PluginRegistrationContext,\
                                             options: PrescientConfig,\
                                             plugin_config: ConfigDict):

        '''
        Register plugins that run before Prescient solves Reliability Unit
        Commitment (RUC) problems.

        Arguments:
            None

        Returns:
            None
        '''

        context.register_before_ruc_solve_callback(self.bid_into_DAM)

    def _register_before_operations_solve_callbacks(self, \
                                                    context: pplugins.PluginRegistrationContext,\
                                                    options: PrescientConfig,\
                                                    plugin_config: ConfigDict):

        '''
        Register plugins that run before Prescient solves Securitiy Constrained
        Economic Dispatch (SCED), aka "operation", problems.

        Arguments:
            None

        Returns:
            None
        '''

        context.register_before_operations_solve_callback(self.bid_into_RTM)

    def _register_after_operations_callbacks(self, \
                                             context: pplugins.PluginRegistrationContext,\
                                             options: PrescientConfig,\
                                             plugin_config: ConfigDict):

        '''
        Register plugins that run after Prescient solves Securitiy Constrained
        Economic Dispatch (SCED), aka "operation", problems.

        Arguments:
            None

        Returns:
            None
        '''

        context.register_after_operations_callback(self.track_sced_signal)

    def _register_update_operations_stats_callbacks(self, \
                                                    context: pplugins.PluginRegistrationContext,\
                                                    options: PrescientConfig,\
                                                    plugin_config: ConfigDict):

        '''
        Register plugins that update stats of Securitiy Constrained Economic
        Dispatch (SCED), aka "operation".

        Arguments:
            None

        Returns:
            None
        '''

        context.register_update_operations_stats_callback(self.update_observed_thermal_dispatch)

    def _register_after_ruc_activation_callbacks(self, \
                                                 context: pplugins.PluginRegistrationContext,\
                                                 options: PrescientConfig,\
                                                 plugin_config: ConfigDict):

        '''
        Register plugins that update stats of Securitiy Constrained Economic
        Dispatch (SCED), aka "operation".

        Arguments:
            None

        Returns:
            None
        '''

        context.register_after_ruc_activation_callback(self.activate_DA_bids)

    def _register_finalization_callbacks(self, \
                                         context: pplugins.PluginRegistrationContext,\
                                         options: PrescientConfig,\
                                         plugin_config: ConfigDict):

        '''
        Register finalization plugins, which run after Prescient simulation
        finishes.

        Arguments:
            None

        Returns:
            None
        '''

        context.register_finalization_callback(self.write_plugin_results)

    def initialize_customized_results(self, options, simulator):

        '''
        This method is outdated.
        '''

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
        sced_horizon = options.sced_horizon

        market_signals = self._assemble_project_tracking_signal(gen_name = gen_name,\
                                                                current_ruc_dispatch = current_ruc_dispatch,\
                                                                hour = hour,\
                                                                sced_horizon = sced_horizon)
        return market_signals

    @staticmethod
    def _assemble_project_tracking_signal(gen_name, current_ruc_dispatch, hour, sced_horizon):

        market_signals = []
        # append corresponding RUC dispatch
        for t in range(hour, hour + sced_horizon):
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
            ruc_hour: the hour RUC is being solved
        Returns:
            full_projected_trajectory: the full projected power dispatch trajectory.
        '''

        current_date = simulator.time_manager.current_time.date
        current_hour = simulator.time_manager.current_time.hour

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
        '''
        Clone the model in tracker and replace that of projection tracker. In this
        way, tracker and projection tracker have the same states before projection.

        Arguments:
            None

        Returns:
            None
        '''
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
        sced_horizon = options.sced_horizon
        current_ruc_dispatch = simulator.data_manager.ruc_market_active.thermal_gen_cleared_DA
        if simulator.data_manager.ruc_market_pending is not None:
            next_ruc_dispatch = simulator.data_manager.ruc_market_pending.thermal_gen_cleared_DA
        else:
            next_ruc_dispatch = None

        market_signals = self._assemble_sced_tracking_market_signals(gen_name = gen_name, \
                                                                     hour = hour, \
                                                                     sced_dispatch = sced_dispatch, \
                                                                     sced_horizon = sced_horizon, \
                                                                     current_ruc_dispatch = current_ruc_dispatch, \
                                                                     next_ruc_dispatch = None)
        return market_signals

    @staticmethod
    def _assemble_sced_tracking_market_signals(gen_name, \
                                               hour, \
                                               sced_dispatch, \
                                               sced_horizon, \
                                               current_ruc_dispatch, \
                                               next_ruc_dispatch = None):

        # append the sced dispatch
        market_signals = [sced_dispatch[0]]

        # append corresponding RUC dispatch
        for t in range(hour+1, hour + sced_horizon):
            if t > 23 and next_ruc_dispatch:
                t = t % 24
                dispatch = next_ruc_dispatch[(gen_name,t)]
            elif t > 23 and next_ruc_dispatch is None:
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
        g = self.bidder.generator
        ops_stats.observed_thermal_dispatch_levels[g] = self.tracker.get_last_delivered_power()

        return

    def activate_DA_bids(self, options, simulator):

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
