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
import pyomo.environ as pyo
from pyomo.common.timing import TicTocTimer
from idaes.core.solvers import get_solver
from idaes.core.util import from_json, to_json
import matplotlib.pyplot as plt
import logging

_logger = logging.getLogger(__name__)


class MultiPeriodModel(pyo.ConcreteModel):
    """
    The `MultiPeriodModel` class helps transfer existing steady-state
    process models to multiperiod versions that contain dynamic time coupling.

    Arguments:
        n_time_points: number of points to use in time horizon
        process_model_func: function that returns a multiperiod capable pyomo model
        linking_variable_func: function that returns a tuple of variable
                               pairs to link between time steps
        periodic_variable_func: a function that returns a tuple of variable
                                pairs to link between last and first time steps
        use_stochastic_build: Uses `build_stochastic_multi_period` method if set to True
        set_days: list containing the set of representative days
        set_years: list containing the set of years
        set_scenarios: list containing the set of scenarios
        initialization_func: function that fixes the degrees of freedom and initializes an
                             instance of the flowsheet.
        unfix_dof_func: function that unfixes a few degrees of freedom for optimization
        flowsheet_options: dictionary containing the arguments needed for `process_model_func`
        initialization_options: dictionary contaning the arguments needed for `initialization_func`
        unfix_dof_options: dictionary containing the arguments needed for `unfix_dof_func`
        solver: pyomo solver object
        outlvl: logging level

    Returns:
        (stochastic) multi-period optimization model
    """

    def __init__(
        self,
        n_time_points,
        process_model_func,
        linking_variable_func,
        periodic_variable_func=None,
        use_stochastic_build=False,
        set_days=None,
        set_years=None,
        set_scenarios=None,
        initialization_func=None,
        unfix_dof_func=None,
        flowsheet_options={},
        initialization_options={},
        unfix_dof_options={},
        solver=None,
        outlvl=logging.WARNING,
    ):  # , state_variable_func=None):

        super().__init__()

        self.n_time_points = n_time_points

        # user provided functions
        self.create_process_model = process_model_func
        self.get_linking_variable_pairs = linking_variable_func
        self.get_periodic_variable_pairs = periodic_variable_func
        self.initialization_func = initialization_func
        self.unfix_dof_func = unfix_dof_func
        # self.get_state_variable_pairs = state_variable_func

        # populated on 'build_multi_period_model'
        self._first_active_time = None

        # Create sets
        if use_stochastic_build:
            self.set_time = pyo.RangeSet(n_time_points)

            if set_days is not None:
                self.set_days = pyo.Set(initialize=set_days)
                self._multiple_days = True

            else:
                self._multiple_days = False

            if set_years is not None:
                self.set_years = pyo.Set(initialize=set_years)
                self._multiyear = True

            else:
                self._multiyear = False

            if set_scenarios is not None:
                self.set_scenarios = pyo.Set(initialize=set_scenarios)
                self._stochastic_model = True

            else:
                self._stochastic_model = False

            if solver is None:
                solver = get_solver()

            _logger = logging.getLogger(__name__)
            _logger.setLevel(outlvl)

            # Build the stochastic multiperiod optimization model
            self.build_stochastic_multi_period(
                flowsheet_options,
                initialization_options,
                unfix_dof_options,
                solver,
            )

        # optional initialzation features
        # self.initialization_points = None   #library of possible initial points
        # self.initialize_func = None         #function to perform the initialize

    def build_multi_period_model(
        self,
        model_data_kwargs=None,
        flowsheet_options={},
        initialization_options={},
        unfix_dof_options={},
        solver=None,
    ):
        """
        Build a multi-period capable model using user-provided functions

        Arguments:
            model_data_kwargs: a dict of dicts with {time:{"key",value}}
                               where `time` is the time in the horizon. each
                               `time` dictionary is passed to the
                               `create_process_model` function
            flowsheet_options: dict containing the arguments needed to build an instance of flowsheet
            initialization_options: dict containing the arguments needed for `initialization_func`
            unfix_dof_options: dict containing the arguments needed for `unfix_dof_func`
            solver: pyomo solver object
        """
        # use default empty dictionaries if no kwargs dict provided
        if model_data_kwargs == None:
            model_data_kwargs = {t: {} for t in range(self.n_time_points)}
        assert list(range(len(model_data_kwargs))) == sorted(model_data_kwargs)

        m = self
        m.TIME = pyo.Set(initialize=range(self.n_time_points))

        # Begin the formulation of the multiperiod optimization problem
        timer = TicTocTimer()  # Create timer object
        timer.toc("Beginning the formulation of the multiperiod optimization problem.")

        # create user defined steady-state models. Each block is a multi-period capable model.
        m.blocks = pyo.Block(m.TIME)
        for t in m.TIME:
            _logger.info(f"...Constructing the flowsheet model for {m.blocks[t].name}")
            m.blocks[t].process = self.create_process_model(**model_data_kwargs[t])

        timer.toc("Completed the formulation of the multiperiod optimization problem.")

        if solver is None:
            solver = get_solver()

        self.initialize_multi_period_model(
            flowsheet_options, initialization_options, solver, False, timer
        )

        self.unfix_degrees_of_freedom(unfix_dof_options, False, timer)

        # link blocks together. loop over every time index except the last one
        for t in m.TIME.data()[: self.n_time_points - 1]:
            link_variable_pairs = self.get_linking_variable_pairs(
                m.blocks[t].process, m.blocks[t + 1].process
            )
            self._create_linking_constraints(m.blocks[t].process, link_variable_pairs)

        if self.get_periodic_variable_pairs is not None:
            N = len(m.blocks)
            periodic_variable_pairs = self.get_periodic_variable_pairs(
                m.blocks[N - 1].process, m.blocks[0].process
            )
            self._create_periodic_constraints(
                m.blocks[N - 1].process, periodic_variable_pairs
            )

        self._first_active_time = m.TIME.first()
        return m

    def advance_time(self, **model_data_kwargs):
        """
        Advance the current model instance to the next time period

        Arguments:
            model_data_kwargs: keyword arguments passed to user provided
                               `create_process_model` function
        """
        m = self
        previous_time = self._first_active_time
        current_time = m.TIME.next(previous_time)

        # deactivate previous time
        m.blocks[previous_time].process.deactivate()

        # track the first time in the problem horizon
        self._first_active_time = current_time

        # populate new time for the end of the horizon
        last_time = m.TIME.last()
        new_time = last_time + 1
        m.TIME.add(new_time)
        m.blocks[new_time].process = self.create_process_model(**model_data_kwargs)

        # sequential time coupling
        link_variable_pairs = self.get_linking_variable_pairs(
            m.blocks[last_time].process, m.blocks[new_time].process
        )
        self._create_linking_constraints(
            m.blocks[last_time].process, link_variable_pairs
        )

        # periodic time coupling
        if self.get_periodic_variable_pairs is not None:
            periodic_variable_pairs = self.get_periodic_variable_pairs(
                m.blocks[new_time].process, m.blocks[current_time].process
            )
            self._create_periodic_constraints(
                m.blocks[new_time].process, periodic_variable_pairs
            )
            # deactivate old periodic constraint
            m.blocks[last_time].process.periodic_constraints.deactivate()

        # TODO: discuss where state goes.
        # sometimes the user might want to fix values based on a 'real' process
        # also TODO: inspect argument and use fix() if possible
        # if self.get_state_variable_pairs is not None:
        #     state_variable_pairs = self.get_state_variable_pairs(
        #                       m.blocks[previous_time].process,
        #                       m.blocks[current_time].process)
        #     self._fix_initial_states(
        #                       m.blocks[current_time].process,
        #                       state_variable_pairs)

    @property
    def pyomo_model(self):
        """
        Retrieve the underlying pyomo model
        """
        return self

    @property
    def current_time(self):
        """
        Retrieve the current multiperiod model time
        """
        return self._first_active_time

    def get_active_process_blocks(self):
        """
        Retrieve the active time blocks of the pyomo model
        """
        return [b.process for b in self.blocks.values() if b.process.active]

    def _create_linking_constraints(self, b1, variable_pairs):
        """
        Create linking constraint on `b1` using `variable_pairs`
        """
        b1.link_constraints = pyo.Constraint(range(len(variable_pairs)))
        for (i, pair) in enumerate(variable_pairs):
            b1.link_constraints[i] = pair[0] == pair[1]

    def _create_periodic_constraints(self, b1, variable_pairs):
        """
        Create periodic linking constraint on `b1` using `variable_pairs`
        """
        b1.periodic_constraints = pyo.Constraint(range(len(variable_pairs)))
        for (i, pair) in enumerate(variable_pairs):
            b1.periodic_constraints[i] = pair[0] == pair[1]

    def build_stochastic_multi_period(
        self,
        flowsheet_options,
        initialization_options,
        unfix_dof_options,
        solver,
    ):
        """
        This function constructs the stochastic multiperiod optimization problem
        """
        # Create set of periods:
        multiple_days = self._multiple_days
        multiyear = self._multiyear

        # Construct the set of time periods
        if multiyear and multiple_days:
            set_period = [
                (t, d, y)
                for y in self.set_years
                for d in self.set_days
                for t in self.set_time
            ]

        elif multiyear and not multiple_days:
            set_period = [(t, y) for y in self.set_years for t in self.set_time]

        elif not multiyear and multiple_days:
            set_period = [(t, d) for d in self.set_days for t in self.set_time]

        else:
            set_period = [t for t in self.set_time]

        self.set_period = pyo.Set(initialize=set_period)

        # Define a function to create a multiperiod model for one scenario
        def _build_scenario_model(m):
            m.period = pyo.Block(self.set_period)

            for i in m.period:
                _logger.info(
                    f"...Constructing the flowsheet model for {m.period[i].name}"
                )
                self.create_process_model(m.period[i], **flowsheet_options)

            # link blocks together. loop over every time index except the last one
            if self.get_linking_variable_pairs is None:
                _logger.warning(
                    f"linking_variable_func is not provided, so variables across"
                    f" time periods are not linked."
                )
                return

            # TODO: Instead of accessing the constraints as
            # m.link_constraints[t, d, y].link_constraints[1], add a Reference for easy access.
            if multiyear and multiple_days:
                m.link_constraints = pyo.Block(
                    self.set_time.data()[:-1], self.set_days, self.set_years
                )

                for y in self.set_years:
                    for d in self.set_days:
                        for t in self.set_time.data()[:-1]:
                            link_variable_pairs = self.get_linking_variable_pairs(
                                m.period[t, d, y], m.period[t + 1, d, y]
                            )
                            self._create_linking_constraints(
                                m.link_constraints[t, d, y], link_variable_pairs
                            )

            elif multiyear and not multiple_days:
                m.link_constraints = pyo.Block(
                    self.set_time.data()[:-1], self.set_years
                )

                for y in self.set_years:
                    for t in self.set_time.data()[:-1]:
                        link_variable_pairs = self.get_linking_variable_pairs(
                            m.period[t, y], m.period[t + 1, y]
                        )
                        self._create_linking_constraints(
                            m.link_constraints[t, y], link_variable_pairs
                        )

            elif not multiyear and multiple_days:
                m.link_constraints = pyo.Block(self.set_time.data()[:-1], self.set_days)

                for d in self.set_days:
                    for t in self.set_time.data()[:-1]:
                        link_variable_pairs = self.get_linking_variable_pairs(
                            m.period[t, d], m.period[t + 1, d]
                        )
                        self._create_linking_constraints(
                            m.link_constraints[t, d], link_variable_pairs
                        )

            else:
                m.link_constraints = pyo.Block(self.set_time.data()[:-1])

                for t in self.set_time.data()[:-1]:
                    link_variable_pairs = self.get_linking_variable_pairs(
                        m.period[t], m.period[t + 1]
                    )
                    self._create_linking_constraints(
                        m.link_constraints[t], link_variable_pairs
                    )

            # Check if a method for periodic constraints is given
            if self.get_periodic_variable_pairs is not None:
                _logger.warning(
                    f"A method is provided for get_periodic_variable_pairs. "
                    f"build_stochastic_multi_period method does not support periodic "
                    f"constraints, so the user needs to add them manually."
                )

        # Begin the formulation of the multiperiod optimization problem
        timer = TicTocTimer()  # Create timer object
        timer.toc("Beginning the formulation of the multiperiod optimization problem.")

        if self._stochastic_model:
            self.scenario = pyo.Block(self.set_scenarios)

            for i in self.scenario:
                _logger.info(f"Constructing the model for scenario {i}")
                _build_scenario_model(self.scenario[i])

        else:
            _build_scenario_model(self)

        timer.toc("Completed the formulation of the multiperiod optimization problem.")

        self.initialize_multi_period_model(
            flowsheet_options, initialization_options, solver, True, timer
        )

        self.unfix_degrees_of_freedom(unfix_dof_options, True, timer)

    def initialize_multi_period_model(
        self,
        flowsheet_options,
        initialization_options,
        solver,
        use_stochastic_build,
        timer,
    ):
        """
        This function initializes the entire multi-period model.

        Args:
            flowsheet_options: dict containing the arguments needed to create an instance of a flowsheet
            initialization_options: dict containing the arguments needed to initialize the flowsheet
            solver: pyomo solver object
            use_stochastic_build: Set it to True if `build_stochastic_multi_period` method is used
            time: Timer object
        """
        if self.initialization_func is None:
            _logger.warning(
                f"Initialization function is not provided. "
                f"Returning the multiperiod model without initialization."
            )
            return

        blk = pyo.ConcreteModel()
        self.create_process_model(blk, **flowsheet_options)
        self.initialization_func(blk, **initialization_options)
        result = solver.solve(blk)

        if not pyo.check_optimal_termination(result):
            raise Exception(
                f"Flowsheet did not converge to optimality after fixing the degrees of freedom. "
                f"To create the multi-period model without initialization, do not provide "
                f"initialization_func argument."
            )

        # Store the initialized model in `init_model` object
        init_model = to_json(blk, return_dict=True)
        timer.toc("Created an instance of the flowsheet and initialized it.")

        # Initialize the multiperiod optimization model
        if use_stochastic_build:
            if self._stochastic_model:
                for s in self.set_scenarios:
                    for p in self.scenario[s].period:
                        from_json(self.scenario[s].period[p], sd=init_model)

            else:
                for p in self.period:
                    from_json(self.period[p], sd=init_model)

        else:
            blks = self.get_active_process_blocks()
            for blk in blks:
                from_json(blk, init_model)

        timer.toc("Initialized the entire multiperiod optimization model.")

    def unfix_degrees_of_freedom(
        self,
        unfix_dof_options,
        use_stochastic_build,
        timer,
    ):
        """
        Unfix the degrees of freedom in each period model for optimization model

        Args:
            unfix_dof_options: dict containing the arguments needed for `unfix_dof_func`
            use_stochastic_build: Set it to True if `build_stochastic_multi_period` method is used
            time: Timer object
        """
        if self.unfix_dof_func is None:
            _logger.warning(
                f"unfix_dof function is not provided. "
                f"Returning the model without unfixing degrees of freedom"
            )
            return

        if use_stochastic_build:
            if self._stochastic_model:
                for s in self.set_scenarios:
                    for p in self.scenario[s].period:
                        self.unfix_dof_func(
                            self.scenario[s].period[p], **unfix_dof_options
                        )

            else:
                for p in self.period:
                    self.unfix_dof_func(self.period[p], **unfix_dof_options)

        else:
            blks = self.get_active_process_blocks()
            for blk in blks:
                self.unfix_dof_func(blk, **unfix_dof_options)

        timer.toc("Unfixed the degrees of freedom from each period model.")

    @staticmethod
    def plot_lmp_signal(
        lmp, time=None, draw_style="steps", x_range=None, y_range=None, grid=None
    ):
        """
        This function plots LMP signals as a function of time.

        Args:
            lmp: list or dict of LMP signals (length of dictionary must be <= 6).
            time: list or dict of time. Taken to be 1:len(lmp) if unspecified.
            draw_style: Plot style. Should be either "steps" or "default".
            x_range: tuple or dict of tuples with x-axis range.
            y_range: tuple or dict of tuples with y-axis range.
            grid: Grid shape of the plot.

        Returns:
            Plot of LMP signals
        """

        if type(lmp) is list:
            grid = (1, 1)

            if time is None:
                plt_time = {1: [i for i in range(1, len(lmp) + 1)]}

            else:
                plt_time = {1: time}

            plt_lmp = {1: lmp}
            plt_title = {1: ""}
            color = {1: "tab:red"}
            plt_x_range = {1: x_range}
            plt_y_range = {1: y_range}

        elif type(lmp) is dict:
            if len(lmp) > 6:
                raise Exception(
                    "Number of LMP signals provided exceeds six: the maximum "
                    "number of subplots the function supports."
                )

            grid_shape = {
                1: (1, 1),
                2: (1, 2),
                3: (2, 2),
                4: (2, 2),
                5: (2, 3),
                6: (2, 3),
            }
            if grid is None:
                grid = grid_shape[len(lmp)]

            color = {
                1: "tab:red",
                2: "tab:red",
                3: "tab:red",
                4: "tab:red",
                5: "tab:red",
                6: "tab:red",
            }

            plt_time = {}
            plt_lmp = {}
            plt_title = {}
            plt_x_range = {}
            plt_y_range = {}

            counter = 1
            for i in lmp:
                plt_lmp[counter] = lmp[i]

                if time is None:
                    plt_time[counter] = [j for j in range(1, len(lmp[i]) + 1)]

                else:
                    plt_time[counter] = time[i]

                plt_title[counter] = str(i)
                plt_x_range[counter] = (
                    None if x_range is None or i not in x_range else x_range[i]
                )
                plt_y_range[counter] = (
                    None if y_range is None or i not in y_range else y_range[i]
                )
                counter = counter + 1

        fig = plt.figure()

        for i in range(1, len(plt_lmp) + 1):
            ax = fig.add_subplot(grid[0], grid[1], i)
            ax.plot(plt_time[i], plt_lmp[i], color=color[i], drawstyle=draw_style)
            ax.set_xlabel("time (hr)")
            ax.set_ylabel("LMP ($/MWh)")
            ax.set_title(plt_title[i])

            if plt_x_range[i] is not None:
                ax.set_xlim(plt_x_range[i][0], plt_x_range[i][1])

            if plt_y_range[i] is not None:
                ax.set_ylim(plt_y_range[i][0], plt_y_range[i][1])

        fig.tight_layout()
        plt.show()

    @staticmethod
    def plot_lmp_and_schedule(
        lmp,
        schedule,
        time=None,
        y_label=None,
        x_range=None,
        lmp_range=None,
        y_range=None,
        x_label="time (hr)",
        lmp_label="LMP ($/MWh)",
        color=None,
        draw_style="steps",
        grid=None,
    ):

        """
        The function plots optimal operation schedule as a function of time.

        Args:
            lmp: list of LMPs.
            schedule: dict of operating schedules {"variable": [profile], ...}.
            time: list of time instances. Taken to be 1:len(lmp) if unspecified.
            y_label: dict of y labels for the schedule. Taken to be keys of schedule if unspecified.
            x_range: tuple containing the range of x-axis.
            lmp_range: tuple containing the range of lmp-axis.
            y_range: dict of tuples containing the range of y-axis.
            x_label: x-axis label. "time (hr)" is the default
            lmp_label: lmp-axis label. "LMP ($/MWh)" is the default.
            color: dict of colors for the plots.
            draw_style: plot sytle. Must be either "steps" or "default".
            grid: grid shape of the plost.

        Returns:
            None
        """
        if len(schedule) > 4:
            raise Exception(
                f"Number of elements in schedule exceeds four: "
                f"the maximum number of subplots the function supports."
            )

        key_list = {index + 1: value for index, value in enumerate(schedule)}

        if time is None:
            time = [i for i in range(1, len(lmp) + 1)]

        if grid is None:
            grid_shape = {1: (1, 1), 2: (1, 2), 3: (2, 2), 4: (2, 2)}
            grid = grid_shape[len(schedule)]

        lmp_color = "tab:red"
        if color is None:
            plt_color = {1: "tab:blue", 2: "magenta", 3: "tab:green", 4: "tab:cyan"}
        else:
            plt_color = {index + 1: color[value] for index, value in enumerate(color)}

        fig = plt.figure()

        for i in range(1, len(schedule) + 1):
            ax = fig.add_subplot(grid[0], grid[1], i)
            ax.set_xlabel(x_label)
            ax.set_ylabel(lmp_label, color=lmp_color)
            ax.plot(time, lmp, color=lmp_color, drawstyle="steps")
            ax.tick_params(axis="y", labelcolor=lmp_color)

            if x_range is not None:
                ax.set_xlim(x_range[0], x_range[1])

            if lmp_range is not None:
                ax.set_ylim(lmp_range[0], lmp_range[1])

            ax1 = ax.twinx()
            ax1.plot(
                time, schedule[key_list[i]], color=plt_color[i], drawstyle=draw_style
            )
            ax1.tick_params(axis="y", labelcolor=plt_color[i])

            if y_label is not None and key_list[i] in y_label:
                ax1.set_ylabel(y_label[key_list[i]], color=plt_color[i])

            else:
                ax1.set_ylabel(str(key_list[i]), color=plt_color[i])

            if y_range is not None and key_list[i] in y_range:
                ax1.set_ylim(y_range[key_list[i]][0], y_range[key_list[i]][1])

        fig.tight_layout()
        plt.show()
