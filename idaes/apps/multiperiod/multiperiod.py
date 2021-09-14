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

class MultiPeriodModel:
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
    """
    def __init__(self, n_time_points, process_model_func, linking_variable_func,
                 periodic_variable_func=None):#, state_variable_func=None):
        self.n_time_points = n_time_points

        #user provided functions
        self.create_process_model = process_model_func
        self.get_linking_variable_pairs = linking_variable_func
        self.get_periodic_variable_pairs = periodic_variable_func
        #self.get_state_variable_pairs = state_variable_func

        #populated on 'build_multi_period_model'
        self._pyomo_model = None
        self._first_active_time = None

        #optional initialzation features
        #self.initialization_points = None   #library of possible initial points
        #self.initialize_func = None         #function to perform the initialize

    def build_multi_period_model(self, model_data_kwargs=None):
        """
            Build a multi-period capable model using user-provided functions

            Arguments:
                model_data_kwargs: a dict of dicts with {time:{"key",value}}
                                   where `time` is the time in the horizon. each
                                   `time` dictionary is passed to the
                                   `create_process_model` function
        """
        #use default empty dictionaries if no kwargs dict provided
        if model_data_kwargs == None:
            model_data_kwargs = {t:{} for t in range(self.n_time_points)}
        assert(list(range(len(model_data_kwargs)))==sorted(model_data_kwargs))

        m = pyo.ConcreteModel()
        m.TIME = pyo.Set(initialize=range(self.n_time_points))

        #create user defined steady-state models. Each block is a multi-period capable model.
        m.blocks = pyo.Block(m.TIME)
        for t in m.TIME:
            m.blocks[t].process = self.create_process_model(**model_data_kwargs[t])

        #link blocks together. loop over every time index except the last one
        for t in m.TIME.data()[:self.n_time_points-1]:
            link_variable_pairs = self.get_linking_variable_pairs(
                                    m.blocks[t].process,
                                    m.blocks[t+1].process)
            self._create_linking_constraints(
                                    m.blocks[t].process,
                                    link_variable_pairs)

        if self.get_periodic_variable_pairs is not None:
            N = len(m.blocks)
            periodic_variable_pairs = self.get_periodic_variable_pairs(
                                    m.blocks[N-1].process,
                                    m.blocks[0].process)
            self._create_periodic_constraints(
                                    m.blocks[N-1].process,
                                    periodic_variable_pairs)

        self._pyomo_model = m
        self._first_active_time = m.TIME.first()
        return m

    def advance_time(self, **model_data_kwargs):
        """
            Advance the current model instance to the next time period

            Arguments:
                model_data_kwargs: keyword arguments passed to user provided
                                   `create_process_model` function
        """
        m = self._pyomo_model
        previous_time = self._first_active_time
        current_time = m.TIME.next(previous_time)

        #deactivate previous time
        m.blocks[previous_time].process.deactivate()

        #track the first time in the problem horizon
        self._first_active_time = current_time

        #populate new time for the end of the horizon
        last_time = m.TIME.last()
        new_time = last_time + 1
        m.TIME.add(new_time)
        m.blocks[new_time].process = self.create_process_model(**model_data_kwargs)

        #sequential time coupling
        link_variable_pairs = self.get_linking_variable_pairs(
                                    m.blocks[last_time].process,
                                    m.blocks[new_time].process)
        self._create_linking_constraints(
                                    m.blocks[last_time].process,
                                    link_variable_pairs)

        #periodic time coupling
        if self.get_periodic_variable_pairs is not None:
            periodic_variable_pairs = self.get_periodic_variable_pairs(
                                    m.blocks[new_time].process,
                                    m.blocks[current_time].process)
            self._create_periodic_constraints(
                                    m.blocks[new_time].process,
                                    periodic_variable_pairs)
            #deactivate old periodic constraint
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
        return self._pyomo_model

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
        return [b.process for b in self._pyomo_model.blocks.values() if b.process.active]

    def _create_linking_constraints(self,b1,variable_pairs):
        """
            Create linking constraint on `b1` using `variable_pairs`
        """
        b1.link_constraints = pyo.Constraint(range(len(variable_pairs)))
        for (i,pair) in enumerate(variable_pairs):
            b1.link_constraints[i] = pair[0]==pair[1]

    def _create_periodic_constraints(self,b1,variable_pairs):
        """
            Create periodic linking constraint on `b1` using `variable_pairs`
        """
        b1.periodic_constraints = pyo.Constraint(range(len(variable_pairs)))
        for (i,pair) in enumerate(variable_pairs):
            b1.periodic_constraints[i] = pair[0]==pair[1]
