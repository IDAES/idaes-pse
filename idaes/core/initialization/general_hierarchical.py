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
"""
Initializer class for implementing Hierarchical initialization routines for
IDAES models with standard forms (e.g. units with 1 control volume)
"""
from pyomo.environ import Block
from pyomo.common.config import ConfigValue

from idaes.core.initialization.initializer_base import InitializerBase
from idaes.core.solvers import get_solver
from idaes.core.util import to_json, from_json, StoreSpec
import idaes.logger as idaeslog

__author__ = "Andrew Lee"


class SingleControlVolumeUnitInitializer(InitializerBase):
    """
    This is a general purpose Initializer object for unit models which
    have a single Control Volume named 'control_volume'.
    """

    CONFIG = InitializerBase.CONFIG()
    CONFIG.declare(
        "solver",
        ConfigValue(
            default=None,
            description="Solver to use for initialization",
        ),
    )
    CONFIG.declare(
        "solver_options",
        ConfigValue(
            default={},
            description="Dict of options to pass to solver",
        ),
    )

    def initialization_routine(self, model: Block):
        """
        Common initialization routine for models with standard form.

        Args:
            model: Pyomo Block to be initialized

        Returns:
            Pyomo solver results object
        """
        if not hasattr(model, "control_volume"):
            raise TypeError(
                f"Model {model.name} does not appear to be a standard form unit model. "
                f"Please use an Initializer specific to the model being initialized."
            )

        # Get logger
        _log = self.get_logger(model)

        # Initialize properties
        try:
            # Guess a 0-D control volume
            prop_in = model.control_volume.properties_in

            # Initialize inlet properties - inlet state should already be fixed
            # TODO: Call an initializer for the StateBlock
            prop_in.initialize(
                solver=self.config.solver,
                optarg=self.config.solver_options,
                outlvl=self.config.output_level,
            )

            # Map solution from inlet properties to outlet properties
            state = to_json(
                model.control_volume.properties_in,
                wts=StoreSpec().value(),
                return_dict=True,
            )
            from_json(
                model.control_volume.properties_out,
                sd=state,
                wts=StoreSpec().value(only_not_fixed=True),
            )
        except AttributeError:
            # Assume it must be a 1-D control volume
            # TODO: Add steps here
            pass

        _log.info_high("Step 1: properties initialization complete")

        # Initialize reactions if they exist
        try:
            model.control_volume.reactions.initialize(
                solver=self.config.solver,
                optarg=self.config.solver_options,
                outlvl=self.config.output_level,
                state_vars_fixed=True,
            )
        except AttributeError:
            pass
        _log.info_high("Step 2: reactions initialization complete")

        # Solve full unit model
        solver = get_solver(self.config.solver, self.config.solver_options)
        solve_log = idaeslog.getSolveLogger(
            model.name, self.config.output_level, tag="unit"
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)

        _log.info_high("Step 3: {}.".format(idaeslog.condition(results)))

        return results
