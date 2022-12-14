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

from idaes.core.initialization.initializer_base import ModularInitializerBase
from idaes.core.solvers import get_solver
from idaes.core.util import to_json, from_json, StoreSpec
from idaes.core.util.exceptions import InitializationError
import idaes.logger as idaeslog

__author__ = "Andrew Lee"


class SingleControlVolumeUnitInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer object for unit models which
    have a single Control Volume named 'control_volume'.
    """

    CONFIG = ModularInitializerBase.CONFIG()
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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._solver = None

    def initialization_routine(self, model: Block, addon_args: dict = {}):
        """
        Common initialization routine for models with standard form.

        Args:
            model: Pyomo Block to be initialized
            addon_args: dict of arguments to be passed to add-on Initializers. Keys should be submodel components.

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

        # Prepare add-ons for initialization
        sub_initializers, addon_args = self._prepare_addons(model, addon_args, _log)

        # Initialize model and sub-models
        results = self._initialize_submodels(model, addon_args, _log)

        # Solve full model including add-ons
        results = self._solve_full_model(model, _log, results)

        # Clean up add-ons
        self._cleanup(model, addon_args, sub_initializers, _log)

        return results

    def _prepare_addons(self, model, addon_args, logger):
        sub_initializers = {}
        addon_args = dict(addon_args)
        for sm in model.initialization_order:
            if sm is not model:
                # Get initializers for add-ons
                sub_initializers[sm] = self.get_submodel_initializer(sm)()

                if sm not in addon_args.keys():
                    addon_args[sm] = {}

                # Call prepare method for add-ons
                # TODO: Need arguments for submodel initialization
                sub_initializers[sm].addon_prepare(sm, **addon_args[sm])

        logger.info_high("Step 1: preparation complete.")

        return sub_initializers, addon_args

    def _initialize_submodels(self, model, addon_args, logger):
        results = None

        for sm in model.initialization_order:
            print(sm.name)
            if sm is model:
                results = self._initialize_main_model(model, logger)
            else:
                # TODO: Need arguments for submodel initialization
                sub_initializers[sm].addon_initialize(sm, **addon_args[sm])

        if results is None:
            raise InitializationError(
                f"Main model ({model.name}) was not initialized (no results returned). "
                f"This is likely due to an error in the model.initialization_order."
            )

        logger.info_high("Step 2: sub-model initialization complete.")

        return results

    def _initialize_main_model(self, model, logger):
        # Initialize properties
        try:
            # Guess a 0-D control volume
            # Initialize inlet properties - inlet state should already be fixed
            prop_init = self.get_submodel_initializer(
                model.control_volume.properties_in
            )

            if prop_init is not None:
                prop_init.initialize(
                    model.control_volume.properties_in,
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
            raise

        logger.info_high("Step 2a: properties initialization complete")

        # Initialize reactions if they exist
        if hasattr(model.control_volume, "reactions"):
            rxn_init = self.get_submodel_initializer(model.control_volume.reactions)

            if rxn_init is not None:
                rxn_init.initialize(
                    model.control_volume.reactions,
                    solver=self.config.solver,
                    optarg=self.config.solver_options,
                    outlvl=self.config.output_level,
                )
        logger.info_high("Step 2b: reactions initialization complete")

        # Solve main model
        solver = get_solver(self.config.solver, self.config.solver_options)
        solve_log = idaeslog.getSolveLogger(
            model.name, self.config.output_level, tag="unit"
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = self._get_solver().solve(model, tee=slc.tee)

        logger.info_high("Step 2c: model {}.".format(idaeslog.condition(results)))

        return results

    def _solve_full_model(self, model, logger, results):
        # Check to see if there are any add-ons
        if len(model.initialization_order) > 1:
            # Solve model with addons
            solve_log = idaeslog.getSolveLogger(
                model.name, self.config.output_level, tag="unit"
            )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = self._get_solver().solve(model, tee=slc.tee)

        logger.info_high(
            "Step 3: full model initialization {}.".format(idaeslog.condition(results))
        )

        return results

    def _cleanup(self, model, addon_args, sub_initializers, logger):
        for sm in reversed(model.initialization_order):
            if sm is not model:
                # TODO: Need arguments for submodel initialization
                sub_initializers[sm].addon_finalize(sm, **addon_args[sm])

        logger.info_high("Step 4: clean up completed.")

    def _get_solver(self):
        if self._solver is None:
            self._solver = get_solver(self.config.solver, self.config.solver_options)

        return self._solver
