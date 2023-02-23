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
from pyomo.common.config import ConfigDict, ConfigValue

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

    For details of the initialization routine, please see the documentation.

    """

    CONFIG = ModularInitializerBase.CONFIG()
    CONFIG.declare(
        "solver",
        ConfigValue(
            default=None,  # TODO: Can we add a square problem solver as the default here?
            # At the moment there is an issue with the scipy solvers not supporting the tee argument.
            description="Solver to use for initialization",
        ),
    )
    CONFIG.declare(
        "solver_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to solver",
        ),
    )

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._solver = None

    def initialization_routine(
        self,
        model: Block,
        plugin_initializer_args: dict = None,
        copy_inlet_state: bool = False,
    ):
        """
        Common initialization routine for models with one control volume.

        Args:
            model: Pyomo Block to be initialized
            plugin_initializer_args: dict-of-dicts containing arguments to be passed to plug-in Initializers.
                Keys should be submodel components.
            copy_inlet_state: bool (default=False). Whether to copy inlet state to other sttes or not.
                Copying will generally be faster, but inlet states may not contain all properties
                required elsewhere.

        Returns:
            Pyomo solver results object
        """
        if not hasattr(model, "control_volume"):
            raise TypeError(
                f"Model {model.name} does not appear to be a standard form unit model. "
                f"Please use an Initializer specific to the model being initialized."
            )
        if plugin_initializer_args is None:
            plugin_initializer_args = {}

        # Get logger
        _log = self.get_logger(model)

        # Prepare plug-ins for initialization
        sub_initializers, plugin_initializer_args = self._prepare_plugins(
            model, plugin_initializer_args, _log
        )

        # Initialize model and sub-models
        results = self._initialize_submodels(
            model, plugin_initializer_args, copy_inlet_state, sub_initializers, _log
        )

        # Solve full model including plug-ins
        results = self._solve_full_model(model, _log, results)

        # Clean up plug-ins
        self._cleanup(model, plugin_initializer_args, sub_initializers, _log)

        return results

    def _prepare_plugins(self, model, plugin_initializer_args, logger):
        sub_initializers = {}
        plugin_initializer_args = dict(plugin_initializer_args)
        for sm in model.initialization_order:
            if sm is not model:
                # Get initializers for plug-ins
                sub_initializers[sm] = self.get_submodel_initializer(sm)()

                if sm not in plugin_initializer_args.keys():
                    plugin_initializer_args[sm] = {}

                # Call prepare method for plug-ins
                sub_initializers[sm].plugin_prepare(sm, **plugin_initializer_args[sm])

        logger.info_high("Step 1: preparation complete.")

        return sub_initializers, plugin_initializer_args

    def _initialize_submodels(
        self, model, plugin_initializer_args, copy_inlet_state, sub_initializers, logger
    ):
        results = None

        for sm in model.initialization_order:
            if sm is model:
                results = self._initialize_main_model(model, copy_inlet_state, logger)
            else:
                sub_initializers[sm].plugin_initialize(
                    sm, **plugin_initializer_args[sm]
                )

        if results is None:
            raise InitializationError(
                f"Main model ({model.name}) was not initialized (no results returned). "
                f"This is likely due to an error in the model.initialization_order."
            )

        logger.info_high("Step 2: sub-model initialization complete.")

        return results

    def _initialize_main_model(self, model, copy_inlet_state, logger):
        # Initialize properties
        try:
            # Guess a 0-D control volume
            self._init_props_0D(model, copy_inlet_state)
        except AttributeError:
            # Assume it must be a 1-D control volume
            self._init_props_1D(model)

        logger.info_high("Step 2a: properties initialization complete")

        # Initialize reactions if they exist
        if hasattr(model.control_volume, "reactions"):
            self._init_rxns(model)
        logger.info_high("Step 2b: reactions initialization complete")

        # Solve main model
        solve_log = idaeslog.getSolveLogger(
            model.name, self.config.output_level, tag="unit"
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = self._get_solver().solve(model, tee=slc.tee)

        logger.info_high(f"Step 2c: model {idaeslog.condition(results)}.")

        return results

    def _init_props_0D(self, model, copy_inlet_state):
        # Initialize inlet properties - inlet state should already be fixed
        prop_init = self.get_submodel_initializer(model.control_volume.properties_in)

        if prop_init is not None:
            prop_init.initialize(
                model.control_volume.properties_in,
                solver=self.config.solver,
                optarg=self.config.solver_options,
                outlvl=self.config.output_level,
            )

        if not copy_inlet_state:
            # Just in case the user set a different initializer for the outlet
            prop_init = self.get_submodel_initializer(
                model.control_volume.properties_out
            )

            if prop_init is not None:
                prop_init.initialize(
                    model.control_volume.properties_out,
                    solver=self.config.solver,
                    optarg=self.config.solver_options,
                    outlvl=self.config.output_level,
                )
        else:
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

    def _init_props_1D(self, model):
        prop_init = self.get_submodel_initializer(model.control_volume.properties)

        if prop_init is not None:
            prop_init.initialize(
                model.control_volume.properties,
                solver=self.config.solver,
                optarg=self.config.solver_options,
                outlvl=self.config.output_level,
            )

    def _init_rxns(self, model):
        rxn_init = self.get_submodel_initializer(model.control_volume.reactions)

        if rxn_init is not None:
            rxn_init.initialize(
                model.control_volume.reactions,
                solver=self.config.solver,
                optarg=self.config.solver_options,
                outlvl=self.config.output_level,
            )

    def _solve_full_model(self, model, logger, results):
        # Check to see if there are any plug-ins
        if len(model.initialization_order) > 1:
            # Solve model with plugins
            solve_log = idaeslog.getSolveLogger(
                model.name, self.config.output_level, tag="unit"
            )
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                results = self._get_solver().solve(model, tee=slc.tee)

        logger.info_high(
            f"Step 3: full model initialization {idaeslog.condition(results)}."
        )

        return results

    def _cleanup(self, model, plugin_initializer_args, sub_initializers, logger):
        for sm in reversed(model.initialization_order):
            if sm is not model:
                sub_initializers[sm].plugin_finalize(sm, **plugin_initializer_args[sm])

        logger.info_high("Step 4: clean up completed.")

    def _get_solver(self):
        if self._solver is None:
            self._solver = get_solver(self.config.solver, self.config.solver_options)

        return self._solver
