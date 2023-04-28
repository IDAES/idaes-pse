#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Initializer class for implementing Hierarchical initialization routines for
IDAES models with standard forms (e.g. units with 1 control volume)
"""
from pyomo.environ import Block
from pyomo.common.config import ConfigValue, Bool

from idaes.core.initialization.initializer_base import ModularInitializerBase
from idaes.core.util import to_json, from_json, StoreSpec
import idaes.logger as idaeslog
from idaes.core.util.model_statistics import variables_in_activated_constraints_set

__author__ = "Andrew Lee"


class SingleControlVolumeUnitInitializer(ModularInitializerBase):
    """
    This is a general purpose Initializer object for unit models which
    have a single Control Volume named 'control_volume'.

    For details of the initialization routine, please see the documentation.

    """

    CONFIG = ModularInitializerBase.CONFIG()
    CONFIG.declare(
        "always_estimate_states",
        ConfigValue(
            default=False,
            domain=Bool,
            doc="Whether initialization routine should estimate values for "
            "state variables that already have values. Note that if set to True, this will "
            "overwrite any initial guesses provided.",
        ),
    )

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
            copy_inlet_state: bool (default=False). Whether to copy inlet state to other states or not
                (0-D control volumes only). Copying will generally be faster, but inlet states may not contain
                all properties required elsewhere.

        Returns:
            Pyomo solver results object
        """
        # The default initialization_routine is sufficient
        return super().initialization_routine(
            model=model,
            plugin_initializer_args=plugin_initializer_args,
            copy_inlet_state=copy_inlet_state,
        )

    def initialize_main_model(self, model: Block, copy_inlet_state: bool = False):
        """
        Initialization routine for models with a single control volume.

        Args:
            model: current model being initialized
            copy_inlet_state: bool (default=False). Whether to copy inlet state to other states or not
                (0-D control volumes only). Copying will generally be faster, but inlet states may not contain
                all properties required elsewhere.

        Returns:
            Pyomo solver results object from solve of main model

        """
        # Check to make sure model has something named "control_volume"
        if not hasattr(model, "control_volume"):
            raise TypeError(
                f"Model {model.name} does not appear to be a standard form unit model. "
                f"Please use an Initializer specific to the model being initialized."
            )

        # Get logger
        _log = self.get_logger(model)

        self.initialize_control_volume(model.control_volume, copy_inlet_state)

        # Solve main model
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = self._get_solver().solve(model, tee=slc.tee)

        _log.info_high(f"Control volume initialization {idaeslog.condition(results)}.")

        return results

    def initialize_control_volume(
        self, control_volume: Block, copy_inlet_state: bool = False
    ):
        """
        Initialization routine for control volumes (supports both 0D and 1D control volumes).

        This routine initialized the state and reaction blocks associated with the control volume.

        Args:
            control_volume: control_volume to be initialized
            copy_inlet_state: bool (default=False). Whether to copy inlet state to other states or not
                (0-D control volumes only). Copying will generally be faster, but inlet states may not contain
                all properties required elsewhere.

        Returns:
            None

        """
        # Get logger
        _log = self.get_logger(control_volume)

        # Initialize properties
        if hasattr(control_volume, "properties_in"):
            # 0-D control volume
            self._init_props_0D(control_volume, copy_inlet_state)
        else:
            # 1-D control volume
            self._init_props_1D(control_volume)

        _log.info_high("Control volume properties initialization complete")

        # Initialize reactions if they exist
        if hasattr(control_volume, "reactions"):
            self._init_rxns(control_volume)
        _log.info_high("Control volume reactions initialization complete")

    def _init_props_0D(self, control_volume, copy_inlet_state):
        # Initialize inlet properties - inlet state should already be fixed
        prop_init = self.get_submodel_initializer(control_volume.properties_in)

        if prop_init is not None:
            prop_init.initialize(
                model=control_volume.properties_in,
                output_level=self.get_output_level(),
            )

        if not copy_inlet_state:
            # Just in case the user set a different initializer for the outlet
            prop_init = self.get_submodel_initializer(control_volume.properties_out)

            # Estimate missing values for outlet state
            control_volume.estimate_outlet_state(
                always_estimate=self.config.always_estimate_states
            )

            if prop_init is not None:
                prop_init.initialize(
                    model=control_volume.properties_out,
                    output_level=self.get_output_level(),
                )
        else:
            # Map solution from inlet properties to outlet properties
            state = to_json(
                control_volume.properties_in,
                wts=StoreSpec().value(),
                return_dict=True,
            )
            from_json(
                control_volume.properties_out,
                sd=state,
                wts=StoreSpec().value(only_not_fixed=True),
            )

    def _init_props_1D(self, control_volume):
        prop_init = self.get_submodel_initializer(control_volume.properties)

        # Estimate missing values for states
        control_volume.estimate_states(
            always_estimate=self.config.always_estimate_states
        )

        if prop_init is not None:
            prop_init.initialize(
                control_volume.properties,
                output_level=self.get_output_level(),
            )

    def _init_rxns(self, control_volume):
        rxn_init = self.get_submodel_initializer(control_volume.reactions)

        if rxn_init is not None:
            # Reaction blocks depend on vars from other blocks
            # Need to make sure these are all fixed before initializing
            fixed_vars = []
            vset = variables_in_activated_constraints_set(control_volume.reactions)
            for v in vset:
                if v.parent_block().parent_component() is not control_volume.reactions:
                    # Variable external to reactions
                    if not v.fixed:
                        v.fix()
                        fixed_vars.append(v)

            rxn_init.initialize(
                control_volume.reactions,
                output_level=self.get_output_level(),
            )

            # Unfix any vars we fixed previously
            for v in fixed_vars:
                v.unfix()
