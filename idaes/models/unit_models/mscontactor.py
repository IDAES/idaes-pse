#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
IDAES model for a generic multiple-stream contactor unit.
"""
from functools import partial

from pandas import DataFrame

# Import Pyomo libraries
from pyomo.environ import (
    Block,
    ComponentMap,
    Constraint,
    Expression,
    RangeSet,
    Reals,
    Set,
    units,
    Var,
)
from pyomo.common.config import ConfigDict, ConfigValue, Bool, In
from pyomo.contrib.incidence_analysis import solve_strongly_connected_components
from pyomo.dae import DerivativeVar

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    FlowDirection,
    MaterialFlowBasis,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
    BurntToast,
    PropertyNotSupportedError,
)
from idaes.core.initialization import ModularInitializerBase
from idaes.core.initialization.initializer_base import StoreState
from idaes.core.util.model_serializer import to_json, from_json
import idaes.logger as idaeslog
from idaes.core.util.units_of_measurement import report_quantity
from idaes.core.scaling import (
    ConstraintScalingScheme,
    CustomScalerBase,
    DefaultScalingRecommendation,
)

__author__ = "Andrew Lee, Douglas Allan"


class MSContactorScaler(CustomScalerBase):
    """
    Scaler for the MSContactor.
    """

    DEFAULT_SCALING_FACTORS = {
        # We could scale volume by magnitude if it were being fixed
        # by the user, but we often have the volume given by an
        # equality constraint involving geometry in the parent
        # unit model.
        "volume": DefaultScalingRecommendation.userInputRequired,
        "volume_frac_stream": 2,
        # Phase fraction may have already been created by the property package.
        # Also, no "phase_fraction" Var exists on the MSContactor, only
        # "stream_name_phase_fraction" Vars. Because we do not know the stream
        # names a priori, this value is copied for those streams as part of
        # the variable_scaling_routine (if the user hasn't already set a value)
        "phase_fraction": 10,
    }

    def _scale_stream_reaction_variables(
        self, model: Block, stream: str, reaction_type: str, overwrite: bool
    ):
        """
        Scaling method for reactions within a single stream (rate reactions, equilibrium reactions
        and inherent reactions). Scales the reaction generation term based on the sensitivity of
        the material balance for that particular phase component pair, then scales the extent of
        reaction based on the the most sensitive generation term.

        Args:
            stream_state: State block for the stream whose reaction variables are being scaled
            rxn_gen: Reaction
        """
        assert reaction_type in {"rate", "equilibrium", "inherent"}
        stream_state = getattr(model, stream)
        rxn_gen = getattr(model, stream + f"_{reaction_type}_reaction_generation")
        rxn_extent = getattr(model, stream + f"_{reaction_type}_reaction_extent")

        # Get stoichiometry dictionary
        if reaction_type == "inherent":
            idx0 = stream_state.index_set().first()
            rxn_index = stream_state.params.inherent_reaction_idx
            stoich = stream_state[idx0].params.inherent_reaction_stoichiometry
        else:
            reaction_block = getattr(model, stream + "_reactions")
            idx0 = reaction_block.index_set().first()
            if reaction_type == "rate":
                rxn_index = reaction_block[idx0].params.rate_reaction_idx
                stoich = reaction_block[idx0].params.rate_reaction_stoichiometry
            else:
                rxn_index = reaction_block[idx0].params.equilibrium_reaction_idx
                stoich = reaction_block[idx0].params.equilibrium_reaction_stoichiometry

        pc_set = stream_state.phase_component_set
        for t, e in stream_state:
            for p, j in pc_set:
                # Scale the reaction generation term to be
                # the same magnitude as the material flow term
                nom = self.get_expression_nominal_value(
                    stream_state[t, e].get_material_flow_terms(p, j)
                )
                self.set_variable_scaling_factor(
                    rxn_gen[t, e, p, j], 1 / nom, overwrite=overwrite
                )
            # We want to scale the rate reaction extent based on the
            # component in that reaction that is tracked with the *greatest*
            # precision. Therefore iterate over components and look for the
            # one with the largest scaling factor (or smallest nominal value)
            # when weighted by stoichiometric coefficient.
            rxn_dict = {}
            for idx, coeff in stoich.items():
                rxn, p, j = idx
                if (p, j) in stream_state.phase_component_set:
                    # Scale the reaction generation term to be
                    # the same magnitude as the material flow term
                    sf = self.get_scaling_factor(rxn_gen[t, e, p, j])
                    if coeff != 0:
                        if rxn in rxn_dict:
                            rxn_dict[rxn] = min(1 / (abs(coeff) * sf), rxn_dict[rxn])
                        else:
                            rxn_dict[rxn] = 1 / (abs(coeff) * sf)
            # Now that we've collected the sensitivities for each reaction
            # we can scale the reaction extent.
            for rxn in rxn_index:
                if rxn not in rxn_dict:
                    raise ConfigurationError(
                        f"Reaction {rxn} does not have any nonzero stoichiometric coefficient."
                    )
                self.set_variable_scaling_factor(
                    rxn_extent[t, e, rxn], 1 / rxn_dict[rxn], overwrite=overwrite
                )

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        """
        Variable scaling routine for the MSContactor.
        Args:
            model: instance of MSContactor to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scaler objects to use for sub-models

        Returns:
            None
        """
        if submodel_scalers is None:
            submodel_scalers = {}

        # Step 1: Property and reaction scaling

        # Step 1a: Propagate any existing scaling from inlet to outlet
        for stream in model.streams:
            # If the stream has an inlet state, we'll propagate scaling
            # factors from that to the stream states and any side stream
            # state that may exist. If no inlet state exists, it is not
            # obvious what, if any, state scaling propagation should occur.
            # In that case, the user can set scaling factors themselves.
            if hasattr(model, stream + "_inlet_state"):
                feed = getattr(model, stream + "_inlet_state")
                side_stream = None
                ss_set = None
                if hasattr(model, stream + "_side_stream_state"):
                    side_stream = getattr(model, stream + "_side_stream_state")
                    ss_set = getattr(model, stream + "_side_stream_set")
                stream_state = getattr(model, stream)
                for t in model.flowsheet().time:
                    for e in model.elements:
                        self.propagate_state_scaling(
                            target_state=stream_state[t, e],
                            source_state=feed[t],
                            overwrite=overwrite,
                        )
                    if side_stream is not None:
                        for s in ss_set:
                            self.propagate_state_scaling(
                                target_state=side_stream[t, s],
                                source_state=feed[t],
                                overwrite=overwrite,
                            )

        # Step 1b: Call Scalers for state blocks and stream-specific reaction blocks
        for stream in model.streams:
            if hasattr(model, stream + "_inlet_state"):
                self.call_submodel_scaler_method(
                    submodel=getattr(model, stream + "_inlet_state"),
                    submodel_scalers=submodel_scalers,
                    method="variable_scaling_routine",
                    overwrite=overwrite,
                )
            if hasattr(model, stream + "_side_stream_state"):
                self.call_submodel_scaler_method(
                    submodel=getattr(model, stream + "_side_stream_state"),
                    submodel_scalers=submodel_scalers,
                    method="variable_scaling_routine",
                    overwrite=overwrite,
                )
            self.call_submodel_scaler_method(
                submodel=getattr(model, stream),
                submodel_scalers=submodel_scalers,
                method="variable_scaling_routine",
                overwrite=overwrite,
            )
            if hasattr(model, stream + "_reactions"):
                self.call_submodel_scaler_method(
                    submodel=getattr(model, stream + "_reactions"),
                    submodel_scalers=submodel_scalers,
                    method="variable_scaling_routine",
                    overwrite=overwrite,
                )
        # Step 1c: Call scaler for heterogeneous reaction block
        if hasattr(model, "heterogeneous_reactions"):
            self.call_submodel_scaler_method(
                submodel=model.heterogeneous_reactions,
                method="variable_scaling_routine",
                submodel_scalers=submodel_scalers,
                overwrite=overwrite,
            )

        # Step 2: Scaling MSContactor level variables
        # Step 2a: Scale holdup variables
        if hasattr(model, "volume"):
            for var in model.volume.values():
                self.scale_variable_by_default(var, overwrite=overwrite)
        if hasattr(model, "volume_frac_stream"):
            for var in model.volume_frac_stream.values():
                self.scale_variable_by_default(var, overwrite=overwrite)
        for stream in model.config.streams.keys():
            # Phase fraction
            if hasattr(model, stream + "_phase_fraction"):
                # Copy default scaling factor for phase_fraction to use
                # for stream_phase_fraction.
                if stream + "_phase_fraction" not in self.default_scaling_factors:
                    self.default_scaling_factors[stream + "_phase_fraction"] = (
                        self.default_scaling_factors["phase_fraction"]
                    )
                phase_frac = getattr(model, stream + "_phase_fraction")
                for var in phase_frac.values():
                    self.scale_variable_by_default(var, overwrite=overwrite)

            # Material holdup
            if hasattr(model, stream + "_material_holdup"):
                material_holdup = getattr(model, stream + "_material_holdup")
                material_holdup_con = getattr(
                    model, stream + "_material_holdup_constraint"
                )
                for idx in material_holdup:
                    self.scale_variable_by_definition_constraint(
                        material_holdup[idx],
                        material_holdup_con[idx],
                        overwrite=overwrite,
                    )

            # Energy Holdup
            # TODO scale_variable_by_definition_constraint behaves badly
            # when internal energy passes through zero
            if hasattr(model, stream + "_energy_holdup"):
                energy_holdup = getattr(model, stream + "_energy_holdup")
                energy_holdup_con = getattr(model, stream + "_energy_holdup_constraint")

                for idx in energy_holdup:
                    self.scale_variable_by_definition_constraint(
                        energy_holdup[idx], energy_holdup_con[idx], overwrite=overwrite
                    )

        # Step 2b: Scale stream interaction terms
        if hasattr(model, "material_transfer_term"):
            for idx in model.stream_component_interactions:
                stream1 = idx[0]
                stream2 = idx[1]
                comp = idx[2]
                stream_state1 = getattr(model, stream1)
                stream_state2 = getattr(model, stream2)

                for t in model.flowsheet().time:
                    for e in model.elements:
                        # Determine the precision at which each stream
                        # is tracking comp in its material balances
                        nom1 = self.get_expression_nominal_value(
                            sum(
                                stream_state1[t, e].get_material_flow_terms(p, comp)
                                for p in stream_state1[t, e].phase_list
                                if (p, comp) in stream_state1[t, e].phase_component_set
                            )
                        )
                        nom2 = self.get_expression_nominal_value(
                            sum(
                                stream_state2[t, e].get_material_flow_terms(p, comp)
                                for p in stream_state2[t, e].phase_list
                                if (p, comp) in stream_state2[t, e].phase_component_set
                            )
                        )
                        # Scale the material transfer term at the same
                        # scale as the stream that's tracking comp with
                        # the *greatest* precision (and therefore the
                        # smallest nominal value)
                        sf = 1 / min(nom1, nom2)
                        self.set_variable_scaling_factor(
                            model.material_transfer_term[t, e, idx],
                            sf,
                            overwrite=overwrite,
                        )
        # Step 2c: Scale heterogeneous reactions
        if hasattr(model, "heterogeneous_reactions"):
            for stream in model.config.streams.keys():
                stream_state = getattr(model, stream)
                pc_set = stream_state.phase_component_set

                het_rxn_gen = getattr(
                    model, stream + "_heterogeneous_reactions_generation"
                )
                for t in model.flowsheet().time:
                    for e in model.elements:
                        for p, j in pc_set:
                            # Scale the heterogeneous reaction term to be
                            # the same magnitude as the material flow term
                            nom = self.get_expression_nominal_value(
                                stream_state[t, e].get_material_flow_terms(p, j)
                            )
                            self.set_variable_scaling_factor(
                                het_rxn_gen[t, e, p, j], 1 / nom, overwrite=overwrite
                            )

            # Get the stoichiometry for the heterogeneous reactions
            t0 = model.flowsheet().time.first()
            e0 = model.elements.first()
            stoich = model.heterogeneous_reactions[t0, e0].params.reaction_stoichiometry

            # Scale heterogeneous reaction extent
            for t in model.flowsheet().time:
                for e in model.elements:
                    # We want to scale the heterogeneous reaction extent based on the
                    # component in that reaction that is tracked with the *greatest*
                    # precision. Therefore iterate over components and look for the
                    # one with the largest scaling factor (or smallest nominal value)
                    # when weighted by stoichiometric coefficient.
                    rxn_dict = {}
                    for idx, coeff in stoich.items():
                        rxn, p, j = idx
                        for stream in model.config.streams.keys():
                            stream_state = getattr(model, stream)[t, e]
                            het_rxn_gen = getattr(
                                model, stream + "_heterogeneous_reactions_generation"
                            )
                            if (p, j) in stream_state.phase_component_set:
                                sf = self.get_scaling_factor(het_rxn_gen[t, e, p, j])
                                if coeff != 0:
                                    if rxn in rxn_dict:
                                        rxn_dict[rxn] = min(
                                            1 / (abs(coeff) * sf), rxn_dict[rxn]
                                        )
                                    else:
                                        rxn_dict[rxn] = 1 / (abs(coeff) * sf)
                    # Now that we've collected the sensitivities for each reaction
                    # we can scale the reaction extent.
                    for rxn in model.config.heterogeneous_reactions.reaction_idx:
                        if rxn not in rxn_dict:
                            raise ConfigurationError(
                                f"Heterogeneous reaction {rxn} does not have any nonzero stoichiometric coefficient."
                            )
                        self.set_variable_scaling_factor(
                            model.heterogeneous_reaction_extent[t, e, rxn],
                            1 / rxn_dict[rxn],
                            overwrite=overwrite,
                        )
        # Step 2d: Scale single-stream reactions
        for stream in model.config.streams.keys():
            stream_state = getattr(model, stream)
            if hasattr(model, stream + "_rate_reaction_extent"):
                self._scale_stream_reaction_variables(
                    model=model,
                    stream=stream,
                    reaction_type="rate",
                    overwrite=overwrite,
                )
            if hasattr(model, stream + "_equilibrium_reaction_extent"):
                self._scale_stream_reaction_variables(
                    model=model,
                    stream=stream,
                    reaction_type="equilibrium",
                    overwrite=overwrite,
                )
            if hasattr(model, stream + "_inherent_reaction_extent"):
                self._scale_stream_reaction_variables(
                    model=model,
                    stream=stream,
                    reaction_type="inherent",
                    overwrite=overwrite,
                )
        # Step 3a: Scale energy transfer between streams
        if hasattr(model, "energy_transfer_term"):
            for idx in model.stream_interactions:
                stream1 = idx[0]
                stream2 = idx[1]
                stream_state1 = getattr(model, stream1)
                stream_state2 = getattr(model, stream2)

                for t in model.flowsheet().time:
                    for e in model.elements:
                        if model.config.streams[stream1].has_energy_balance:
                            nom1 = self.get_expression_nominal_value(
                                sum(
                                    stream_state1[t, e].get_enthalpy_flow_terms(p)
                                    for p in stream_state1[t, e].phase_list
                                )
                            )
                        else:
                            nom1 = None
                        if model.config.streams[stream2].has_energy_balance:
                            nom2 = self.get_expression_nominal_value(
                                sum(
                                    stream_state2[t, e].get_enthalpy_flow_terms(p)
                                    for p in stream_state2[t, e].phase_list
                                )
                            )
                        else:
                            nom2 = None
                        if nom1 is None and nom2 is None:
                            raise BurntToast(
                                "Energy transfer term should not be constructed if "
                                "neither stream has an energy balance, please report "
                                "this problem to the IDAES developers."
                            )
                        elif nom1 is None:
                            sf = 1 / nom2
                        elif nom2 is None:
                            sf = 1 / nom1
                        else:
                            sf = 1 / min(nom1, nom2)
                        self.set_variable_scaling_factor(
                            model.energy_transfer_term[t, e, idx],
                            sf,
                            overwrite=overwrite,
                        )

        # Step 3b: Scale heat transfer between stream and environment
        for stream in model.config.streams.keys():
            if hasattr(model, stream + "_heat"):
                heat_var = getattr(model, stream + "_heat")
                for t in model.flowsheet().time:
                    for e in model.elements:
                        nom = self.get_expression_nominal_value(
                            sum(
                                stream_state[t, e].get_enthalpy_flow_terms(p)
                                for p in stream_state[t, e].phase_list
                            )
                        )
                        self.set_variable_scaling_factor(
                            heat_var[t, e],
                            1 / nom,
                            overwrite=overwrite,
                        )
        # Step 4a: Scale pressure difference
        for stream in model.config.streams.keys():
            if hasattr(model, stream + "_deltaP"):
                stream_state = getattr(model, stream)
                deltaP = getattr(model, stream + "_deltaP")
                for idx, vardata in deltaP.items():
                    self.scale_variable_by_component(
                        vardata, stream_state[idx].pressure, overwrite=overwrite
                    )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to constraints in model.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            None
        """
        # Step 1a: Call Scalers for state blocks and stream reaction blocks
        for stream in model.streams:
            if hasattr(model, stream + "_inlet_state"):
                self.call_submodel_scaler_method(
                    submodel=getattr(model, stream + "_inlet_state"),
                    submodel_scalers=submodel_scalers,
                    method="constraint_scaling_routine",
                    overwrite=overwrite,
                )
            if hasattr(model, stream + "_side_stream_state"):
                self.call_submodel_scaler_method(
                    submodel=getattr(model, stream + "_side_stream_state"),
                    submodel_scalers=submodel_scalers,
                    method="constraint_scaling_routine",
                    overwrite=overwrite,
                )
            self.call_submodel_scaler_method(
                submodel=getattr(model, stream),
                submodel_scalers=submodel_scalers,
                method="constraint_scaling_routine",
                overwrite=overwrite,
            )
            if hasattr(model, stream + "_reactions"):
                self.call_submodel_scaler_method(
                    submodel=getattr(model, stream + "_reactions"),
                    submodel_scalers=submodel_scalers,
                    method="constraint_scaling_routine",
                    overwrite=overwrite,
                )

        # Step 1b: Call Scalers for heterogeneous blocks
        if hasattr(model, "heterogeneous_reactions"):
            self.call_submodel_scaler_method(
                submodel=model.heterogeneous_reactions,
                method="constraint_scaling_routine",
                submodel_scalers=submodel_scalers,
                overwrite=overwrite,
            )

        # Step 2: Scaling MSContactor level variables
        for stream in model.config.streams.keys():
            stream_state = getattr(model, stream)

            # Step 2a: Holdup equations
            if hasattr(model, "sum_volume_frac"):
                sum_volume_frac_eqn = getattr(model, "sum_volume_frac")
                # Register the fact this constraint is well-scaled by default
                for condata in sum_volume_frac_eqn.values():
                    self.set_component_scaling_factor(condata, 1, overwrite=overwrite)

            if hasattr(model, stream + "_sum_phase_fractions"):
                sum_phase_frac_eqn = getattr(model, stream + "_sum_phase_fractions")
                for condata in sum_phase_frac_eqn.values():
                    # Register the fact this constraint is well-scaled by default
                    self.set_component_scaling_factor(condata, 1, overwrite=overwrite)
            if hasattr(model, stream + "_material_holdup_constraint"):
                holdup = getattr(model, stream + "_material_holdup")
                holdup_eqn = getattr(model, stream + "_material_holdup_constraint")
                for idx in holdup_eqn:
                    self.scale_constraint_by_component(
                        holdup_eqn[idx], holdup[idx], overwrite=overwrite
                    )
            if hasattr(model, stream + "_energy_holdup_constraint"):
                energy_holdup = getattr(model, stream + "_energy_holdup")
                energy_holdup_eqn = getattr(model, stream + "_energy_holdup_constraint")
                for idx in energy_holdup_eqn:
                    self.scale_constraint_by_component(
                        energy_holdup_eqn[idx],
                        energy_holdup[idx],
                        overwrite=overwrite,
                    )
            # Step 2b: Single stream reaction equations
            if hasattr(model, stream + "_rate_reaction_constraint"):
                rate_rxn_gen = getattr(model, stream + "_rate_reaction_generation")
                rate_rxn_con = getattr(model, stream + "_rate_reaction_constraint")
                for idx, con in rate_rxn_con.items():
                    self.scale_constraint_by_component(
                        con, rate_rxn_gen[idx], overwrite=overwrite
                    )

            if hasattr(model, stream + "_equilibrium_reaction_constraint"):
                equil_rxn_gen = getattr(
                    model, stream + "_equilibrium_reaction_generation"
                )
                equil_rxn_con = getattr(
                    model, stream + "_equilibrium_reaction_constraint"
                )
                for idx, con in equil_rxn_con.items():
                    self.scale_constraint_by_component(
                        con, equil_rxn_gen[idx], overwrite=overwrite
                    )

            if hasattr(model, stream + "_inherent_reaction_constraint"):
                inherent_reaction_generation = getattr(
                    model, stream + "_inherent_reaction_generation"
                )
                inherent_reaction_eqn = getattr(
                    model, stream + "_inherent_reaction_constraint"
                )
                for idx in inherent_reaction_eqn:
                    self.scale_constraint_by_component(
                        inherent_reaction_eqn[idx],
                        inherent_reaction_generation[idx],
                        overwrite=overwrite,
                    )
            # Step 2c: Heterogeneous reaction equations
            if hasattr(model, "heterogeneous_reactions"):
                heterogeneous_reaction_generation = getattr(
                    model, stream + "_heterogeneous_reactions_generation"
                )
                heterogeneous_reaction_eqn = getattr(
                    model, stream + "_heterogeneous_reaction_constraint"
                )
                for idx in heterogeneous_reaction_eqn:
                    self.scale_constraint_by_component(
                        heterogeneous_reaction_eqn[idx],
                        heterogeneous_reaction_generation[idx],
                        overwrite=overwrite,
                    )
            # Step 2d: Material balance
            if hasattr(model, stream + "_material_balance"):
                mbal = getattr(model, stream + "_material_balance")
                for t, e, comp in mbal:
                    nom = self.get_expression_nominal_value(
                        sum(
                            stream_state[t, e].get_material_flow_terms(p, comp)
                            for p in stream_state[t, e].phase_list
                            if (p, comp) in stream_state[t, e].phase_component_set
                        )
                    )
                    self.set_component_scaling_factor(
                        mbal[t, e, comp],
                        1 / nom,
                        overwrite=overwrite,
                    )
            # Step 2e: Energy balance
            if hasattr(model, stream + "_energy_balance"):
                ebal = getattr(model, stream + "_energy_balance")
                for t, e in ebal:
                    nom = self.get_expression_nominal_value(
                        sum(
                            stream_state[t, e].get_enthalpy_flow_terms(p)
                            for p in stream_state[t, e].phase_list
                        )
                    )
                    self.set_component_scaling_factor(
                        ebal[t, e],
                        1 / nom,
                        overwrite=overwrite,
                    )
            # Step 2f: Pressure balance
            if hasattr(model, stream + "_pressure_balance"):
                pbal = getattr(model, stream + "_pressure_balance")
                for condata in pbal.values():
                    self.scale_constraint_by_nominal_value(
                        condata,
                        scheme=ConstraintScalingScheme.inverseMaximum,
                        overwrite=overwrite,
                    )
            if hasattr(model, stream + "_side_stream_pressure_balance"):
                spbal = getattr(model, stream + "_side_stream_pressure_balance")
                for condata in spbal.values():
                    self.scale_constraint_by_nominal_value(
                        condata,
                        scheme=ConstraintScalingScheme.inverseMaximum,
                        overwrite=overwrite,
                    )


class MSContactorInitializer(ModularInitializerBase):
    """
    This is a general purpose sequential-modular Initializer object for
    multi-stream contactor unit models.

    This routine starts by deactivating any constraints that are not
    part of the base model and fixing all inter-stream transfer variables.
    The model is then solved using the Pyomo ssc_solver function
    to initialize each stream separately.

    The inter-stream transfer variables are then unfixed and the additional
    constraints reactivated, and the full model solved using the user-specified
    solver.

    """

    CONFIG = ModularInitializerBase.CONFIG()

    CONFIG.declare(
        "ssc_solver_options",
        ConfigDict(
            implicit=True,
            description="Dict of arguments for solver calls by ssc_solver",
        ),
    )
    CONFIG.declare(
        "calculate_variable_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to 1x1 block solver",
            doc="Dict of options to pass to calc_var_kwds argument in "
            "scc_solver method.",
        ),
    )

    def initialization_routine(
        self,
        model: Block,
    ):
        """
        Initialization routine for MSContactor Blocks.

        Args:
            model: model to be initialized

        Returns:
            None
        """
        # Get loggers
        init_log = idaeslog.getInitLogger(
            model.name, self.get_output_level(), tag="unit"
        )
        solve_log = idaeslog.getSolveLogger(
            model.name, self.get_output_level(), tag="unit"
        )

        # Get current model state
        initial_state = to_json(model, wts=StoreState, return_dict=True)

        # Isolate streams by fixing inter-stream variables
        model.material_transfer_term.fix()

        # Deactivate any constraints that are not part of the base model
        # First, build list of names for known constraints
        const_names = []
        for s in model.config.streams.keys():
            const_names.append(model.name + "." + s + "_rate_reaction_constraint")
            const_names.append(
                model.name + "." + s + "_equilibrium_reaction_constraint"
            )
            const_names.append(model.name + "." + s + "_inherent_reaction_constraint")
            const_names.append(
                model.name + "." + s + "_heterogeneous_reaction_constraint"
            )
            const_names.append(model.name + "." + s + "_material_balance")
            const_names.append(model.name + "." + s + "_energy_balance")
            const_names.append(model.name + "." + s + "_pressure_balance")
            const_names.append(model.name + "." + s + "_side_stream_pressure_balance")

            try:
                # If has rate reactions ,fi extent to 0 for first pass
                getattr(model, s + "_rate_reaction_extent").fix(0)
            except AttributeError:
                pass

        # Fix extents for heterogeneous reactions to 0 for first pass if present
        if hasattr(model, "heterogeneous_reaction_extent"):
            model.heterogeneous_reaction_extent.fix(0)

        # Iterate through all constraints attached to model - do not search sub-blocks
        for c in model.component_objects(Constraint, descend_into=False):
            # Deactivate constraint if its name is not in list of known names
            if c.name not in const_names:
                c.deactivate()

        # Call css_solver
        solver = self._get_solver()
        solve_strongly_connected_components(
            model,
            solver=solver,
            solve_kwds=self.config.ssc_solver_options,
            calc_var_kwds=self.config.calculate_variable_options,
        )
        init_log.info("Stream Initialization Completed.")

        # Revert state
        from_json(model, sd=initial_state, wts=StoreState)

        # Solve full model
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(model, tee=slc.tee)
        init_log.info(f"Initialization Completed, {idaeslog.condition(res)}")

        return res


# Config dict to contain information on each stream
STREAM_CONFIG = ConfigDict()
STREAM_CONFIG.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for given stream",
        doc="""Property parameter object used to define property calculations for given stream,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
    ),
)
STREAM_CONFIG.declare(
    "property_package_args",
    ConfigDict(
        implicit=True,
        description="Dict of arguments to use for constructing property package",
        doc="""A ConfigDict with arguments to be passed to property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)
STREAM_CONFIG.declare(
    "reaction_package",
    ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for stream",
        doc="""Reaction parameter object used to define reaction calculations
for stream. **default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
    ),
)
STREAM_CONFIG.declare(
    "reaction_package_args",
    ConfigDict(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
    ),
)
STREAM_CONFIG.declare(
    "flow_direction",
    ConfigValue(
        default=FlowDirection.forward,
        domain=In(FlowDirection),
        doc="Direction of flow for stream",
        description="FlowDirection Enum indicating direction of "
        "flow for given stream. Default=FlowDirection.forward.",
    ),
)
STREAM_CONFIG.declare(
    "has_feed",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether stream has a feed.",
        description="Bool indicating whether stream has a feed Port and inlet "
        "state, or if all flow is provided via mass transfer. Default=True.",
    ),
)
STREAM_CONFIG.declare(
    "has_rate_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether rate-based reactions occur in stream.",
    ),
)
STREAM_CONFIG.declare(
    "has_equilibrium_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether equilibrium-based reactions occur in stream.",
    ),
)
STREAM_CONFIG.declare(
    "has_energy_balance",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether to include energy balance for stream. Default=True.",
    ),
)
STREAM_CONFIG.declare(
    "has_heat_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include external heat transfer terms in energy "
        "balance for stream. Default=False.",
    ),
)
STREAM_CONFIG.declare(
    "has_heat_of_reaction",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether heat of reaction terms should be included in energy "
        "balance for stream (required reactions). Default=False.",
    ),
)
STREAM_CONFIG.declare(
    "has_pressure_balance",
    ConfigValue(
        default=True,
        domain=Bool,
        doc="Bool indicating whether to include pressure balance for stream. Default=True.",
    ),
)
STREAM_CONFIG.declare(
    "has_pressure_change",
    ConfigValue(
        default=False,
        domain=Bool,
        doc="Bool indicating whether to include deltaP terms in pressure balance for "
        "stream. Default=False.",
    ),
)
STREAM_CONFIG.declare(
    "side_streams",
    ConfigValue(
        default=None,
        domain=list,
        doc="List of finite elements at which a side stream should be included.",
    ),
)


@declare_process_block_class("MSContactor")
class MSContactorData(UnitModelBlockData):
    """
    Multi-Stream Contactor Unit Model Class
    """

    default_initializer = MSContactorInitializer
    default_scaler = MSContactorScaler

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "streams",
        ConfigDict(
            implicit=True,
            implicit_domain=STREAM_CONFIG,
            description="Dict of streams and associated property packages",
            doc="ConfigDict with keys indicating names for each stream in system and "
            "values indicating property package and associated arguments.",
        ),
    )
    CONFIG.declare(
        "number_of_finite_elements",
        ConfigValue(domain=int, description="Number of finite elements to use"),
    )
    CONFIG.declare(
        "interacting_streams",
        ConfigValue(
            domain=list,
            doc="List of interacting stream pairs as 2-tuples ('stream1', 'stream2').",
        ),
    )
    # TODO: Consider a base call for heterogeneous reactions and set domain
    CONFIG.declare(
        "heterogeneous_reactions",
        ConfigValue(
            default=None,
            # domain=list,
            description="Heterogeneous reaction package to use in contactor.",
            doc="Heterogeneous reaction package to use in contactor. Heterogeneous "
            "reaction packages are expected to have a certain structure and methods; "
            "please refer to the documentation for more details.",
        ),
    )
    CONFIG.declare(
        "heterogeneous_reactions_args",
        ConfigDict(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="ConfigBlock with arguments to be passed to heterogeneous reaction block(s)",
        ),
    )

    def build(self):
        """
        Begin building model.

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build
        super().build()

        # Placeholders for things we will get from first StateBlock
        self.flow_basis = None
        self.uom = None

        self._verify_inputs()
        self._build_state_blocks()

        if self.config.heterogeneous_reactions is not None:
            self._build_heterogeneous_reaction_blocks()

        self._add_geometry()

        self._build_material_balance_constraints()
        self._build_energy_balance_constraints()
        self._build_pressure_balance_constraints()
        self._build_ports()

    def _verify_inputs(self):
        # Check that at least two streams were declared
        if len(self.config.streams) < 2:
            raise ConfigurationError(
                f"MSContactor models must define at least two streams; received "
                f"{list(self.config.streams.keys())}"
            )

        # Build indexing sets
        self.elements = RangeSet(
            1,
            self.config.number_of_finite_elements,
            doc="Set of finite elements in cascade (1 to number of elements)",
        )
        self.streams = Set(
            initialize=[k for k in self.config.streams.keys()],
            doc="Set of streams in unit",
        )

        interacting_streams = self.config.interacting_streams
        # If user did not provide interacting streams list, assume all steams interact
        if interacting_streams is None:
            interacting_streams = []
            for s1 in self.config.streams:
                for s2 in self.config.streams:
                    if not s1 == s2 and (s2, s1) not in interacting_streams:
                        interacting_streams.append((s1, s2))
        self.stream_interactions = Set(
            initialize=interacting_streams, doc="Set of interacting streams."
        )

        self.stream_component_interactions = Set(
            doc="Set of interacting components between streams."
        )
        for stream1, stream2 in self.stream_interactions:
            for j in self.config.streams[stream1].property_package.component_list:
                if (
                    j in self.config.streams[stream2].property_package.component_list
                    and (stream2, stream1, j) not in self.stream_component_interactions
                ):
                    # Common component, assume interaction
                    self.stream_component_interactions.add((stream1, stream2, j))
        if (
            len(self.stream_component_interactions) == 0
            and self.config.heterogeneous_reactions is None
        ):
            raise ConfigurationError(
                "No common components found in property packages and no heterogeneous reactions "
                "specified. The MSContactor model assumes that mass transfer occurs between "
                "components with the same name in different streams or due to heterogeneous reactions."
            )

        # Check that reaction block was provided if reactions requested
        for s, sconfig in self.config.streams.items():
            if (
                sconfig.has_rate_reactions or sconfig.has_equilibrium_reactions
            ) and sconfig.reaction_package is None:
                raise ConfigurationError(
                    f"Stream {s} was set to include reactions, but no reaction package was provided."
                )

    def _build_state_blocks(self):
        # Build state blocks
        for stream, pconfig in self.config.streams.items():
            ppack = pconfig.property_package

            arg_dict1 = dict(**pconfig.property_package_args)
            arg_dict1["defined_state"] = False

            state = ppack.build_state_block(
                self.flowsheet().time,
                self.elements,
                doc=f"States for stream {stream} in each stage.",
                **arg_dict1,
            )
            self.add_component(stream, state)

            # Add feed state if required
            if pconfig.has_feed:
                arg_dict0 = dict(**pconfig.property_package_args)
                arg_dict0["defined_state"] = True

                inlet_state = ppack.build_state_block(
                    self.flowsheet().time,
                    doc=f"Inlet states for stream {stream}.",
                    **arg_dict0,
                )
                self.add_component(stream + "_inlet_state", inlet_state)

            # Add side streams if required
            if pconfig.side_streams is not None:
                # First, verify side stream set is a sub-set of elements
                for ss in pconfig.side_streams:
                    if ss not in self.elements:
                        raise ConfigurationError(
                            f"side_streams must be a sub-set of the set of elements. "
                            f"Found {ss} in side_streams which is not in elements."
                        )
                # Add indexing Set for side streams
                ss_set = Set(initialize=pconfig.side_streams)
                self.add_component(stream + "_side_stream_set", ss_set)

                ss_state = ppack.build_state_block(
                    self.flowsheet().time,
                    ss_set,
                    doc=f"States for {stream} side streams in each stage.",
                    **arg_dict1,
                )
                self.add_component(stream + "_side_stream_state", ss_state)

            tref = self.flowsheet().time.first()
            sref = self.elements.first()
            if self.flow_basis is None:
                # Set unit level flow basis and units from first stream

                self.flow_basis = state[tref, sref].get_material_flow_basis()
                self.uom = state[tref, sref].params.get_metadata().derived_units
            else:
                # Check that flow bases are consistent
                if not state[tref, sref].get_material_flow_basis() == self.flow_basis:
                    raise ConfigurationError(
                        f"Property packages use different flow bases: ExtractionCascade "
                        f"requires all property packages to use the same basis. "
                        f"{stream} uses {state[tref, sref].get_material_flow_basis()}, "
                        f"whilst first stream uses {self.flow_basis}."
                    )

            # Build reactions blocks if provided
            if pconfig.reaction_package is not None:
                tmp_dict = dict(**pconfig.reaction_package_args)
                tmp_dict["state_block"] = state
                tmp_dict["has_equilibrium"] = pconfig.has_equilibrium_reactions

                reactions = pconfig.reaction_package.build_reaction_block(
                    self.flowsheet().time,
                    self.elements,
                    doc=f"Reaction properties for stream {stream}",
                    **tmp_dict,
                )
                self.add_component(stream + "_reactions", reactions)

    def _build_heterogeneous_reaction_blocks(self):
        rpack = self.config.heterogeneous_reactions
        rpack_args = self.config.heterogeneous_reactions_args

        try:
            self.heterogeneous_reactions = rpack.build_reaction_block(
                self.flowsheet().time,
                self.elements,
                doc="Heterogeneous reaction block for contactor.",
                **rpack_args,
            )
        except AttributeError:
            raise ConfigurationError(
                "Heterogeneous reaction package has not implemented a "
                "build_reaction_block method. Please ensure that your "
                "reaction block conforms to the required standards."
            )

        if not hasattr(self.config.heterogeneous_reactions, "reaction_idx"):
            raise PropertyNotSupportedError(
                "Heterogeneous reaction package does not contain a list of "
                "reactions (reaction_idx)."
            )

    def _add_geometry(self):
        if self.config.has_holdup:
            # Add volume for each element
            # TODO: Assuming constant volume for now
            self.volume = Var(
                self.elements,
                initialize=1,
                units=self.uom.VOLUME,
                doc="Volume of element",
            )
            self.volume_frac_stream = Var(
                self.flowsheet().time,
                self.elements,
                self.streams,
                initialize=1 / len(self.streams),
                units=units.dimensionless,
                doc="Volume fraction of each stream in element",
            )

            @self.Constraint(
                self.flowsheet().time,
                self.elements,
                doc="Sum of volume fractions constraint",
            )
            def sum_volume_frac(b, t, e):
                return 1 == sum(b.volume_frac_stream[t, e, s] for s in b.streams)

            for stream in self.config.streams.keys():
                phase_list = getattr(self, stream).phase_list
                _add_phase_fractions(self, stream, phase_list)

    def _build_material_balance_constraints(self):
        # Get units for transfer terms
        if self.flow_basis is MaterialFlowBasis.molar:
            mb_units = self.uom.FLOW_MOLE
            hu_units = self.uom.AMOUNT
        elif self.flow_basis is MaterialFlowBasis.mass:
            mb_units = self.uom.FLOW_MASS
            hu_units = self.uom.MASS
        else:
            # Flow type other, so cannot determine units
            mb_units = None
            hu_units = None

        # Material transfer terms are indexed by stream pairs and components.
        # Convention is that a positive material flow term indicates flow into
        # the first stream from the second stream for a given component.
        self.material_transfer_term = Var(
            self.flowsheet().time,
            self.elements,
            self.stream_component_interactions,
            initialize=0,
            units=mb_units,
            doc="Inter-stream mass transfer term",
        )

        if hasattr(self, "heterogeneous_reactions"):
            # Add extents of reaction and stoichiometric constraints
            # We will assume the user will define how extent will be calculated
            self.heterogeneous_reaction_extent = Var(
                self.flowsheet().time,
                self.elements,
                self.config.heterogeneous_reactions.reaction_idx,
                domain=Reals,
                initialize=0.0,
                doc="Extent of heterogeneous reactions",
                units=mb_units,
            )

        # Build balance equations
        for stream, sconfig in self.config.streams.items():
            state_block = getattr(self, stream)
            component_list = state_block.component_list
            pc_set = state_block.phase_component_set

            # Get reaction block if present
            if hasattr(self, stream + "_reactions"):
                reaction_block = getattr(self, stream + "_reactions")

            # Material holdup and accumulation
            if self.config.has_holdup:
                material_holdup = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=1.0,
                    doc="Material holdup of stream in element",
                    units=hu_units,
                )
                self.add_component(
                    stream + "_material_holdup",
                    material_holdup,
                )

                holdup_eq = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Holdup constraint for stream {stream}",
                    rule=partial(
                        _holdup_rule,
                        stream=stream,
                    ),
                )
                self.add_component(
                    stream + "_material_holdup_constraint",
                    holdup_eq,
                )

            if self.config.dynamic:
                material_accumulation = DerivativeVar(
                    material_holdup,
                    wrt=self.flowsheet().time,
                    doc="Material accumulation for in element",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_material_accumulation",
                    material_accumulation,
                )

            # Add homogeneous rate reaction terms (if required)
            if sconfig.has_rate_reactions:
                if not hasattr(sconfig.reaction_package, "rate_reaction_idx"):
                    raise PropertyNotSupportedError(
                        f"Reaction package for {stream} does not contain a list of "
                        "rate reactions (rate_reaction_idx), thus "
                        "does not support rate-based reactions."
                    )
                # Add extents of reaction and stoichiometric constraints
                # We will assume the user will define how extent will be calculated
                rate_reaction_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    sconfig.reaction_package.rate_reaction_idx,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Extent of rate-based reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_rate_reaction_extent",
                    rate_reaction_extent,
                )

                rate_reaction_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Generation due to rate-based reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_rate_reaction_generation",
                    rate_reaction_generation,
                )

                rate_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Rate-based reaction stoichiometry for stream {stream}",
                    rule=partial(
                        _rate_reaction_rule,
                        stream=stream,
                        # pylint: disable-next=possibly-used-before-assignment
                        rblock=reaction_block,
                        generation=rate_reaction_generation,
                        extent=rate_reaction_extent,
                    ),
                )
                self.add_component(
                    stream + "_rate_reaction_constraint",
                    rate_reaction_constraint,
                )

            # Add equilibrium reaction terms (if required)
            if sconfig.has_equilibrium_reactions:
                if not hasattr(sconfig.reaction_package, "equilibrium_reaction_idx"):
                    raise PropertyNotSupportedError(
                        f"Reaction package for {stream} does not contain a list of "
                        "equilibrium reactions (equilibrium_reaction_idx), thus "
                        "does not support equilibrium-based reactions."
                    )
                # Add extents of reaction and stoichiometric constraints
                equilibrium_reaction_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    sconfig.reaction_package.equilibrium_reaction_idx,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Extent of equilibrium reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_equilibrium_reaction_extent",
                    equilibrium_reaction_extent,
                )

                equilibrium_reaction_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Generation due to equilibrium reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_equilibrium_reaction_generation",
                    equilibrium_reaction_generation,
                )

                equilibrium_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Equilibrium reaction stoichiometry for stream {stream}",
                    rule=partial(
                        _equilibrium_reaction_rule,
                        stream=stream,
                        rblock=reaction_block,
                        generation=equilibrium_reaction_generation,
                        extent=equilibrium_reaction_extent,
                    ),
                )
                self.add_component(
                    stream + "_equilibrium_reaction_constraint",
                    equilibrium_reaction_constraint,
                )

            # Inherent reaction terms (if required)
            if state_block.include_inherent_reactions:
                # Add extents of reaction and stoichiometric constraints
                inherent_reaction_extent = Var(
                    self.flowsheet().time,
                    self.elements,
                    state_block.params.inherent_reaction_idx,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Extent of inherent reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_inherent_reaction_extent",
                    inherent_reaction_extent,
                )

                inherent_reaction_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc=f"Generation due to inherent reactions in stream {stream}",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_inherent_reaction_generation",
                    inherent_reaction_generation,
                )

                inherent_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc=f"Inherent reaction stoichiometry for stream {stream}",
                    rule=partial(
                        _inherent_reaction_rule,
                        stream=stream,
                        generation=inherent_reaction_generation,
                        extent=inherent_reaction_extent,
                    ),
                )
                self.add_component(
                    stream + "_inherent_reaction_constraint",
                    inherent_reaction_constraint,
                )

            # Add heterogeneous reaction terms (if required)
            if hasattr(self, "heterogeneous_reactions"):
                heterogeneous_reactions_generation = Var(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    domain=Reals,
                    initialize=0.0,
                    doc="Generation due to heterogeneous reactions",
                    units=mb_units,
                )
                self.add_component(
                    stream + "_heterogeneous_reactions_generation",
                    heterogeneous_reactions_generation,
                )

                heterogeneous_reaction_constraint = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    pc_set,
                    doc="Heterogeneous reaction stoichiometry constraint",
                    rule=partial(
                        _heterogeneous_reaction_rule,
                        pc_set=pc_set,
                        generation=heterogeneous_reactions_generation,
                    ),
                )
                self.add_component(
                    stream + "_heterogeneous_reaction_constraint",
                    heterogeneous_reaction_constraint,
                )

            # Material balance for stream
            mbal = Constraint(
                self.flowsheet().time,
                self.elements,
                component_list,
                rule=partial(
                    _material_balance_rule,
                    stream=stream,
                    mb_units=mb_units,
                ),
            )
            self.add_component(stream + "_material_balance", mbal)

    def _build_energy_balance_constraints(self):
        # Energy Balances
        for stream, pconfig in self.config.streams.items():
            if pconfig.has_energy_balance:
                state_block = getattr(self, stream)
                phase_list = state_block.phase_list

                # Material holdup and accumulation
                if self.config.has_holdup:
                    energy_holdup = Var(
                        self.flowsheet().time,
                        self.elements,
                        phase_list,
                        domain=Reals,
                        initialize=1.0,
                        doc="Energy holdup of stream in element",
                        units=self.uom.ENERGY,
                    )
                    self.add_component(
                        stream + "_energy_holdup",
                        energy_holdup,
                    )

                    energy_holdup_eq = Constraint(
                        self.flowsheet().time,
                        self.elements,
                        phase_list,
                        doc=f"Energy holdup constraint for stream {stream}",
                        rule=partial(
                            _energy_holdup_rule,
                            stream=stream,
                        ),
                    )
                    self.add_component(
                        stream + "_energy_holdup_constraint",
                        energy_holdup_eq,
                    )

                if self.config.dynamic:
                    energy_accumulation = DerivativeVar(
                        energy_holdup,
                        wrt=self.flowsheet().time,
                        doc="Energy accumulation for in element",
                        units=self.uom.POWER,
                    )
                    self.add_component(
                        stream + "_energy_accumulation",
                        energy_accumulation,
                    )

                if pconfig.has_heat_transfer:
                    heat = Var(
                        self.flowsheet().time,
                        self.elements,
                        initialize=0,
                        units=self.uom.POWER,
                        doc=f"External heat transfer term for stream {stream}",
                    )
                    self.add_component(stream + "_heat", heat)

                ebal = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    rule=partial(
                        _energy_balance_rule,
                        stream=stream,
                        uom=self.uom,
                    ),
                )
                self.add_component(stream + "_energy_balance", ebal)

    def _build_pressure_balance_constraints(self):
        # Pressure Balances
        for stream, pconfig in self.config.streams.items():
            if pconfig.has_pressure_balance:

                if pconfig.has_pressure_change:
                    deltaP = Var(
                        self.flowsheet().time,
                        self.elements,
                        initialize=0,
                        units=self.uom.PRESSURE,
                        doc=f"DeltaP term for stream {stream}",
                    )
                    self.add_component(stream + "_deltaP", deltaP)

                pbal = Constraint(
                    self.flowsheet().time,
                    self.elements,
                    rule=partial(
                        _pressure_balance_rule,
                        stream=stream,
                        uom=self.uom,
                    ),
                )
                self.add_component(stream + "_pressure_balance", pbal)

                # Add side stream pressure equality if required
                if hasattr(self, stream + "_side_stream_state"):
                    side_set = getattr(self, stream + "_side_stream_set")

                    side_pbal = Constraint(
                        self.flowsheet().time,
                        side_set,
                        rule=partial(_side_stream_pressure_rule, stream=stream),
                    )
                    self.add_component(
                        stream + "_side_stream_pressure_balance", side_pbal
                    )

    def _build_ports(self):
        # Add Ports
        for stream, pconfig in self.config.streams.items():
            sblock = getattr(self, stream)
            flow_dir = pconfig.flow_direction

            if pconfig.has_feed:
                inlet_state = getattr(self, stream + "_inlet_state")
                in_port, _ = inlet_state.build_port(
                    f"{stream} Inlet", slice_index=(slice(None))
                )
                self.add_component(stream + "_inlet", in_port)

            if flow_dir == FlowDirection.forward:
                outlet = self.elements.last()
            elif flow_dir == FlowDirection.backward:
                outlet = self.elements.first()
            else:
                raise BurntToast("If/else overrun when constructing Ports")

            out_port, _ = sblock.build_port(
                f"{stream} Outlet", slice_index=(slice(None), outlet)
            )
            self.add_component(stream + "_outlet", out_port)

    # TODO: Initialization - use the new framework
    def initialize(self, **kwargs):
        raise NotImplementedError(
            "The MSContactor unit model does not support the old initialization API. "
            "Please use the new API (InitializerObjects) instead."
        )

    def _get_performance_contents(self, time_point=0):
        # Due to the flexibility of the MSContactor and the number of possible terms
        # that could be included here, we will leave this up to the user to define.
        return {}

    def _get_stream_table_contents(self, time_point=0):
        stream_attributes = {}
        stream_attributes["Units"] = {}

        sblocks = {}
        for stream, pconfig in self.config.streams.items():
            sblock = getattr(self, stream)
            flow_dir = pconfig.flow_direction

            if pconfig.has_feed:
                inlet_state = getattr(self, stream + "_inlet_state")
                sblocks[stream + " Inlet"] = inlet_state[time_point]

            if flow_dir == FlowDirection.forward:
                outlet = self.elements.last()
            elif flow_dir == FlowDirection.backward:
                outlet = self.elements.first()
            else:
                raise BurntToast("If/else overrun when constructing stream table")

            sblocks[stream + " Outlet"] = sblock[time_point, outlet]

        for n, v in sblocks.items():
            dvars = v.define_display_vars()

            stream_attributes[n] = {}

            for k in dvars:
                for i in dvars[k].keys():
                    stream_key = k if i is None else f"{k} {i}"

                    quant = report_quantity(dvars[k][i])

                    stream_attributes[n][stream_key] = quant.m
                    stream_attributes["Units"][stream_key] = quant.u

        return DataFrame.from_dict(stream_attributes, orient="columns")


def _get_state_blocks(blk, t, s, stream):
    """
    Utility method for collecting states representing flows into and out of
    a stage for a given stream.
    """
    state_block = getattr(blk, stream)

    if blk.config.streams[stream].flow_direction == FlowDirection.forward:
        if s == blk.elements.first():
            if not blk.config.streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(blk, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, blk.elements.prev(s)]
    elif blk.config.streams[stream].flow_direction == FlowDirection.backward:
        if s == blk.elements.last():
            if not blk.config.streams[stream].has_feed:
                in_state = None
            else:
                in_state = getattr(blk, stream + "_inlet_state")[t]
        else:
            in_state = state_block[t, blk.elements.next(s)]
    else:
        raise BurntToast("If/else overrun when constructing balances")

    out_state = state_block[t, s]

    # Look for side state
    side_state = None
    if hasattr(blk, stream + "_side_stream_state"):
        try:
            side_state = getattr(blk, stream + "_side_stream_state")[t, s]
        except KeyError:
            pass

    return in_state, out_state, side_state


def _rate_reaction_rule(blk, t, s, p, j, stream, rblock, generation, extent):
    sconfig = blk.config.streams[stream]
    sblock = getattr(blk, stream)
    pc_set = sblock.phase_component_set

    if (p, j) in pc_set:
        return generation[t, s, p, j] == (
            sum(
                rblock[t, s].params.rate_reaction_stoichiometry[r, p, j]
                * extent[t, s, r]
                for r in sconfig.reaction_package.rate_reaction_idx
            )
        )
    return Constraint.Skip


def _equilibrium_reaction_rule(blk, t, s, p, j, stream, rblock, generation, extent):
    sconfig = blk.config.streams[stream]
    sblock = getattr(blk, stream)
    pc_set = sblock.phase_component_set

    if (p, j) in pc_set:
        return generation[t, s, p, j] == (
            sum(
                rblock[t, s].params.equilibrium_reaction_stoichiometry[r, p, j]
                * extent[t, s, r]
                for r in sconfig.reaction_package.equilibrium_reaction_idx
            )
        )
    return Constraint.Skip


def _inherent_reaction_rule(blk, t, s, p, j, stream, generation, extent):
    sconfig = blk.config.streams[stream]
    sblock = getattr(blk, stream)
    pc_set = sblock.phase_component_set

    if (p, j) in pc_set:
        return generation[t, s, p, j] == (
            sum(
                sblock[t, s].params.inherent_reaction_stoichiometry[r, p, j]
                * extent[t, s, r]
                for r in sconfig.property_package.inherent_reaction_idx
            )
        )
    return Constraint.Skip


def _heterogeneous_reaction_rule(blk, t, s, p, j, pc_set, generation):
    if (p, j) in pc_set:
        return generation[t, s, p, j] == (
            sum(
                blk.heterogeneous_reactions[t, s].params.reaction_stoichiometry[r, p, j]
                * blk.heterogeneous_reaction_extent[t, s, r]
                for r in blk.config.heterogeneous_reactions.reaction_idx
                if (r, p, j)
                in blk.heterogeneous_reactions[t, s].params.reaction_stoichiometry
            )
        )
    return Constraint.Skip


def _material_balance_rule(blk, t, s, j, stream, mb_units):
    sconfig = blk.config.streams[stream]
    state_block = getattr(blk, stream)
    phase_list = state_block.phase_list
    pc_set = state_block.phase_component_set

    in_state, out_state, side_state = _get_state_blocks(blk, t, s, stream)

    if in_state is not None:
        rhs = sum(
            in_state.get_material_flow_terms(p, j)
            for p in phase_list
            if (p, j) in pc_set
        ) - sum(
            out_state.get_material_flow_terms(p, j)
            for p in phase_list
            if (p, j) in pc_set
        )
    else:
        rhs = -sum(
            out_state.get_material_flow_terms(p, j)
            for p in phase_list
            if (p, j) in pc_set
        )

    # Add side stream energy flow if required
    if side_state is not None:
        rhs += sum(
            side_state.get_material_flow_terms(p, j)
            for p in phase_list
            if (p, j) in pc_set
        )

    # As overall units come from the first StateBlock constructed, we
    # cannot guarantee that any units are consistent, so convert all flow terms
    if mb_units is not None:
        rhs = units.convert(rhs, mb_units)

    # Add mass transfer terms
    for k in blk.stream_component_interactions:
        if k[0] == stream and k[2] == j:
            # Positive mass transfer
            rhs += blk.material_transfer_term[t, s, k]
        elif k[1] == stream and k[2] == j:
            # Negative mass transfer
            rhs += -blk.material_transfer_term[t, s, k]

    # Add rate reactions (if required)
    if sconfig.has_rate_reactions:
        rhs += sum(
            getattr(
                blk,
                stream + "_rate_reaction_generation",
            )[t, s, p, j]
            for p in phase_list
        )

    # Add equilibrium reactions (if required)
    if sconfig.has_equilibrium_reactions:
        rhs += sum(
            getattr(
                blk,
                stream + "_equilibrium_reaction_generation",
            )[t, s, p, j]
            for p in phase_list
        )

    # Add inherent reactions (if required)
    if state_block.include_inherent_reactions:
        rhs += sum(
            getattr(
                blk,
                stream + "_inherent_reaction_generation",
            )[t, s, p, j]
            for p in phase_list
        )

    # Add heterogeneous reactions (if required)
    if blk.config.heterogeneous_reactions is not None:
        rhs += sum(
            getattr(
                blk,
                stream + "_heterogeneous_reactions_generation",
            )[t, s, p, j]
            for p in phase_list
        )

    if not blk.config.dynamic:
        lhs = 0 * mb_units
    else:
        acc = getattr(blk, stream + "_material_accumulation")
        lhs = sum(acc[t, s, p, j] for p in phase_list)
        if mb_units is not None:
            lhs = units.convert(lhs, mb_units)

    return lhs == rhs


def _get_energy_transfer_term(blk, uom):
    try:
        return blk.energy_transfer_term
    except AttributeError:
        # Assume that if energy balances are enabled that energy transfer
        # occurs between all interacting phases.
        # For now, we will not distinguish different types of energy transfer.
        # Convention is that a positive energy flow term indicates flow into
        # the first stream from the second stream.
        blk.energy_transfer_term = Var(
            blk.flowsheet().time,
            blk.elements,
            blk.stream_interactions,
            initialize=0,
            units=uom.POWER,
            doc="Inter-stream energy transfer term",
        )

    return blk.energy_transfer_term


def _energy_balance_rule(blk, t, s, stream, uom):
    pconfig = blk.config.streams[stream]
    state_block = getattr(blk, stream)
    phase_list = state_block.phase_list

    in_state, out_state, side_state = _get_state_blocks(blk, t, s, stream)

    if in_state is not None:
        rhs = sum(in_state.get_enthalpy_flow_terms(p) for p in phase_list) - sum(
            out_state.get_enthalpy_flow_terms(p) for p in phase_list
        )
    else:
        rhs = -sum(out_state.get_enthalpy_flow_terms(p) for p in phase_list)

    # Add side stream energy flow if required
    if side_state is not None:
        rhs += sum(side_state.get_enthalpy_flow_terms(p) for p in phase_list)

    # As overall units come from the first StateBlock constructed, we
    # cannot guarantee that any units are consistent, so convert all flow terms
    rhs = units.convert(rhs, uom.POWER)

    # Add interstream transfer terms
    for k in blk.stream_interactions:
        ett = _get_energy_transfer_term(blk, uom)
        if k[0] == stream:
            # Positive energy transfer
            rhs += ett[t, s, k]
        elif k[1] == stream:
            # Negative energy transfer
            rhs += -ett[t, s, k]

    # Add external heat term if required
    if pconfig.has_heat_transfer:
        rhs += getattr(blk, stream + "_heat")[t, s]

    # Add heat of reaction terms if required
    if pconfig.has_heat_of_reaction:
        if not (
            hasattr(blk, stream + "_rate_reaction_extent")
            or hasattr(blk, stream + "_equilibrium_reaction_extent")
        ):
            raise ConfigurationError(
                f"Stream {stream} was set to include heats of reaction, "
                "but no extent of reactions terms could be found. "
                "Please ensure that you defined a reaction package for this "
                "stream and that the material balances were set to include "
                "reactions."
            )
        reactions = getattr(blk, stream + "_reactions")

        if hasattr(blk, stream + "_rate_reaction_extent"):
            rate_extent = getattr(blk, stream + "_rate_reaction_extent")
            rhs += -sum(
                rate_extent[t, s, r] * reactions[t, s].dh_rxn[r]
                for r in pconfig.reaction_package.rate_reaction_idx
            )

        if hasattr(blk, stream + "_equilibrium_reaction_extent"):
            equil_extent = getattr(blk, stream + "_equilibrium_reaction_extent")
            rhs += -sum(
                equil_extent[t, s, e] * reactions[t, s].dh_rxn[e]
                for e in pconfig.reaction_package.equilibrium_reaction_idx
            )

    if not blk.config.dynamic:
        lhs = 0 * uom.POWER
    else:
        acc = getattr(blk, stream + "_energy_accumulation")
        lhs = units.convert(sum(acc[t, s, p] for p in phase_list), uom.POWER)

    return lhs == rhs


def _pressure_balance_rule(blk, t, s, stream, uom):
    pconfig = blk.config.streams[stream]
    in_state, out_state, _ = _get_state_blocks(blk, t, s, stream)

    if in_state is None:
        # If there is no feed, then there is no need for a pressure balance
        return Constraint.Skip

    rhs = in_state.pressure - out_state.pressure

    # As overall units come from the first StateBlock constructed, we
    # cannot guarantee that any units are consistent, so convert all flow terms
    rhs = units.convert(rhs, uom.PRESSURE)

    # Add deltaP term if required
    if pconfig.has_pressure_change:
        rhs += getattr(blk, stream + "_deltaP")[t, s]

    return 0 * uom.PRESSURE == rhs


def _side_stream_pressure_rule(b, t, s, stream):
    stage_state = getattr(b, stream)[t, s]
    side_state = getattr(b, stream + "_side_stream_state")[t, s]

    return stage_state.pressure == side_state.pressure


def _holdup_rule(b, t, e, p, j, stream):
    holdup = getattr(b, stream + "_material_holdup")
    stage_state = getattr(b, stream)[t, e]
    phase_frac = getattr(b, stream + "_phase_fraction")

    return holdup[t, e, p, j] == (
        b.volume[e]
        * b.volume_frac_stream[t, e, stream]
        * phase_frac[t, e, p]
        * stage_state.get_material_density_terms(p, j)
    )


def _energy_holdup_rule(b, t, e, p, stream):
    holdup = getattr(b, stream + "_energy_holdup")
    stage_state = getattr(b, stream)[t, e]
    phase_frac = getattr(b, stream + "_phase_fraction")

    return holdup[t, e, p] == (
        b.volume[e]
        * b.volume_frac_stream[t, e, stream]
        * phase_frac[t, e, p]
        * stage_state.get_energy_density_terms(p)
    )


def _add_phase_fractions(b, stream, phase_list):
    if len(phase_list) > 1:
        phase_fraction = Var(
            b.flowsheet().time,
            b.elements,
            phase_list,
            initialize=1 / len(phase_list),
            doc=f"Volume fraction of holdup by phase in stream {stream}",
        )
        b.add_component(stream + "_phase_fraction", phase_fraction)

        sum_of_phase_fractions = Constraint(
            b.flowsheet().time,
            b.elements,
            rule=partial(
                _sum_phase_frac_rule, phase_frac=phase_fraction, phase_list=phase_list
            ),
            doc=f"Sum of phase fractions constraint for stream {stream}",
        )
        b.add_component(stream + "_sum_phase_fractions", sum_of_phase_fractions)

    else:

        def phase_frac_rule(b, t, x, p):
            return 1

        phase_fraction = Expression(
            b.flowsheet().time,
            b.elements,
            phase_list,
            rule=phase_frac_rule,
            doc=f"Volume fraction of holdup by phase in stream {stream}",
        )

        b.add_component(stream + "_phase_fraction", phase_fraction)


def _sum_phase_frac_rule(b, t, x, phase_frac, phase_list):
    return 1 == sum(phase_frac[t, x, p] for p in phase_list)
