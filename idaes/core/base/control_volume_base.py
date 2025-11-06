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
Base class for control volumes
"""

# Import Python libraries
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import value
from pyomo.common.collections import ComponentMap
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES cores
from idaes.core import (
    ProcessBlockData,
    MaterialFlowBasis,
    useDefault,
    declare_process_block_class,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
    DefaultBool,
)
from idaes.core.util.exceptions import (
    BurntToast,
    ConfigurationError,
    PropertyNotSupportedError,
)
from idaes.core.scaling import (
    ConstraintScalingScheme,
    CustomScalerBase,
    get_scaling_factor,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__author__ = "Andrew Lee, Douglas Allan"


# Enumerate options for material balances
class MaterialBalanceType(Enum):
    """
    Enum for material balance types.
    """

    useDefault = -1
    none = 0
    componentPhase = 1
    componentTotal = 2
    elementTotal = 3
    total = 4


# Enumerate options for energy balances
class EnergyBalanceType(Enum):
    """
    Enum for energy balance types.
    """

    useDefault = -1
    none = 0
    enthalpyPhase = 1
    enthalpyTotal = 2
    energyPhase = 3
    energyTotal = 4
    isothermal = 5


# Enumerate options for momentum balances
class MomentumBalanceType(Enum):
    """
    Enum for momentum/pressure balance types.
    """

    none = 0
    pressureTotal = 1
    pressurePhase = 2
    momentumTotal = 3
    momentumPhase = 4


# Enumerate options for flow direction
class FlowDirection(Enum):
    """
    Enum indicating direction of flow.
    """

    notSet = 0
    forward = 1
    backward = 2


class ControlVolumeScalerBase(CustomScalerBase):
    """
    Scaler object for elements common to the ControlVolume0D and ControlVolume1D
    """

    # Attribute name to use as a weight when scaling material and energy
    # terms. Presently (9/25/25), this attribute exists to take into
    # account the fact that all the material and energy terms in the
    # ControlVolume1D are given on the basis of material or energy per
    # unit length, so we want to weight them accordingly.
    _weight_attr_name = None

    def _get_reference_state_block(self, model):
        """
        This method gives the parent class ControlVolumeScalerBase
        methods a state block with the same index as the material
        and energy balances to get scaling information from
        """
        raise NotImplementedError(
            "This method is intended to be overridden by other scaler "
            "objects inheriting from it."
        )

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers=None
    ):
        """
        Routine to apply scaling factors to variables in model.

        Derived classes must overload this method.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scalers to use for sub-models

        Returns:
            None
        """
        props = self._get_reference_state_block(model)
        phase_list = props.phase_list
        phase_component_set = props.phase_component_set

        phase_equilibrium_idx = getattr(
            props.params,
            "phase_equilibrium_idx",
            None,  # Default value if attr does not exist
        )
        phase_equilibrium_list = getattr(
            props.params,
            "phase_equilibrium_list",
            None,  # Default value if attr does not exist
        )

        if self._weight_attr_name is None:
            weight = 1
        else:
            # For ControlVolume1D, the weight is L and the
            # terms have units of material per length or
            # energy per length. The scaling factor of L has
            # units of 1 / L, so we want to divide by the scaling
            # factor to render the material or energy terms
            # dimensionless.
            weight_attr = getattr(model, self._weight_attr_name)
            weight = 1 / self.get_scaling_factor(weight_attr, default=1, warning=True)

        idx0 = props.index_set().first()
        params = props[idx0].params
        if hasattr(model, "reactions"):
            self.call_submodel_scaler_method(
                model.reactions,
                submodel_scalers=submodel_scalers,
                method="variable_scaling_routine",
                overwrite=overwrite,
            )

        # Set scaling factors for common material balance variables
        if hasattr(model, "material_holdup"):
            for idx, v in model.material_holdup.items():
                self.scale_variable_by_definition_constraint(
                    v, model.material_holdup_calculation[idx], overwrite=overwrite
                )

        # Material accumulation should be scaled by a global method for scaling
        # time DerivativeVars
        if hasattr(model, "material_accumulation"):
            pass

        # Using a heuristic to estimate scaling factors for reaction generation
        # terms and extents of reaction. This heuristic may not be suitable when
        # reactions generate and consume large quantities of an intermediate
        # species

        # Rate reactions
        if hasattr(model, "rate_reaction_generation"):
            stoich = model.reactions[idx0].params.rate_reaction_stoichiometry
            rate_rxn_gen = getattr(model, "rate_reaction_generation")
            rate_rxn_idx = model.config.reaction_package.rate_reaction_idx
            # Material generation scaling is based on the magnitude of
            # the material flow terms
            for idx in rate_rxn_gen:
                prop_idx = idx[:-2]
                p = idx[-2]
                j = idx[-1]
                nom = self.get_expression_nominal_value(
                    props[prop_idx].get_material_flow_terms(p, j)
                )
                self.set_component_scaling_factor(
                    rate_rxn_gen[idx], weight / nom, overwrite=overwrite
                )

            # Extent of reaction scaling is based on the species in the
            # reaction which has the largest scaling factor (and thus is
            # the species whose concentration is the most sensitive).
            # This scaling may not work for systems with highly reactive
            # intermediates, in which multiple extents of reaction
            # cancel each other out (but if that's the case, we either
            # should be able to eliminate the highly reactive species from
            # the reaction system by combining reactions or we need to keep
            # track of the concentration of this highly reactive intermediate.)
            for prop_idx in props:
                for rxn in rate_rxn_idx:
                    nom_rxn = float("inf")
                    for p, j in phase_component_set:
                        sf_pc = get_scaling_factor(rate_rxn_gen[prop_idx, p, j])
                        coeff = stoich[rxn, p, j]
                        if coeff != 0:
                            nom_rxn = min(abs(coeff) / sf_pc, nom_rxn)
                    if nom_rxn == float("inf"):
                        raise ConfigurationError(
                            f"Reaction {rxn} has no nonzero stoichiometric coefficient."
                        )
                    # Note this scaling works only if we don't
                    # have multiple reactions cancelling each other out.
                    # No need to weight here because the rate_rxn_gen
                    # scaling factor is already weighted
                    self.set_component_scaling_factor(
                        model.rate_reaction_extent[prop_idx, rxn],
                        1 / nom_rxn,
                        overwrite=overwrite,
                    )

        # Equilibrium reaction
        if hasattr(model, "equilibrium_reaction_generation"):
            stoich = model.reactions[idx0].params.equilibrium_reaction_stoichiometry
            equil_rxn_gen = getattr(model, "equilibrium_reaction_generation")
            equil_rxn_idx = model.config.reaction_package.equilibrium_reaction_idx
            # Material generation scaling is based on the magnitude of
            # the material flow terms
            for idx in equil_rxn_gen:
                prop_idx = idx[:-2]
                p = idx[-2]
                j = idx[-1]
                nom = self.get_expression_nominal_value(
                    props[prop_idx].get_material_flow_terms(p, j)
                )
                self.set_component_scaling_factor(
                    equil_rxn_gen[idx], weight / nom, overwrite=overwrite
                )

            # Extent of reaction scaling is based on the species in the
            # reaction which has the largest scaling factor (and thus is
            # the species whose concentration is the most sensitive).
            # This scaling may not work for systems with highly reactive
            # intermediates, in which multiple extents of reaction
            # cancel each other out (but if that's the case, we either
            # should be able to eliminate the highly reactive species from
            # the reaction system by combining reactions or we need to keep
            # track of the concentration of this highly reactive intermediate.)
            for prop_idx in props:
                for rxn in equil_rxn_idx:
                    nom_rxn = float("inf")
                    for p, j in phase_component_set:
                        sf_pc = get_scaling_factor(equil_rxn_gen[prop_idx, p, j])
                        coeff = stoich[rxn, p, j]
                        if coeff != 0:
                            nom_rxn = min(abs(coeff) / sf_pc, nom_rxn)
                    if nom_rxn == float("inf"):
                        raise ConfigurationError(
                            f"Reaction {rxn} has no nonzero stoichiometric coefficient."
                        )
                    # Note this scaling works only if we don't
                    # have multiple reactions cancelling each other out
                    # No need to weight here because the equil_rxn_gen
                    # scaling factor is already weighted
                    self.set_component_scaling_factor(
                        model.equilibrium_reaction_extent[prop_idx, rxn],
                        1 / nom_rxn,
                        overwrite=overwrite,
                    )

        # Inherent reaction
        if hasattr(model, "inherent_reaction_generation"):
            stoich = params.inherent_reaction_stoichiometry
            inh_rxn_gen = getattr(model, "inherent_reaction_generation")
            inh_rxn_idx = params.inherent_reaction_idx
            # Material generation scaling is based on the magnitude of
            # the material flow terms
            for idx in inh_rxn_gen:
                prop_idx = idx[:-2]
                p = idx[-2]
                j = idx[-1]
                nom = self.get_expression_nominal_value(
                    props[prop_idx].get_material_flow_terms(p, j)
                )
                self.set_component_scaling_factor(
                    inh_rxn_gen[idx], weight / nom, overwrite=overwrite
                )

            # Extent of reaction scaling is based on the species in the
            # reaction which has the largest scaling factor (and thus is
            # the species whose concentration is the most sensitive).
            # This scaling may not work for systems with highly reactive
            # intermediates, in which multiple extents of reaction
            # cancel each other out (but if that's the case, we either
            # should be able to eliminate the highly reactive species from
            # the reaction system by combining reactions or we need to keep
            # track of the concentration of this highly reactive intermediate.)
            for prop_idx in props:
                for rxn in inh_rxn_idx:
                    nom_rxn = float("inf")
                    for p, j in phase_component_set:
                        sf_pc = get_scaling_factor(inh_rxn_gen[prop_idx, p, j])
                        coeff = stoich[rxn, p, j]
                        if coeff != 0:
                            nom_rxn = min(abs(coeff) / sf_pc, nom_rxn)
                    if nom_rxn == float("inf"):
                        raise ConfigurationError(
                            f"Reaction {rxn} has no nonzero stoichiometric coefficient."
                        )
                    # Note this scaling works only if we don't
                    # have multiple reactions cancelling each other out
                    # No need to weight here because the inh_rxn_gen
                    # scaling factor is already weighted
                    self.set_component_scaling_factor(
                        model.inherent_reaction_extent[prop_idx, rxn],
                        1 / nom_rxn,
                        overwrite=overwrite,
                    )

        if hasattr(model, "mass_transfer_term"):
            for prop_idx in props:
                for p, j in phase_component_set:
                    nom = self.get_expression_nominal_value(
                        props[prop_idx].get_material_flow_terms(p, j)
                    )
                    self.set_component_scaling_factor(
                        model.mass_transfer_term[prop_idx, p, j],
                        weight / nom,
                        overwrite=overwrite,
                    )

        if hasattr(model, "phase_equilibrium_generation"):
            for prop_idx in props:
                for pe_idx in phase_equilibrium_idx:
                    j, pp = phase_equilibrium_list[pe_idx]
                    nom1 = self.get_expression_nominal_value(
                        props[prop_idx].get_material_flow_terms(pp[0], j)
                    )
                    nom2 = self.get_expression_nominal_value(
                        props[prop_idx].get_material_flow_terms(pp[1], j)
                    )
                    nom = min(nom1, nom2)
                    self.set_component_scaling_factor(
                        model.phase_equilibrium_generation[prop_idx, pe_idx],
                        1 / nom,
                        overwrite=False,
                    )

        # Set scaling factors for element balance variables
        # TODO
        # if hasattr(model, "elemental_flow_out"):
        #     for prop_idx in props:
        #         for e in params.element_list:
        #             flow_basis = model.properties_out[t].get_material_flow_basis()
        #             for p, j in phase_component_set:

        #             sf = iscale.min_scaling_factor(
        #                 [
        #                     model.properties_out[t].get_material_density_terms(p, j)
        #                     for (p, j) in phase_component_set
        #                 ],
        #                 default=1,
        #                 warning=True,
        #             )
        #             if flow_basis == MaterialFlowBasis.molar:
        #                 sf *= 1
        #             elif flow_basis == MaterialFlowBasis.mass:
        #                 # MW scaling factor is the inverse of its value
        #                 sf *= value(model.properties_out[t].mw_comp[j])

        #             iscale.set_scaling_factor(v, sf)
        #             iscale.set_scaling_factor(model.elemental_flow_in[t, p, e], sf)

        # if hasattr(model, "element_holdup"):
        #     for (t, e), v in model.element_holdup.items():
        #         flow_basis = model.properties_out[t].get_material_flow_basis()
        #         sf_list = []
        #         for p, j in phase_component_set:
        #             if flow_basis == MaterialFlowBasis.molar:
        #                 sf = 1
        #             elif flow_basis == MaterialFlowBasis.mass:
        #                 # MW scaling factor is the inverse of its value
        #                 sf = value(model.properties_out[t].mw_comp[j])
        #             sf *= get_scaling_factor(model.phase_fraction[t, p])
        #             sf *= get_scaling_factor(
        #                 model.properties_out[t].get_material_density_terms(p, j),
        #                 default=1,
        #                 warning=True,
        #             )
        #             sf *= value(model.properties_out[t].params.element_comp[j][e]) ** -1
        #             sf_list.append(sf)
        #         sf_h = min(sf_list) * get_scaling_factor(model.volume[t])
        #         iscale.set_scaling_factor(v, sf_h)

        # if hasattr(model, "element_accumulation"):
        #     for (t, e), v in model.element_accumulation.items():
        #         if get_scaling_factor(v) is None:
        #             sf = iscale.min_scaling_factor(
        #                 model.elemental_flow_out[t, ...], default=1, warning=True
        #             )
        #             iscale.set_scaling_factor(v, sf)

        # if hasattr(model, "elemental_mass_transfer_term"):
        #     for (t, e), v in model.elemental_mass_transfer_term.items():
        #         # minimum scaling factor for elemental_flow terms
        #         sf_list = []
        #         flow_basis = model.properties_out[t].get_material_flow_basis()
        #         if get_scaling_factor(v) is None:
        #             sf = iscale.min_scaling_factor(
        #                 model.elemental_flow_out[t, ...], default=1, warning=True
        #             )
        #             iscale.set_scaling_factor(v, sf)

        # Set scaling factors for energy balance variables
        if hasattr(model, "energy_holdup"):
            for idx in model.energy_holdup:
                self.scale_variable_by_definition_constraint(
                    model.energy_holdup[idx], model.energy_holdup_calculation[idx]
                )
        # Energy accumulation should be scaled by a global method for scaling
        # time derivative variables
        if hasattr(model, "energy_accumulation"):
            pass

        # Energy transfer terms
        if (
            hasattr(model, "heat")
            or hasattr(model, "work")
            or hasattr(model, "enthalpy_transfer")
        ):
            for prop_idx in props:
                nom_list = []
                for p in phase_list:
                    nom_list.append(
                        self.get_expression_nominal_value(
                            props[prop_idx].get_enthalpy_flow_terms(p)
                        )
                    )
                # TODO we need to do some validation so that nom isn't zero or near-zero
                nom = max(nom_list)
                if hasattr(model, "heat"):
                    self.set_component_scaling_factor(
                        model.heat[prop_idx], weight / nom, overwrite=overwrite
                    )
                if hasattr(model, "work"):
                    self.set_component_scaling_factor(
                        model.work[prop_idx], weight / nom, overwrite=overwrite
                    )
                if hasattr(model, "enthalpy_transfer"):
                    self.set_component_scaling_factor(
                        model.enthalpy_transfer[prop_idx],
                        weight / nom,
                        overwrite=overwrite,
                    )

        # Set scaling for momentum balance variables
        if hasattr(model, "deltaP"):
            for prop_idx in props:
                sf_P = get_scaling_factor(
                    props[prop_idx].pressure, default=1e-5, warning=True
                )
                self.set_component_scaling_factor(
                    model.deltaP[prop_idx], weight * sf_P, overwrite=overwrite
                )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        """
        Routine to apply scaling factors to constraints in model.

        Derived classes must overload this method.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: ComponentMap of Scalers to use for sub-models

        Returns:
            None
        """
        props = self._get_reference_state_block(model)
        phase_list = props.phase_list
        pc_set = props.phase_component_set

        if hasattr(model, "reactions"):
            self.call_submodel_scaler_method(
                model.reactions,
                submodel_scalers=submodel_scalers,
                method="constraint_scaling_routine",
                overwrite=overwrite,
            )

        if hasattr(model, "material_holdup_calculation"):
            for idx in model.material_holdup_calculation:
                self.scale_constraint_by_component(
                    model.material_holdup_calculation[idx],
                    model.material_holdup[idx],
                    overwrite=overwrite,
                )

        if hasattr(model, "rate_reaction_stoichiometry_constraint"):
            for idx in model.rate_reaction_stoichiometry_constraint:
                self.scale_constraint_by_component(
                    model.rate_reaction_stoichiometry_constraint[idx],
                    model.rate_reaction_generation[idx],
                    overwrite=overwrite,
                )

        if hasattr(model, "equilibrium_reaction_stoichiometry_constraint"):
            for idx in model.equilibrium_reaction_stoichiometry_constraint:
                self.scale_constraint_by_component(
                    model.equilibrium_reaction_stoichiometry_constraint[idx],
                    model.equilibrium_reaction_generation[idx],
                    overwrite=overwrite,
                )

        inh_rxn_con = None
        if hasattr(model, "inherent_reaction_stoichiometry_constraint"):
            # ControlVolume0D and ControlVolume1D
            inh_rxn_con = model.inherent_reaction_stoichiometry_constraint
        elif hasattr(model, "inherent_reaction_constraint"):
            # Mixer
            inh_rxn_con = model.inherent_reaction_constraint

        if inh_rxn_con is not None:
            for idx in inh_rxn_con:
                self.scale_constraint_by_component(
                    inh_rxn_con[idx],
                    model.inherent_reaction_generation[idx],
                    overwrite=overwrite,
                )

        mb_eqn = None
        if hasattr(model, "material_balances"):
            # ControlVolume0D and ControlVolume1D
            mb_eqn = model.material_balances
        elif hasattr(model, "material_mixing_equations"):
            # Mixer
            mb_eqn = model.material_mixing_equations

        if mb_eqn is not None:
            mb_type = model._constructed_material_balance_type  # pylint: disable=W0212
            if mb_type == MaterialBalanceType.componentTotal:
                for idx in mb_eqn:
                    c = idx[-1]
                    nom_list = []
                    for p in phase_list:
                        if (p, c) in props.phase_component_set:
                            nom_list.append(
                                self.get_expression_nominal_value(
                                    props[idx[:-1]].get_material_flow_terms(p, c)
                                )
                            )
                    nom = max(nom_list)
                    self.set_component_scaling_factor(
                        mb_eqn[idx], 1 / nom, overwrite=overwrite
                    )
            elif mb_type == MaterialBalanceType.componentPhase:
                for idx in mb_eqn:
                    p = idx[-2]
                    c = idx[-1]
                    nom = self.get_expression_nominal_value(
                        props[idx[:-2]].get_material_flow_terms(p, c)
                    )
                    self.set_component_scaling_factor(
                        mb_eqn[idx], 1 / nom, overwrite=overwrite
                    )
            elif mb_type == MaterialBalanceType.total:
                for idx in mb_eqn:
                    nom_list = []
                    for p, c in pc_set:
                        nom_list.append(
                            self.get_expression_nominal_value(
                                props[idx[:-1]].get_material_flow_terms(p, c)
                            )
                        )
                    nom = max(nom_list)
                    self.set_component_scaling_factor(
                        mb_eqn[idx], 1 / nom, overwrite=overwrite
                    )
            else:
                # There are some other material balance types but they create
                # constraints with different names.
                _log.warning(
                    f"Unknown material balance type {mb_type}. It cannot be "
                    "automatically scaled."
                )

        # TODO element balances
        # if hasattr(self, "element_balances"):
        #     for (t, e), c in self.element_balances.items():
        #         sf = iscale.min_scaling_factor(
        #             [self.elemental_flow_out[t, p, e] for p in phase_list]
        #         )
        #         iscale.constraint_scaling_transform(c, sf, overwrite=False)

        # if hasattr(self, "elemental_holdup_calculation"):
        #     for (t, e), c in self.elemental_holdup_calculation.items():
        #         sf = iscale.get_scaling_factor(self.element_holdup[t, e])
        #         iscale.constraint_scaling_transform(c, sf, overwrite=False)

        eb_eqn = None
        if hasattr(model, "enthalpy_balances"):
            eb_eqn = model.enthalpy_balances
        elif hasattr(model, "enthalpy_mixing_equations"):
            # Mixer
            eb_eqn = model.enthalpy_mixing_equations

        if eb_eqn is not None:
            # Phase enthalpy balances are not implemented
            # as of 9/26/25
            for idx in eb_eqn:
                nom_list = []
                for p in phase_list:
                    nom_list.append(
                        self.get_expression_nominal_value(
                            props[idx].get_enthalpy_flow_terms(p)
                        )
                    )
                nom = max(nom_list)
                self.set_component_scaling_factor(
                    eb_eqn[idx], 1 / nom, overwrite=overwrite
                )

        if hasattr(model, "energy_holdup_calculation"):
            for idx in model.energy_holdup_calculation:
                self.scale_constraint_by_component(
                    model.energy_holdup_calculation[idx],
                    model.energy_holdup[idx],
                    overwrite=overwrite,
                )

        pb_eqn = None
        if hasattr(model, "pressure_balance"):
            # ControlVolume0D and ControlVolume1D
            pb_eqn = model.pressure_balance

        if pb_eqn is not None:
            for idx, con in pb_eqn.items():
                self.scale_constraint_by_component(
                    con,
                    props[idx].pressure,
                    overwrite=overwrite,
                )

        if hasattr(model, "sum_of_phase_fractions"):
            for con in model.sum_of_phase_fractions.values():
                self.scale_constraint_by_nominal_value(
                    con,
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )

        # Scaling for discretization equations
        # These equations should be scaled by a global method to scale time discretization equations
        if hasattr(model, "material_accumulation_disc_eq"):
            pass

        if hasattr(model, "energy_accumulation_disc_eq"):
            pass

        if hasattr(model, "element_accumulation_disc_eq"):
            pass


# Set up example ConfigBlock that will work with ControlVolume autobuild method
CONFIG_Template = ProcessBlockData.CONFIG()
CONFIG_Template.declare(
    "dynamic",
    ConfigValue(
        default=useDefault,
        domain=DefaultBool,
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic,
**default** - useDefault.
**Valid values:** {
**useDefault** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}""",
    ),
)
CONFIG_Template.declare(
    "has_holdup",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
    ),
)
CONFIG_Template.declare(
    "material_balance_type",
    ConfigValue(
        default=MaterialBalanceType.componentPhase,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
    ),
)
CONFIG_Template.declare(
    "energy_balance_type",
    ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
    ),
)
CONFIG_Template.declare(
    "momentum_balance_type",
    ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
    ),
)
CONFIG_Template.declare(
    "has_rate_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Rate reaction construction flag",
        doc="""Indicates whether terms for rate controlled reactions should be
constructed,
**default** - False.
**Valid values:** {
**True** - include kinetic reaction terms,
**False** - exclude kinetic reaction terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_equilibrium_reactions",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Equilibrium reaction construction flag",
        doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - False.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_phase_equilibrium",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Phase equilibrium construction flag",
        doc="""Indicates whether terms for phase equilibrium should be
constructed,
**default** = False.
**Valid values:** {
**True** - include phase equilibrium terms
**False** - exclude phase equilibrium terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_mass_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Mass transfer term construction flag",
        doc="""Indicates whether terms for mass transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include mass transfer terms,
**False** - exclude mass transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_heat_of_reaction",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Heat of reaction term construction flag",
        doc="""Indicates whether terms for heat of reaction should be constructed,
**default** - False.
**Valid values** {
**True** - include heat of reaction terms,
**False** - exclude heat of reaction terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_heat_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_work_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Work transfer term construction flag",
        doc="""Indicates whether terms for work transfer should be constructed,
**default** - False.
**Valid values** {
**True** - include work transfer terms,
**False** - exclude work transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_enthalpy_transfer",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Enthalpy transfer term construction flag",
        doc="""Indicates whether terms for enthalpy transfer due to mass transfer
should be constructed, **default** - False.
**Valid values** {
**True** - include enthalpy transfer terms,
**False** - exclude enthalpy transfer terms.}""",
    ),
)
CONFIG_Template.declare(
    "has_pressure_change",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
    ),
)
CONFIG_Template.declare(
    "property_package",
    ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
    ),
)
CONFIG_Template.declare(
    "property_package_args",
    ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
    ),
)
CONFIG_Template.declare(
    "reaction_package",
    ConfigValue(
        default=None,
        domain=is_reaction_parameter_block,
        description="Reaction package to use for control volume",
        doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
    ),
)
CONFIG_Template.declare(
    "reaction_package_args",
    ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing reaction packages",
        doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
    ),
)


@declare_process_block_class(
    "ControlVolume",
    doc="This class is not usually used directly. "
    "Use ControlVolume0DBlock or ControlVolume1DBlock"
    " instead.",
)
class ControlVolumeBlockData(ProcessBlockData):
    """
    The ControlVolumeBlockData Class forms the base class for all IDAES
    ControlVolume models. The purpose of this class is to automate the tasks
    common to all control volume blocks and ensure that the necessary
    attributes of a control volume block are present.

    The most signfiicant role of the ControlVolumeBlockData class is to set up
    the construction arguments for the control volume block, automatically link
    to the time domain of the parent block, and to get the information about
    the property and reaction packages.
    """

    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=DefaultBool,
            default=useDefault,
            description="Dynamic model flag",
            doc="""Indicates whether this model will be dynamic,
**default** - useDefault.
**Valid values:** {
**useDefault** - get flag from parent,
**True** - set as a dynamic model,
**False** - set as a steady-state model}""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these, **default** - None. **Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package",
        ConfigValue(
            default=None,
            domain=is_reaction_parameter_block,
            description="Reaction package to use for control volume",
            doc="""Reaction parameter object used to define reaction calculations,
**default** - None.
**Valid values:** {
**None** - no reaction package,
**ReactionParameterBlock** - a ReactionParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "reaction_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing reaction packages",
            doc="""A ConfigBlock with arguments to be passed to a reaction block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see reaction package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "auto_construct",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Argument indicating whether ControlVolume should "
            "automatically construct balance equations",
            doc="""If set to True, this argument will trigger the auto_construct
method which will attempt to construct a set of material, energy and momentum
balance equations based on the parent unit's config block. The parent unit must
have a config block which derives from CONFIG_Base,
**default** - False.
**Valid values:** {
**True** - use automatic construction,
**False** - do not use automatic construction.}""",
        ),
    )

    def build(self):
        """
        General build method for Control Volumes blocks. This method calls a
        number of sub-methods which automate the construction of expected
        attributes of all ControlVolume blocks.

        Inheriting models should call `super().build`.

        Args:
            None

        Returns:
            None
        """
        super(ControlVolumeBlockData, self).build()

        # Setup dynamics flag and time domain
        self._setup_dynamics()

        # Get property package details
        self._get_property_package()

        # Get indexing sets
        self._get_indexing_sets()

        # Get reaction package details (as necessary)
        self._get_reaction_package()

        if self.config.auto_construct is True:
            self._auto_construct()

    def add_geometry(self, *args, **kwargs):
        """
        Method for defining the geometry of the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_geometry. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_material_balances(
        self, balance_type=MaterialBalanceType.useDefault, **kwargs
    ):
        """
        General method for adding material balances to a control volume.
        This method makes calls to specialised sub-methods for each type of
        material balance.

        Args:
            balance_type - MaterialBalanceType Enum indicating which type of
                    material balance should be constructed.
            has_rate_reactions - whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions - whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium - whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer - whether generic mass transfer terms should be
                    included in material balances
            custom_molar_term - a Pyomo Expression representing custom terms to
                    be included in material balances on a molar basis.
            custom_mass_term - a Pyomo Expression representing custom terms to
                    be included in material balances on a mass basis.

        Returns:
            Constraint objects constructed by sub-method
        """
        # Check if balance_type is useDefault, and get default if necessary
        if balance_type == MaterialBalanceType.useDefault:
            try:
                blk = self._get_representative_property_block()
                balance_type = blk.default_material_balance_type()
            except NotImplementedError:
                raise ConfigurationError(
                    "{} property package has not implemented a "
                    "default_material_balance_type, thus cannot use "
                    "MaterialBalanceType.useDefault when constructing "
                    "material balances. Please contact the developer of "
                    "your property package to implement the necessary "
                    "default attributes.".format(self.name)
                )

        self._constructed_material_balance_type = balance_type
        if balance_type == MaterialBalanceType.none:
            mb = None
        elif balance_type == MaterialBalanceType.componentPhase:
            mb = self.add_phase_component_balances(**kwargs)
        elif balance_type == MaterialBalanceType.componentTotal:
            mb = self.add_total_component_balances(**kwargs)
        elif balance_type == MaterialBalanceType.elementTotal:
            mb = self.add_total_element_balances(**kwargs)
        elif balance_type == MaterialBalanceType.total:
            mb = self.add_total_material_balances(**kwargs)
        else:
            raise ConfigurationError(
                "{} invalid balance_type for add_material_balances."
                "Please contact the unit model developer with this bug.".format(
                    self.name
                )
            )

        return mb

    def add_energy_balances(self, balance_type=EnergyBalanceType.useDefault, **kwargs):
        """
        General method for adding energy balances to a control volume.
        This method makes calls to specialised sub-methods for each type of
        energy balance.

        Args:
            balance_type (EnergyBalanceType): Enum indicating which type of
                energy balance should be constructed.
            has_heat_of_reaction (bool): whether terms for heat of reaction
                should be included in energy balance
            has_heat_transfer (bool): whether generic heat transfer terms
                should be included in energy balances
            has_work_transfer (bool): whether generic mass transfer terms
                should be included in energy balances
            has_enthalpy_transfer (bool): whether generic enthalpy transfer
                terms should be included in energy balances
            custom_term (Expression): a Pyomo Expression representing custom
                terms to be included in energy balances

        Returns:
            Constraint objects constructed by sub-method
        """
        # Check if balance_type is useDefault, and get default if necessary
        if balance_type == EnergyBalanceType.useDefault:
            try:
                blk = self._get_representative_property_block()
                balance_type = blk.default_energy_balance_type()
            except NotImplementedError:
                raise ConfigurationError(
                    "{} property package has not implemented a "
                    "default_energy_balance_type, thus cannot use "
                    "EnergyBalanceType.useDefault when constructing "
                    "energy balances. Please contact the developer of "
                    "your property package to implement the necessary "
                    "default attributes.".format(self.name)
                )

        self._constructed_energy_balance_type = balance_type
        if balance_type == EnergyBalanceType.none:
            eb = None
        elif balance_type == EnergyBalanceType.enthalpyTotal:
            eb = self.add_total_enthalpy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.enthalpyPhase:
            eb = self.add_phase_enthalpy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.energyTotal:
            eb = self.add_total_energy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.energyPhase:
            eb = self.add_phase_energy_balances(**kwargs)
        elif balance_type == EnergyBalanceType.isothermal:
            eb = self.add_isothermal_constraint(**kwargs)
        else:
            raise ConfigurationError(
                "{} invalid balance_type for add_energy_balances."
                "Please contact the unit model developer with this bug.".format(
                    self.name
                )
            )

        return eb

    def add_momentum_balances(
        self, balance_type=MomentumBalanceType.pressureTotal, **kwargs
    ):
        """
        General method for adding momentum balances to a control volume.
        This method makes calls to specialised sub-methods for each type of
        momentum balance.

        Args:
            balance_type (MomentumBalanceType): Enum indicating which type of
                momentum balance should be constructed. Default =
                MomentumBalanceType.pressureTotal.
            has_pressure_change (bool): whether default generation terms for
                pressure change should be included in momentum balances
            custom_term (Expression): a Pyomo Expression representing custom
                terms to be included in momentum balances

        Returns:
            Constraint objects constructed by sub-method
        """
        self._constructed_momentum_balance_type = balance_type
        if balance_type == MomentumBalanceType.none:
            mb = None
        elif balance_type == MomentumBalanceType.pressureTotal:
            mb = self.add_total_pressure_balances(**kwargs)
        elif balance_type == MomentumBalanceType.pressurePhase:
            mb = self.add_phase_pressure_balances(**kwargs)
        elif balance_type == MomentumBalanceType.momentumTotal:
            mb = self.add_total_momentum_balances(**kwargs)
        elif balance_type == MomentumBalanceType.momentumPhase:
            mb = self.add_phase_momentum_balances(**kwargs)
        else:
            raise ConfigurationError(
                "{} invalid balance_type for add_momentum_balances."
                "Please contact the unit model developer with this bug.".format(
                    self.name
                )
            )

        return mb

    def _auto_construct(self):
        """
        Placeholder _auto_construct method to ensure a useful exception is
        returned if auto_build is set to True but something breaks in the
        process. Derived ControlVolume classes should overload this.

        Args:
            None

        Returns:
            None
        """
        parent = self.parent_block()

        self.add_geometry()
        self.add_state_blocks()
        self.add_reaction_blocks()

        self.add_material_balances(
            material_balance_type=parent.config.material_balance_type,
            has_rate_reactions=parent.config.has_rate_reactions,
            has_equilibrium_reactions=parent.config.has_equilibrium_reactions,
            has_phase_equilibrium=parent.config.has_phase_equilibrium,
            has_mass_transfer=parent.config.has_mass_transfer,
        )

        self.add_energy_balances(
            energy_balance_type=parent.config.energy_balance_type,
            has_heat_of_reaction=parent.config.has_heat_of_reaction,
            has_heat_transfer=parent.config.has_heat_transfer,
            has_work_transfer=parent.config.has_work_transfer,
            has_enthalpy_transfer=parent.config.has_enthalpy_transfer,
        )

        self.add_momentum_balances(
            has_pressure_change=parent.config.has_pressure_change
        )

        try:
            self.apply_transformation()
        except AttributeError:
            pass

    # Add placeholder methods for adding property and reaction packages
    def add_state_blocks(self, *args, **kwargs):
        """
        Method for adding StateBlocks to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_state_blocks. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_reaction_blocks(self, *args, **kwargs):
        """
        Method for adding ReactionBlocks to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_reaction_blocks. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    # Add placeholder methods for all types of material, energy and momentum
    # balance equations which return NotImplementedErrors
    def add_phase_component_balances(self, *args, **kwargs):
        """
        Method for adding material balances indexed by phase and component to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_component_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_component_balances(self, *args, **kwargs):
        """
        Method for adding material balances indexed by component to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_component_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_element_balances(self, *args, **kwargs):
        """
        Method for adding total elemental material balances indexed to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_element_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_material_balances(self, *args, **kwargs):
        """
        Method for adding a total material balance to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_material_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_enthalpy_balances(self, *args, **kwargs):
        """
        Method for adding enthalpy balances indexed by phase to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_enthalpy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_enthalpy_balances(self, *args, **kwargs):
        """
        Method for adding a total enthalpy balance to
        the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_enthalpy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_energy_balances(self, *args, **kwargs):
        """
        Method for adding energy balances (including kinetic energy) indexed by
        phase to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_energy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_energy_balances(self, *args, **kwargs):
        """
        Method for adding a total energy balance (including kinetic energy)
        to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_energy_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_isothermal_constraint(self, *args, **kwargs):
        """
        Method for adding an isothermal constraint to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            f"{self.name} control volume class has not implemented a method for "
            "add_isothermal_constraint. Please contact the "
            "developer of the ControlVolume class you are using."
        )

    def add_phase_pressure_balances(self, *args, **kwargs):
        """
        Method for adding pressure balances indexed by
        phase to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_pressure_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_pressure_balances(self, *args, **kwargs):
        """
        Method for adding a total pressure balance to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_pressure_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_phase_momentum_balances(self, *args, **kwargs):
        """
        Method for adding momentum balances indexed by phase to the control
        volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_phase_momentum_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def add_total_momentum_balances(self, *args, **kwargs):
        """
        Method for adding a total momentum balance to the control volume.

        See specific control volume documentation for details.
        """
        raise NotImplementedError(
            "{} control volume class has not implemented a method for "
            "add_total_momentum_balances. Please contact the "
            "developer of the ControlVolume class you are using.".format(self.name)
        )

    def _rxn_rate_conv(b, t, x, j):
        """
        Method to determine conversion term for reaction rate terms in material
        balance equations. This method gets the basis of the material flow
        and reaction rate terms and determines the correct conversion factor.
        """
        if x is None:
            # 0D control volume
            flow_basis = b.properties_out[t].get_material_flow_basis()
            prop = b.properties_out[t]
            rxn_basis = b.reactions[t].get_reaction_rate_basis()
        else:
            # 1D control volume
            flow_basis = b.properties[t, x].get_material_flow_basis()
            prop = b.properties[t, x]
            rxn_basis = b.reactions[t, x].get_reaction_rate_basis()

        # Check for undefined basis
        if flow_basis == MaterialFlowBasis.other:
            raise ConfigurationError(
                "{} contains reaction terms, but the property package "
                "used an undefined basis (MaterialFlowBasis.other). "
                "Rate based reaction terms require the property "
                "package to define the basis of the material flow "
                "terms.".format(b.name)
            )
        if rxn_basis == MaterialFlowBasis.other:
            raise ConfigurationError(
                "{} contains reaction terms, but the reaction package "
                "used an undefined basis (MaterialFlowBasis.other). "
                "Rate based reaction terms require the reaction "
                "package to define the basis of the reaction rate "
                "terms.".format(b.name)
            )

        try:
            if flow_basis == rxn_basis:
                return 1
            elif (
                flow_basis == MaterialFlowBasis.mass
                and rxn_basis == MaterialFlowBasis.molar
            ):
                return prop.mw_comp[j]
            elif (
                flow_basis == MaterialFlowBasis.molar
                and rxn_basis == MaterialFlowBasis.mass
            ):
                return 1 / prop.mw_comp[j]
            else:
                raise BurntToast(
                    "{} encountered unrecognized combination of bases "
                    "for reaction rate terms. Please contact the IDAES"
                    " developers with this bug.".format(b.name)
                )
        except AttributeError:
            raise PropertyNotSupportedError(
                "{} property package does not support "
                "molecular weight (mw), which is required for "
                "using property and reaction packages with "
                "different bases.".format(b.name)
            )

    def _get_representative_property_block(self):
        try:
            t_ref = self.flowsheet().time.first()
        except AttributeError:
            raise ConfigurationError(
                "{} control volume does not appear to be part of a "
                "flowsheet (could not find a time attribute).".format(self.name)
            )

        try:
            rep_blk = self.properties_out[t_ref]
        except AttributeError:
            try:
                d_ref = self.length_domain.first()
                rep_blk = self.properties[t_ref, d_ref]
            except AttributeError:
                raise BurntToast(
                    "{} Something went wrong when trying to find "
                    "a representative StateBlock. Please contact "
                    "the IDAES developers with this bug.".format(self.name)
                )
        return rep_blk

    def _estimate_next_state(self, state1, state2, index, always_estimate=False):
        """
        Common method to estimate values for state variables in one state based on
        previous state. This method will not change values of fixed variables.
        Works for both 0D and 1D control volumes.

        Args:
            state1 - StateBlockData to use as the source for values
            state2 - StateBlockData on which to set estimated values
            index - index to use for states and other indexed CV level variables
            always_estimate - bool indicating whether method should overwrite existing values
                on unfixed variables

        Returns:
            None

        """
        state_vars = state2.define_state_vars()

        for v2 in state_vars.values():
            v1 = state1.find_component(v2.parent_component().local_name)
            if v2.is_indexed():
                for k in v2.keys():
                    self._estimate_state_var(v1[k], v2[k], index, always_estimate)
            else:
                self._estimate_state_var(v1, v2, index, always_estimate)

    def _estimate_state_var(self, v1, v2, index, always_estimate=False):
        """
        Method to set value of a given state variable (scalar or indexed) based
        on value from another state variable.

        This method contains some logic for incorporating control volume level
        transfer terms into the estimates. Transfer terms supported include deltaT,
        deltaP, heat and work. Material flow and concentration terms do not include
        additional terms due to complexity in determining how to interpret these.

        Args:
            v1 - state variable to use as source for values
            v2 - state variable to set estimate values on
            index - index for getting values from control volume level terms
            always_estimate - bool indicating whether method should overwrite existing values
                on unfixed variables

        Returns:
            None

        """
        if v2.fixed:
            # Do not touch fixed Vars
            pass
        elif v2.value is not None and not always_estimate:
            # Var has a value and we are not always estimating - do nothing
            pass
        else:
            # Estimate value for v2 from v1
            # If no better guess, use v1
            v2val = v1

            # Check for known special cases
            n = v2.local_name
            # Material flows are too hard to deal with generically, as they
            # can be defined as bilinear terms, can have various indexing and
            # can be on different bases.
            # Similarly, mass and mole fractions are hard to deal with.
            # Energy terms also have problems, as they are defined on an intensive
            # basis, so we would need to convert heat and work
            if n.startswith("pressure"):
                # Pressure, see if there is a fixed deltaP
                if hasattr(self, "deltaP") and self.deltaP[index].fixed:
                    v2val = value(v1 + self.deltaP[index])
            elif n.startswith("temperature"):
                # Temperature, see if there is a fixed deltaT
                if hasattr(self, "deltaT") and self.deltaT[index].fixed:
                    v2val = value(v1 + self.deltaT[index])

            v2.set_value(v2val)
