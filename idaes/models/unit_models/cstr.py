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
Standard IDAES CSTR model.
"""

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.environ import Reference, Constraint

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
)
from idaes.core.scaling import CustomScalerBase

__author__ = "Andrew Lee, Vibhav Dabadghao, Brandon Paul"


class CSTRScaler(CustomScalerBase):
    """
    Default modular scaler for CSTR unit model.

    This Scaler relies on modular the associated property and reaction packages,
    either through user provided options (submodel_scalers argument) or by default
    Scalers assigned to the packages.

    Reaction generation terms are scaled based on component flow rates, whilst
    extents of reaction are unscaled. Heat duty is scaled to kW and pressure drop
    to 0.1 bar. All constraints are scaled using the inverse maximum scheme.
    """

    DEFAULT_SCALING_FACTORS = {
        # "QuantityName: (reference units, scaling factor)
        "deltaP": 1e-5,
        "volume": 1e3,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to variables in model.

        Submodel Scalers are called for the property and reaction blocks.
        Reaction generation terms are scaled based on component flow rates, whilst
        extents of reaction are unscaled. Heat duty is scaled to kW and pressure drop
        to 0.1 bar.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            None
        """
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            model.control_volume,
            method="variable_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scaling control volume variables

        # Pressure drop
        if hasattr(model.control_volume, "deltaP"):
            for t in model.flowsheet().time:
                self.scale_variable_by_default(
                    model.control_volume.deltaP[t], overwrite=overwrite
                )

        # Heat transfer - optional
        # Scale heat based on enthalpy flow entering reactor
        if hasattr(model.control_volume, "heat"):
            for t in model.flowsheet().time:
                h_in = 0
                for p in model.control_volume.properties_in.phase_list:
                    # The expression for enthalpy flow might include multiple terms,
                    # so we will sum over all the terms provided
                    h_in += sum(
                        self.get_expression_nominal_values(
                            model.control_volume.properties_in[
                                t
                            ].get_enthalpy_flow_terms(p)
                        )
                    )
                # Scale for heat is general one order of magnitude less than enthalpy flow
                self.set_variable_scaling_factor(
                    model.control_volume.heat[t], abs(1 / (0.1 * h_in))
                )

        # Volume
        if hasattr(model.control_volume, "volume"):
            for t in model.flowsheet().time:
                # self.scale_variable_by_default(
                #     model.control_volume.volume[t], overwrite=overwrite
                # )
                self.set_variable_scaling_factor(model.control_volume.volume[t], 1e3)

        # Reaction rate
        for t in model.flowsheet().time:
            if hasattr(model.control_volume.reactions[t], "reaction_rate"):
                for v in model.control_volume.reactions[t].reaction_rate.values():
                    self.set_variable_scaling_factor(v, 1e-3)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to constraints in model.

        Submodel Scalers are called for the property and reaction blocks. All other constraints
        are scaled using the inverse maximum shceme.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            None
        """
        # Call scaling methods for sub-models
        self.call_submodel_scaler_method(
            model.control_volume,
            method="constraint_scaling_routine",
            submodel_scalers=submodel_scalers,
            overwrite=overwrite,
        )

        # Scale control volume constraints
        for c in model.control_volume.component_data_objects(
            Constraint, descend_into=False
        ):
            self.scale_constraint_by_nominal_value(
                c,
                scheme="inverse_maximum",
                overwrite=overwrite,
            )

        # Scale unit level constraints
        if hasattr(model, "cstr_performance_eqn"):
            for c in model.cstr_performance_eqn.values():
                self.scale_constraint_by_nominal_value(
                    c,
                    scheme="inverse_maximum",
                    overwrite=overwrite,
                )


@declare_process_block_class("CSTR")
class CSTRData(UnitModelBlockData):
    """
    Standard CSTR Unit Model Class
    """

    default_scaler = CSTRScaler

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
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
    CONFIG.declare(
        "has_equilibrium_reactions",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Equilibrium reaction construction flag",
            doc="""Indicates whether terms for equilibrium controlled reactions
should be constructed,
**default** - True.
**Valid values:** {
**True** - include equilibrium reaction terms,
**False** - exclude equilibrium reaction terms.}""",
        ),
    )
    CONFIG.declare(
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
    CONFIG.declare(
        "has_heat_of_reaction",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat of reaction term construction flag",
            doc="""Indicates whether terms for heat of reaction terms should be
constructed,
**default** - False.
**Valid values:** {
**True** - include heat of reaction terms,
**False** - exclude heat of reaction terms.}""",
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
**PhysicalParameterObject** - a PhysicalParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
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

    def build(self):
        """
        Begin building model (pre-DAE transformation).
        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(CSTRData, self).build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            reaction_package=self.config.reaction_package,
            reaction_package_args=self.config.reaction_package_args,
        )

        self.control_volume.add_geometry()

        self.control_volume.add_state_blocks(
            has_phase_equilibrium=self.config.has_phase_equilibrium
        )

        self.control_volume.add_reaction_blocks(
            has_equilibrium=self.config.has_equilibrium_reactions
        )

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_rate_reactions=True,
            has_equilibrium_reactions=self.config.has_equilibrium_reactions,
            has_phase_equilibrium=self.config.has_phase_equilibrium,
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_of_reaction=self.config.has_heat_of_reaction,
            has_heat_transfer=self.config.has_heat_transfer,
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Add object references
        self.volume = Reference(self.control_volume.volume[:])

        # Add CSTR performance equation
        @self.Constraint(
            self.flowsheet().time,
            self.config.reaction_package.rate_reaction_idx,
            doc="CSTR performance equation",
        )
        def cstr_performance_eqn(b, t, r):
            return b.control_volume.rate_reaction_extent[t, r] == (
                b.volume[t] * b.control_volume.reactions[t].reaction_rate[r]
            )

        # Set references to balance terms at unit level
        if (
            self.config.has_heat_transfer is True
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            self.heat_duty = Reference(self.control_volume.heat[:])

        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.control_volume.deltaP[:])

    def _get_performance_contents(self, time_point=0):
        var_dict = {"Volume": self.volume[time_point]}
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}
