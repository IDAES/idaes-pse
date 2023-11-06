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
Base model for solid-liquid separations.

This model is intended to form the basis for solid-liquid separations where a mixed
stream enters the unit and is separated into a concentrated solid-liquid stream and
a stream of pure liquid (e.g. filters, thickeners, etc.)

This model assumes separate property packages and ports for the solid and liquid
streams and thus two inlet Ports (solid and liquid) and three outlet Ports (solids,
liquid with solids, separated liquids).

"""
# Import Python libraries
import logging
from pandas import DataFrame

# Import Pyomo libraries
from pyomo.environ import Reference
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.models.unit_models.separator import (
    Separator,
    SplittingType,
    EnergySplittingType,
)
from idaes.core.initialization import BlockTriangularizationInitializer
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.units_of_measurement import report_quantity


__author__ = "Andrew Lee"


# Set up logger
logger = logging.getLogger("idaes.unit_model")


@declare_process_block_class("SLSeparator")
class SLSeparatorData(UnitModelBlockData):
    """
    Standard Solid-Liquid Separator Unit Model Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Flash units do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Flash units do not have defined volume, thus
this must be False.""",
        ),
    )
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
        "energy_split_basis",
        ConfigValue(
            default=EnergySplittingType.equal_temperature,
            domain=EnergySplittingType,
            description="Type of constraint to write for energy splitting",
            doc="""Argument indicating basis to use for splitting energy this is
not used for when ideal_separation == True.
**default** - EnergySplittingType.equal_temperature.
**Valid values:** {
**EnergySplittingType.equal_temperature** - outlet temperatures equal inlet
**EnergySplittingType.equal_molar_enthalpy** - outlet molar enthalpies equal
inlet,
**EnergySplittingType.enthalpy_split** - apply split fractions to enthalpy
flows.}""",
        ),
    )
    CONFIG.declare(
        "solid_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for solid phase",
            doc="""Property parameter object used to define solid phase property
calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "solid_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing solid phase property packages",
            doc="""A ConfigBlock with arguments to be passed to a solid phase property
block(s) and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "liquid_property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for liquid phase",
            doc="""Property parameter object used to define liquid phase property
    calculations,
    **default** - useDefault.
    **Valid values:** {
    **useDefault** - use default package from parent model or flowsheet,
    **PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "liquid_property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing liquid phase property packages",
            doc="""A ConfigBlock with arguments to be passed to a liquid phase property
    block(s) and used when constructing these,
    **default** - None.
    **Valid values:** {
    see property package for documentation.}""",
        ),
    )

    default_initializer = BlockTriangularizationInitializer

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build()

        # Build Solid Phase - inlet=outlet so only need a StateBlock
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.solid_property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.solid_state = self.config.solid_property_package.build_state_block(
            self.flowsheet().time, doc="Solid properties in separator", **tmp_dict
        )

        # Add Solids Ports
        self.add_port(
            name="solid_inlet", block=self.solid_state, doc="Solid inlet to separator"
        )
        self.add_port(
            name="solid_outlet",
            block=self.solid_state,
            doc="Solid outlet from separator",
        )

        # Add Liquid Inlet
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.liquid_property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.liquid_inlet_state = self.config.liquid_property_package.build_state_block(
            self.flowsheet().time,
            doc="Liquid properties at inlet to separator",
            **tmp_dict,
        )

        self.split = Separator(
            property_package=self.config.liquid_property_package,
            property_package_args=self.config.liquid_property_package_args,
            outlet_list=["recovered", "retained"],
            split_basis=SplittingType.totalFlow,
            ideal_separation=False,
            mixed_state_block=self.liquid_inlet_state,
            has_phase_equilibrium=False,
            material_balance_type=self.config.material_balance_type,
            momentum_balance_type=self.config.momentum_balance_type,
            energy_split_basis=self.config.energy_split_basis,
        )

        # Add liquid ports
        self.add_port(
            name="liquid_inlet",
            block=self.liquid_inlet_state,
            doc="Liquid inlet to separator",
        )
        self.recovered_liquid_outlet = Port(extends=self.split.recovered)
        self.retained_liquid_outlet = Port(extends=self.split.retained)

        # Add liquid recovery
        self.liquid_recovery = Reference(self.split.split_fraction[:, "recovered"])

    def initialize(self, **kwargs):
        raise NotImplementedError(
            "The SLSeparator unit model does not support the old initialization API. "
            "Please use the new API (InitializerObjects) instead."
        )

    def _get_performance_contents(self, time_point=0):

        return {"vars": {"Liquid Recovery": self.liquid_recovery[time_point]}}

    def _get_stream_table_contents(self, time_point=0):
        stream_attributes = {}
        stream_attributes["Units"] = {}

        sblocks = {
            "Solid Inlet": self.solid_state,
            "Liquid Inlet": self.liquid_inlet_state,
            "Solid Outlet": self.solid_state,
            "Liquid in Solids Outlet": self.split.retained_state,
            "Recovered Liquid Outlet": self.split.recovered_state,
        }

        for n, v in sblocks.items():
            dvars = v[time_point].define_display_vars()

            stream_attributes[n] = {}

            for k in dvars:
                for i in dvars[k].keys():
                    stream_key = k if i is None else f"{k} {i}"

                    quant = report_quantity(dvars[k][i])

                    stream_attributes[n][stream_key] = quant.m
                    stream_attributes["Units"][stream_key] = quant.u

        return DataFrame.from_dict(stream_attributes, orient="columns")
