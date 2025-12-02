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
Thickener unit model.

Unit model is derived from:

R. Burger, F. Concha, K.H. Karlsen, A. Narvaez,
Numerical simulation of clarifier-thickener units treating ideal
suspensions with a flux density function having two inflection points,
Mathematical and Computer Modelling 44 (2006) 255–275
doi:10.1016/j.mcm.2005.11.008

Settling velocity function from:

N.G. Barton, C.-H. Li, S.J. Spencer, Control of a surface of discontinuity in continuous thickeners,
Journal of the Australian Mathematical Society Series B 33 (1992) 269–289

"""
# Import Python libraries
import logging
from pandas import DataFrame

# Import Pyomo libraries
from pyomo.environ import Expr_if, inequality, units, Var
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
from idaes.core.util.constants import Constants as CONST


__author__ = "Andrew Lee"


# Set up logger
logger = logging.getLogger("idaes.unit_model")


# TODO Evaluate whether this unit model would be better in WaterTap or PrOMMiS
@declare_process_block_class("Thickener0D")
class Thickener0DData(UnitModelBlockData):
    """
    Thickener0D Unit Model Class
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
        # Call super().build to setup dynamics
        super().build()

        # Build Solid Phase
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.solid_property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.solid_inlet_state = self.config.solid_property_package.build_state_block(
            self.flowsheet().time, doc="Solid properties in separator", **tmp_dict
        )

        # Add solid splitter
        self.solid_split = Separator(
            property_package=self.config.solid_property_package,
            property_package_args=self.config.solid_property_package_args,
            outlet_list=["underflow", "overflow"],
            split_basis=SplittingType.totalFlow,
            ideal_separation=False,
            mixed_state_block=self.solid_inlet_state,
            has_phase_equilibrium=False,
            material_balance_type=self.config.material_balance_type,
            momentum_balance_type=self.config.momentum_balance_type,
            energy_split_basis=self.config.energy_split_basis,
        )

        # Add solid ports
        self.add_port(
            name="solid_inlet",
            block=self.solid_inlet_state,
            doc="Solid inlet to thickener",
        )
        self.solid_underflow = Port(extends=self.solid_split.underflow)
        self.solid_overflow = Port(extends=self.solid_split.overflow)

        # Build liquid Phase
        # Setup StateBlock argument dict
        tmp_dict = dict(**self.config.liquid_property_package_args)
        tmp_dict["has_phase_equilibrium"] = False
        tmp_dict["defined_state"] = True

        self.liquid_inlet_state = self.config.liquid_property_package.build_state_block(
            self.flowsheet().time, doc="liquid properties in separator", **tmp_dict
        )

        # Add liquid splitter
        self.liquid_split = Separator(
            property_package=self.config.liquid_property_package,
            property_package_args=self.config.liquid_property_package_args,
            outlet_list=["underflow", "overflow"],
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
            doc="liquid inlet to thickener",
        )
        self.liquid_underflow = Port(extends=self.liquid_split.underflow)
        self.liquid_overflow = Port(extends=self.liquid_split.overflow)

        # Add additional variables and constraints
        uom = self.solid_inlet_state.params.get_metadata().derived_units

        self.area = Var(
            initialize=1,
            units=uom.AREA,
            doc="Cross sectional area of thickener",
        )

        # Volumetric Flowrates
        self.flow_vol_feed = Var(
            self.flowsheet().time,
            initialize=0.7,
            units=uom.FLOW_VOL,
            bounds=(0, None),
            doc="Total volumetric flowrate of feed",
        )
        self.flow_vol_overflow = Var(
            self.flowsheet().time,
            initialize=0.7,
            units=uom.FLOW_VOL,
            bounds=(0, None),
            doc="Total volumetric flowrate of overflow",
        )
        self.flow_vol_underflow = Var(
            self.flowsheet().time,
            initialize=0.7,
            units=uom.FLOW_VOL,
            bounds=(0, None),
            doc="Total volumetric flowrate of underflow",
        )

        # Solid Fractions
        self.solid_fraction_feed = Var(
            self.flowsheet().time,
            initialize=0.7,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Volume fraction of solids in feed",
        )
        self.solid_fraction_underflow = Var(
            self.flowsheet().time,
            initialize=0.7,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Volume fraction of solids in underflow",
        )
        self.solid_fraction_overflow = Var(
            self.flowsheet().time,
            initialize=0.7,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Volume fraction of solids in overflow",
        )

        # Flux densities
        self.flux_density_underflow = Var(
            self.flowsheet().time,
            initialize=0,
            units=uom.VELOCITY,
            doc="Kynch flux density in underflow",
        )
        self.flux_density_overflow = Var(
            self.flowsheet().time,
            initialize=0,
            units=uom.VELOCITY,
            doc="Kynch flux density in overflow",
        )

        # Parameters
        self.particle_size = Var(
            self.flowsheet().time,
            initialize=1e-5,
            units=uom.LENGTH,
            doc="Characteristic length of particle",
        )
        self.v0 = Var(
            self.flowsheet().time,
            initialize=1e-4,
            units=uom.VELOCITY,
            doc="Stokes velocity of individual particle",
        )
        self.v1 = Var(
            initialize=1e-5,
            units=uom.VELOCITY,
            doc="Superficial velocity of a Darcy type flow of liquid",
        )
        self.C = Var(
            initialize=5,
            units=units.dimensionless,
            bounds=(0, None),
            doc="Settling velocity exponent",
        )
        self.solid_fraction_max = Var(
            initialize=0.9,
            units=units.dimensionless,
            bounds=(0, 1),
            doc="Maximum achievable solids volume fraction",
        )

        # ---------------------------------------------------------------------------------------------
        # Constraints
        @self.Constraint(self.flowsheet().time)
        def feed_flowrate(b, t):
            return b.flow_vol_feed[t] == (
                b.solid_inlet_state[t].flow_vol
                + units.convert(b.liquid_inlet_state[t].flow_vol, to_units=uom.FLOW_VOL)
            )

        @self.Constraint(self.flowsheet().time)
        def overflow_flowrate(b, t):
            return b.flow_vol_overflow[t] == (
                b.solid_split.overflow_state[t].flow_vol
                + units.convert(
                    b.liquid_split.overflow_state[t].flow_vol, to_units=uom.FLOW_VOL
                )
            )

        @self.Constraint(self.flowsheet().time)
        def underflow_flowrate(b, t):
            return b.flow_vol_underflow[t] == (
                b.solid_split.underflow_state[t].flow_vol
                + units.convert(
                    b.liquid_split.underflow_state[t].flow_vol, to_units=uom.FLOW_VOL
                )
            )

        # Eqn 2.8 from [1]
        @self.Constraint(self.flowsheet().time)
        def flux_density_function_overflow(b, t):
            u = b.solid_fraction_overflow
            return b.flux_density_overflow[t] == Expr_if(
                IF=inequality(0, u[t], b.solid_fraction_max),
                THEN=b.v0[t] * u[t] * (1 - u[t] / b.solid_fraction_max) ** b.C
                + b.v1 * u[t] ** 2 * (b.solid_fraction_max - u[t]),
                ELSE=0 * units.m * units.s**-1,
            )

        # Eqn 2.8 from [1]
        @self.Constraint(self.flowsheet().time)
        def flux_density_function_underflow(b, t):
            u = b.solid_fraction_underflow
            return b.flux_density_underflow[t] == Expr_if(
                IF=inequality(0, u[t], b.solid_fraction_max),
                THEN=b.v0[t] * u[t] * (1 - u[t] / b.solid_fraction_max) ** b.C
                + b.v1 * u[t] ** 2 * (b.solid_fraction_max - u[t]),
                ELSE=0 * units.m * units.s**-1,
            )

        # Modified from 2.23 and 2.33 from [1]
        @self.Constraint(self.flowsheet().time)
        def solids_continuity(b, t):
            return b.flow_vol_feed[t] * b.solid_fraction_feed[t] == (
                b.area * (b.flux_density_overflow[t] + b.flux_density_underflow[t])
                - b.flow_vol_overflow[t]
                * (b.solid_fraction_overflow[t] - b.solid_fraction_feed[t])
                + b.flow_vol_underflow[t]
                * (b.solid_fraction_underflow[t] - b.solid_fraction_feed[t])
            )

        @self.Constraint(self.flowsheet().time)
        def solids_conservation(b, t):
            return b.flow_vol_feed[t] * b.solid_fraction_feed[t] == (
                +b.flow_vol_overflow[t] * b.solid_fraction_overflow[t]
                + b.flow_vol_underflow[t] * b.solid_fraction_underflow[t]
            )

        @self.Constraint(self.flowsheet().time)
        def maximum_underflow_volume_fraction(b, t):
            return b.solid_fraction_underflow[t] <= b.solid_fraction_max

        @self.Constraint(self.flowsheet().time)
        def maximum_overflow_volume_fraction(b, t):
            return b.solid_fraction_overflow[t] <= b.solid_fraction_max

        @self.Constraint(self.flowsheet().time)
        def inlet_volume_fraction(b, t):
            return b.solid_inlet_state[t].flow_vol == (
                b.solid_fraction_feed[t]
                * (
                    b.solid_inlet_state[t].flow_vol
                    + units.convert(
                        b.liquid_inlet_state[t].flow_vol, to_units=uom.FLOW_VOL
                    )
                )
            )

        @self.Constraint(self.flowsheet().time)
        def underflow_volume_fraction(b, t):
            return b.solid_split.underflow_state[t].flow_vol == (
                b.solid_fraction_underflow[t]
                * (
                    b.solid_split.underflow_state[t].flow_vol
                    + units.convert(
                        b.liquid_split.underflow_state[t].flow_vol,
                        to_units=uom.FLOW_VOL,
                    )
                )
            )

        @self.Constraint(self.flowsheet().time)
        def stokes_law(b, t):
            # Assuming constant properties, source from feed states
            return 18 * b.v0[t] * units.convert(
                b.liquid_inlet_state[t].visc_d, to_units=uom.DYNAMIC_VISCOSITY
            ) == (
                (
                    b.solid_inlet_state[t].dens_mass
                    - units.convert(
                        b.liquid_inlet_state[t].dens_mass, to_units=uom.DENSITY_MASS
                    )
                )
                * units.convert(CONST.acceleration_gravity, to_units=uom.ACCELERATION)
                * b.particle_size[t] ** 2
            )

    def _get_performance_contents(self, time_point=0):
        return {
            "vars": {
                "Area": self.area,
                "Liquid Recovery": self.liquid_split.split_fraction[
                    time_point, "overflow"
                ],
                "Feed Solid Fraction": self.solid_fraction_feed[time_point],
                "Underflow Solid Fraction": self.solid_fraction_underflow[time_point],
                "particle size": self.particle_size[time_point],
                "v0": self.v0[time_point],
                "v1": self.v1,
                "C": self.C,
                "solid_fraction_max": self.solid_fraction_max,
            },
        }

    def _get_stream_table_contents(self, time_point=0):
        stream_attributes = {}
        stream_attributes["Units"] = {}

        sblocks = {
            "Feed Solid": self.solid_inlet_state,
            "Feed Liquid": self.liquid_inlet_state,
            "Underflow Solid": self.solid_split.underflow_state,
            "Underflow Liquid": self.liquid_split.underflow_state,
            "Overflow Solid": self.solid_split.overflow_state,
            "Overflow Liquid": self.liquid_split.overflow_state,
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

    def initialize(self, **kwargs):
        raise NotImplementedError(
            "The Thickener0D unit model does not support the old initialization API. "
            "Please use the new API (InitializerObjects) instead."
        )
