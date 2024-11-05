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
One-dimensional membrane class for CO2 gas separation
"""

# pylint: disable=unused-import
from enum import Enum
from pyomo.common.config import Bool, ConfigDict, ConfigValue, In
from pyomo.environ import (
    Constraint,
    Param,
    Var,
    units,
    Expression,
)
from pyomo.network import Port

from idaes.core import (
    FlowDirection,
    UnitModelBlockData,
    declare_process_block_class,
    useDefault,
    MaterialFlowBasis,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.models.unit_models.mscontactor import MSContactor
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import create_stream_table_dataframe

__author__ = "Maojian Wang"


class MembraneFlowPattern(Enum):
    """
    Enum of supported flow patterns for membrane.
    So far only support countercurrent and cocurrent flow
    """

    COUNTERCURRENT = 1
    COCURRENT = 2


@declare_process_block_class("Membrane1D")
class Membrane1DData(UnitModelBlockData):

    """Standard Membrane 1D Unit Model Class."""

    CONFIG = UnitModelBlockData.CONFIG()

    Stream_Config = ConfigDict()

    Stream_Config.declare(
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
    Stream_Config.declare(
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

    Stream_Config.declare(
        "has_energy_balance",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Bool indicating whether to include energy balance for stream. Default=True.",
        ),
    )
    Stream_Config.declare(
        "has_pressure_balance",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Bool indicating whether to include pressure balance for stream. Default=True.",
        ),
    )

    CONFIG.declare(
        "sweep_flow",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Bool indicating whether there is a sweep flow in the permeate side.",
            description="Bool indicating whether stream has a feed Port and inlet "
            "state, or if all flow is provided via mass transfer. Default=True.",
        ),
    )
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=5,
            domain=int,
            description="Number of finite elements in length domain",
            doc="""Number of finite elements to use when discretizing length
                domain (default=5)""",
        ),
    )
    CONFIG.declare(
        "flow_type",
        ConfigValue(
            default=MembraneFlowPattern.COUNTERCURRENT,
            domain=In(MembraneFlowPattern),
            description="Flow configuration of membrane",
            doc="""Flow configuration of membrane
               - MembraneFlowPattern.COCURRENT: feed and sweep flows from 0 to 1
               - MembraneFlowPattern.COUNTERCURRENT: feed side flows from 0 to 1
                                                    sweep side flows from 1 to 0  (default)""",
        ),
    )

    for side_name in ["feed", "sweep"]:
        CONFIG.declare(
            side_name + "_side",
            Stream_Config(),
        )

    def build(self):
        """
        This is a one-dimensional model for gas separation in CO₂ capture applications.
        The model will be discretized in the flow direction, and it supports two flow patterns:
        counter-current flow and co-current flow. The model was customized for gas-phase separation
        in CO₂ capture with a single-layer design. If a multi-layer design is needed, multiple units
        can be connected for this application. The two sides of the membrane are called the feed side
        and sweep side. The sweep stream inlet is optional. The driving force across the membrane is the
        partial pressure difference in this gas separation application. Additionally, the energy balance
        assumes that temperature remains constant on each side of the membrane.

        """
        super().build()

        feed_dict = dict(self.config.feed_side)
        sweep_dict = dict(self.config.sweep_side)

        feed_dict["flow_direction"] = FlowDirection.forward
        if self.config.flow_type == MembraneFlowPattern.COCURRENT:
            sweep_dict["flow_direction"] = FlowDirection.forward
        elif self.config.flow_type == MembraneFlowPattern.COUNTERCURRENT:
            sweep_dict["flow_direction"] = FlowDirection.backward
        else:
            raise ConfigurationError(
                f"{self.name} Membrane1D only supports cocurrent and "
                "countercurrent flow patterns, but flow_type configuration"
                " argument was set to {config.flow_type}."
            )

        if self.config.sweep_flow is False:
            sweep_dict["has_feed"] = False

        streams_dict = {"feed_side": feed_dict, "sweep_side": sweep_dict}
        self.mscontactor = MSContactor(
            streams=streams_dict,
            number_of_finite_elements=self.config.finite_elements,
        )

        self.feed_side_inlet = Port(extends=self.mscontactor.feed_side_inlet)
        self.feed_side_outlet = Port(extends=self.mscontactor.feed_side_outlet)
        if self.config.sweep_flow is True:
            self.sweep_side_inlet = Port(extends=self.mscontactor.sweep_side_inlet)
        self.sweep_side_outlet = Port(extends=self.mscontactor.sweep_side_outlet)

        self._make_geometry()
        self._make_performance()

    def _make_geometry(self):

        self.area = Var(initialize=100, units=units.cm**2, doc="The membrane area")

        self.length = Var(initialize=100, units=units.cm, doc="The membrane length")
        self.cell_length = Expression(expr=self.length / self.config.finite_elements)

        self.cell_area = Var(initialize=100, units=units.cm**2, doc="The membrane area")

        @self.Constraint()
        def area_per_cell(self):
            return self.cell_area == self.area / self.config.finite_elements

    def _make_performance(self):
        feed_side_units = (
            self.config.feed_side.property_package.get_metadata().derived_units
        )

        self.permeance = Var(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.mscontactor.feed_side.component_list,
            initialize=1,
            doc="Values in Gas Permeance Unit(GPU)",
            units=units.dimensionless,
        )

        self.gpu_factor = Param(
            default=10e-8 / 13333.2239,
            units=units.m / units.s / units.Pa,
            mutable=True,
            doc=" This is a coefficient that will convert the unit of permeability from GPU to SI units for further calculation",
        )

        """
        Selectivity is defined for cases where some permeabilities are unavailable. Only define the selectivity
         between different components; others should be set to 1. If the permeabilities of all components are defined
         there is no need to define selectivity, as this may overdefine the problem.
        """
        self.selectivity = Var(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.mscontactor.feed_side.component_list,
            self.mscontactor.feed_side.component_list,
            initialize=1,
            units=units.dimensionless,
        )

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            self.mscontactor.feed_side.component_list,
            self.mscontactor.feed_side.component_list,
            doc="permeance calculation",
        )
        def permeance_calculation(self, t, e, a, b):
            return (
                self.permeance[t, e, a] * self.selectivity[t, e, a, b]
                == self.permeance[t, e, b]
            )

        p_units = feed_side_units.PRESSURE

        crossover_component_list = list(
            set(self.mscontactor.feed_side.component_list)
            & set(self.mscontactor.sweep_side.component_list)
        )

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            crossover_component_list,
            doc="permeability calculation",
        )
        def permeability_calculation(self, t, s, m):
            feed_side_state = self.mscontactor.feed_side[t, s]
            if feed_side_state.get_material_flow_basis() is MaterialFlowBasis.molar:
                mb_units = feed_side_units.FLOW_MOLE
                rho = self.mscontactor.feed_side[t, s].dens_mol
            elif feed_side_state.get_material_flow_basis() is MaterialFlowBasis.mass:
                mb_units = feed_side_units.FLOW_MASS
                rho = self.mscontactor.feed_side[t, s].dens_mass
            else:
                raise TypeError(
                    "This model only supports MaterialFlowBasis equal to molar or mass"
                )

            return self.mscontactor.material_transfer_term[
                t, s, "feed_side", "sweep_side", m
            ] == -units.convert(
                (
                    rho
                    * self.gpu_factor
                    * self.permeance[t, s, m]
                    * self.cell_area
                    * (
                        self.mscontactor.feed_side[t, s].pressure
                        * self.mscontactor.feed_side[t, s].mole_frac_comp[m]
                        - units.convert(
                            self.mscontactor.sweep_side[t, s].pressure, to_units=p_units
                        )
                        * self.mscontactor.sweep_side[t, s].mole_frac_comp[m]
                    )
                ),
                to_units=mb_units,
            )

        @self.Constraint(
            self.flowsheet().time,
            self.mscontactor.elements,
            doc="Energy balance",
        )
        def energy_transfer(self, t, s):
            return (
                self.mscontactor.feed_side[t, s].temperature
                == self.mscontactor.sweep_side[t, s].temperature
            )

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Feed Inlet": self.feed_side_inlet,
                "Feed Outlet": self.feed_side_outlet,
                "Permeate Inlet": self.sweep_side_inlet,
                "Permeate Outlet": self.sweep_side_outlet,
            },
            time_point=time_point,
        )
