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
from pyomo.common.config import ConfigValue
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.var import Var
from pyomo.core.base.param import Param
from pyomo.core.base.expression import Expression
from pyomo.core.base.units_container import units as pyunits
from pyomo.core.base.reference import Reference
from pyomo.dae.contset import ContinuousSet
from pyomo.dae.diffvar import DerivativeVar
from pyomo.core.expr.current import log10, sqrt
from pyomo.core.expr.numvalue import value as pyo_value

from idaes.core.unit_model import UnitModelBlockData
from idaes.core.process_block import declare_process_block_class
from idaes.core.control_volume_base import MaterialBalanceType
from idaes.core.property_base import StateBlock
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_transformation_method,
    is_transformation_scheme,
)
from idaes.core.util.constants import Constants
from idaes.core.control_volume1d import ControlVolume1DBlock

"""
A simple pipeline unit model.
"""


@declare_process_block_class("GasPipeline")
class GasPipelineData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(default=None, domain=is_physical_parameter_block),
    )
    CONFIG.declare("transformation_method", ConfigValue(
        default="dae.finite_difference",
        domain=is_transformation_method,
        description="DAE transformation method",
        doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory."""))
    CONFIG.declare("transformation_scheme", ConfigValue(
        default="FORWARD",
        domain=is_transformation_scheme,
        description="DAE transformation scheme",
        doc="""Scheme to use when transforming domain. See Pyomo
documentation for supported schemes."""))
    CONFIG.declare("finite_elements", ConfigValue(
        default=1,
        domain=int,
        description="Number of finite elements",
        doc="""Number of finite elements to use in transformation (equivalent
to Pyomo nfe argument)."""))
    CONFIG.declare("collocation_points", ConfigValue(
        default=None,
        domain=int,
        description="Number of collocation points",
        doc="""Number of collocation points to use (equivalent to Pyomo ncp
argument)."""))

    def build(self):
        super(GasPipelineData, self).build()

        # self.config is the ConfigBlock "instantiated" from self.CONFIG
        # in ProcessBlockData.
        config = self.config
        time = self.flowsheet().time
        property_package = config.property_package

        if len(property_package.phase_list) != 1:
            raise ValueError(
                "%s can only be constructed with a "
                "single-phase property package.\n"
                "Got phases %s."
                % (self.__class__, [p for p in property_package.phase_list])
            )
        self.phase = next(iter(property_package.phase_list))
        if self.phase != "Vap":
            raise ValueError(
                "%s can only be constructed with a single phase, \"Vap\"."
                "Got phase %s."
                % (self.__class__, self.phase)
            )

        #
        # Make sure property package has the properties that this unit
        # model requires.
        #
        property_dict = property_package.get_metadata().properties
        if "pressure" not in property_dict:
            raise ValueError(
                "Property package supplied to pipeline must have a property "
                "for 'pressure', which was not found in %s."
                % type(property_package)
            )

        #
        # Set up control volume:
        #
        #cv_config = {
        #    "transformation_method": config.transformation_method,
        #    "transformation_scheme": config.transformation_scheme,
        #    "finite_elements": config.finite_elements,
        #    "collocation_points": config.collocation_points,
        #    "has_holdup": config.has_holdup,
        #}
        cv_config = config()
        self.control_volume = ControlVolume1DBlock(default=cv_config)
        self.control_volume.add_geometry()
        self.control_volume.add_state_blocks(
            has_phase_equilibrium=False,
        )
        self.control_volume.add_phase_component_balances()

        cv = self.control_volume
        self.add_diameter()
        self.add_friction_factor()

        #
        # Add constraints indexed by space before applying transformation
        #
        self.add_flow_mass_linking_constraint()
        self.add_momentum_balance_equation()

        #
        # Apply spatial discretization transformation
        #
        self.control_volume.apply_transformation()

        # NOTE: Later we may try to support switches in direction of flow.
        # At this time these may need to be renamed...
        x0 = self.control_volume.length_domain.first()
        xf = self.control_volume.length_domain.last()
        inlet_state = Reference(
            self.control_volume.properties[:, x0],
            ctype=StateBlock,
        )
        outlet_state = Reference(
            self.control_volume.properties[:, xf],
            ctype=StateBlock,
        )
        # I want my ports to be indexed only by time. The natural solution
        # is to build time-indexed state blocks called inlet_state and
        # outlet_state, build ports from these blocks, and connect their
        # state variables to my control volume state blocks with equality
        # constraints. I would rather not add these equality constraints,
        # however, so I use references.
        self.add_port(
            name="inlet_port",
            block=inlet_state,
            doc="The inlet to the pipeline",
        )
        self.add_port(
            name="outlet_port",
            block=outlet_state,
            doc="The outlet from the pipeline",
        )

        # TODO: How should the energy balance in this pipeline be handled?
        self.add_isothermal_constraint()

    def add_friction_factor(self, rugosity=None):
        if rugosity is None:
            rugosity = 0.025 * pyunits.mm
        rug_value = pyo_value(rugosity)
        rug_units = pyunits.get_units(rugosity)
        self.rugosity = Param(
            initialize=rug_value,
            units=rug_units,
            #doc=TODO,
            mutable=True,
        )

        diam_mm = pyunits.convert(self.diameter, pyunits.mm)
        friction_factor_expr = (2*log10(3.7*diam_mm/self.rugosity))**-2

        self.friction_factor = Var(
            initialize=pyo_value(friction_factor_expr),
            units=pyunits.dimensionless,
            doc="This is fraction factor lambda_l in V. Zavala's paper",
        )

        self.friction_factor_eqn = Constraint(expr=(
            self.friction_factor == friction_factor_expr
        ))

    def add_diameter(self):
        # We make diameter a variable as it may be more convenient
        # to fix diameter than to fix area.
        diam_units = pyunits.m
        self.diameter = Var(
            initialize=1.0,
            units=diam_units,
            doc="Diameter of the pipeline",
        )

        area = self.control_volume.area
        def diameter_rule(b):
            return (
                pyunits.convert(area, diam_units**2) 
                == Constants.pi * self.diameter**2 / 4.0
            )
        self.diameter_eqn = Constraint(
            rule=diameter_rule,
            doc="Equation linking control volume area and pipeline diameter",
        )

    def add_flow_mass_linking_constraint(self):
        # We know flow_mass exists in the property package, but it might be
        # an expression. We need variables to construct the DerivativeVar.
        time = self.flowsheet().time
        space = self.control_volume.length_domain
        kghr = pyunits.kg/pyunits.hr
        cv = self.control_volume
        cv.flow_mass = Var(
            time, space, units=kghr, doc="Mass flow rate through the pipeline"
        )
        state = cv.properties
        def flow_mass_linking_rule(b, t, x):
            return (
                b.flow_mass[t, x]
                - pyunits.convert(state[t, x].flow_mass, kghr)
                == 0
            )
        cv.flow_mass_linking_constraint = Constraint(
            time, space, rule=flow_mass_linking_rule
        )

    def get_friction_term(self, t, x):
        """
        Equation (2.5b) in Zavala (2014)
        """
        cv = self.control_volume
        state = cv.properties[t, x]
        area = cv.area
        flow = state.flow_mass
        pressure = state.pressure
        diameter = self.diameter
        lam = self.friction_factor
        #nu = self.speed_of_sound[t, x]
        nu = state.speed_of_sound
        friction_term = (
            8 * lam * nu**2 / Constants.pi**2 / diameter**5
            * flow * abs(flow) / pressure
        )
        return friction_term

    def add_flow_mass_dt(self):
        time = self.flowsheet().time
        # Should I attach this variable to the unit model or the
        # control volume?
        kghr2 = pyunits.kg / pyunits.hr**2
        cv = self.control_volume
        cv.flow_mass_dt = DerivativeVar(
            cv.flow_mass,
            units=kghr2,
            wrt=time,
        )

    def add_pressure_dx(self):
        # We add pressure_dx here rather than using the control volume
        # add_total_pressure_balances method because we don't want to add
        # a deltaP variable and pressure_balance equation.
        cv = self.control_volume
        cv.pressure = Reference(cv.properties[:, :].pressure)
        pressure_data = next(iter(cv.pressure.values()))
        space = cv.length_domain
        # Refrain from using metadata.derived_units here as I want the
        # units of my pressure variable (bar, hopefully) to override
        # the amalgamation of whatever base units I happen to use.
        p_units = pressure_data.get_units()
        cv.pressure_dx = DerivativeVar(
            cv.pressure,
            # Units are same as pressure because we assume length domain
            # is normalized.
            units=p_units,
            wrt=space,
        )

    def add_momentum_balance_equation(self):
        # This is Equation 2.5 in the paper
        time = self.flowsheet().time
        dynamic = self.flowsheet().config.dynamic
        space = self.control_volume.length_domain
        if dynamic:
            self.add_flow_mass_dt()
        self.add_pressure_dx()
        cv = self.control_volume
        # TODO: Units of these equations should probably not
        # be hard-coded. These should be configurable.
        kgm2hr2 = pyunits.kg / pyunits.m**2 / pyunits.hr**2
        def momentum_balance_rule(b, t, x):
            # TODO: Should probably avoid having unit conversion calls
            # inside this function, as it will be called for every time/space
            # index.
            if dynamic:
                accum_expr = pyunits.convert(
                    cv.flow_mass_dt[t, x] / cv.area,
                    kgm2hr2,
                )
            else:
                accum_expr = 0.0
            flux_expr = pyunits.convert(
                cv.pressure_dx[t, x] / cv.length, kgm2hr2
            )
            friction_expr = pyunits.convert(
                self.get_friction_term(t, x), kgm2hr2
            )
            return accum_expr + flux_expr + friction_expr == 0
        cv.momentum_balance = Constraint(
            time,
            space,
            rule=momentum_balance_rule,
        )

    def add_isothermal_constraint(self):
        time = self.flowsheet().time
        length = self.control_volume.length_domain
        x0 = length.first()
        state = self.control_volume.properties
        def isothermal_rule(b, t, x):
            # NOTE: This constraint makes the pipeline non-templatizable
            # with respect to length. Also, we might want to skip at
            # xf depending on the direction of flow?
            if x != x0:
                x_prev = length.prev(x)
                return state[t, x].temperature == state[t, x_prev].temperature
            else:
                return Constraint.Skip
        self.state_isothermal_eqn = Constraint(
            time, length, rule=isothermal_rule
        )
