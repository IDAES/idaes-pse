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
from pyomo.core.base.units_container import units as pyunits
from pyomo.core.base.reference import Reference
from pyomo.dae.diffvar import DerivativeVar
from pyomo.core.expr.current import log10
from pyomo.core.expr.numvalue import value as pyo_value

from idaes.core import (
    declare_process_block_class,
    ControlVolume1DBlock,
    StateBlock,
    UnitModelBlockData,
)
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_transformation_method,
    is_transformation_scheme,
)
from idaes.core.util.constants import Constants

"""
A simple pipeline unit model. We include a bulk momentum balance and leverage
the 1D control volume for a phase material balance. For now, we assume the
pipeline is isothermal.

Data sources:
    [1] Stochastic Optimal Control Model for Natural Gas Network
        Operations. V. Zavala, 2014, Comp. Chem. Eng.

"""

EXPLICIT_DISCRETIZATION_SCHEMES = {
    "FORWARD",
}
IMPLICIT_DISCRETIZATION_SCHEMES = {
    "BACKWARD",
    "LAGRANGE-RADAU",
}
UNSUPPORTED_DISCRETIZATION_SCHEMES = {
    "CENTRAL",
    "LAGRANGE-LEGENDRE",
}


@declare_process_block_class("GasPipeline")
class GasPipelineData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(default=None, domain=is_physical_parameter_block),
    )
    CONFIG.declare(
        "transformation_method",
        ConfigValue(
            default="dae.finite_difference",
            domain=is_transformation_method,
            description="DAE transformation method",
            doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory.""",
        ),
    )
    CONFIG.declare(
        "transformation_scheme",
        ConfigValue(
            default="FORWARD",
            domain=is_transformation_scheme,
            description="DAE transformation scheme",
            doc="""Scheme to use when transforming domain. See Pyomo
documentation for supported schemes.""",
        ),
    )
    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=1,
            domain=int,
            description="Number of finite elements",
            doc="""Number of finite elements to use in transformation (equivalent
to Pyomo nfe argument).""",
        ),
    )
    CONFIG.declare(
        "collocation_points",
        ConfigValue(
            default=None,
            domain=int,
            description="Number of collocation points",
            doc="""Number of collocation points to use (equivalent to Pyomo ncp
argument).""",
        ),
    )

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
                '%s can only be constructed with a single phase, "Vap".'
                "Got phase %s." % (self.__class__, self.phase)
            )

        #
        # Make sure property package has the properties that this unit
        # model requires.
        #
        property_dict = property_package.get_metadata().properties
        if "pressure" not in property_dict:
            raise ValueError(
                "Property package supplied to pipeline must have a property "
                "for 'pressure', which was not found in %s." % type(property_package)
            )

        #
        # Set up control volume:
        #
        # Anything not valid for CV config will have to be removed here.
        cv_config = config()
        self.control_volume = ControlVolume1DBlock(**cv_config)
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

        #
        # Deactivate momentum balance at point where discretization equation
        # doesn't exist.
        #
        t0 = self.flowsheet().time.first()
        x0 = self.control_volume.length_domain.first()
        xf = self.control_volume.length_domain.last()
        dynamic = self.flowsheet().config.dynamic
        if dynamic:
            # In a dynamic model with an implicit time discretization,
            # at t0 and the point in space where the spatial discretization
            # equation is not defined, the momentum balance contains two
            # variables, the time and space derivatives, that are not present
            # in any other equations. For the model to have a perfect matching,
            # this equation must be deactivated.
            scheme = self.config.transformation_scheme
            momentum_bal = self.control_volume.momentum_balance
            if scheme in IMPLICIT_DISCRETIZATION_SCHEMES:
                momentum_bal[t0, x0].deactivate()
            elif scheme in EXPLICIT_DISCRETIZATION_SCHEMES:
                momentum_bal[t0, xf].deactivate()
            elif scheme in UNSUPPORTED_DISCRETIZATION_SCHEMES:
                raise ValueError(
                    "Discretization scheme %s is not supported for "
                    "dynamic pipelines" % scheme
                )

        # NOTE: Later we may try to support switches in direction of flow.
        # At this time these may need to be renamed...
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

        self.add_isothermal_constraint()

    def add_friction_factor(self, rugosity=None):
        if rugosity is None:
            rugosity = 0.025 * pyunits.mm
        rug_value = pyo_value(rugosity)
        rug_units = pyunits.get_units(rugosity)
        self.rugosity = Param(
            initialize=rug_value,
            units=rug_units,
            # doc=TODO,
            mutable=True,
        )

        diam_mm = pyunits.convert(self.diameter, pyunits.mm)
        friction_factor_expr = (2 * log10(3.7 * diam_mm / self.rugosity)) ** -2

        self.friction_factor = Var(
            initialize=pyo_value(friction_factor_expr),
            units=pyunits.dimensionless,
            doc="This is fraction factor lambda_l in V. Zavala's paper",
        )

        self.friction_factor_eqn = Constraint(
            expr=(self.friction_factor == friction_factor_expr)
        )

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
        kghr = pyunits.kg / pyunits.hr
        cv = self.control_volume
        cv.flow_mass = Var(
            time, space, units=kghr, doc="Mass flow rate through the pipeline"
        )
        state = cv.properties

        def flow_mass_linking_rule(b, t, x):
            return b.flow_mass[t, x] - pyunits.convert(state[t, x].flow_mass, kghr) == 0

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
        nu = state.speed_of_sound
        friction_term = (
            8
            * lam
            * nu**2
            / Constants.pi**2
            / diameter**5
            * flow
            * abs(flow)
            / pressure
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
        """
        Bulk momentum balance. This is Equation 2.5 in [1]
        """
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
            flux_expr = pyunits.convert(cv.pressure_dx[t, x] / cv.length, kgm2hr2)
            friction_expr = pyunits.convert(self.get_friction_term(t, x), kgm2hr2)
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

        self.state_isothermal_eqn = Constraint(time, length, rule=isothermal_rule)
