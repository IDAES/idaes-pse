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
Simple model to represent resistance originating from the contact between the
flow mesh and electrode or flow mesh and interconnect. An equality constraint is
created to add a resistance heating term from the heat flux coming from one side
to that coming out the other.

Boundary variables:
    - ``temperature_deviation[t, iz]``
    - ``heat_flux_x0[t, iz]``
    - ``heat_flux_x1[t, iz]``

Instances of ``Var`` that must be fixed:
    - ``log_preexponential_factor``: Natural logarithm of area squared resistance preexponential factor in ohm * m**2
    - ``thermal_exponent_dividend``: Parameter divided by temperature in area squared resistance equation, in K.
      Would be something like (reduced) activation energy, but it can be both negative and positive.
    - ``contact_fraction``: Fraction of area at which both surfaces touch. If unknown, can fix at one.
"""
__author__ = "Douglas Allan"

from pyomo.common.config import ConfigBlock, ConfigValue, In
import pyomo.environ as pyo


from idaes.core import declare_process_block_class, UnitModelBlockData
import idaes.models_extra.power_generation.unit_models.soc_submodels.common as common
import idaes.core.util.scaling as iscale
from idaes.core.solvers import get_solver

import idaes.logger as idaeslog


@declare_process_block_class("SocContactResistor")
class SocContactResistorData(UnitModelBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag",
            doc="No capacities or holdups, so no internal dynamics",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(domain=In([False]), default=False),
    )
    common._submodel_boilerplate_config(CONFIG)
    common._thermal_boundary_conditions_config(CONFIG, thin=True)

    def build(self):
        super().build()

        # Set up some sets for the space and time indexing
        tset = self.flowsheet().config.time
        # Set up node and face sets and get integer indices for them
        izfaces, iznodes = common._face_initializer(
            self, self.config.control_volume_zfaces, "z"
        )
        common._submodel_boilerplate_create_if_none(self)
        common._create_thermal_boundary_conditions_if_none(self, thin=True)

        # Preexponential factor needs to be given in units of ohm*m**2
        self.log_preexponential_factor = pyo.Var(
            initialize=-50, units=pyo.units.dimensionless
        )
        self.thermal_exponent_dividend = pyo.Var(initialize=0, units=pyo.units.K)
        self.contact_fraction = pyo.Var(
            initialize=1, units=pyo.units.dimensionless, bounds=(0, 1)
        )

        @self.Expression(tset, iznodes)
        def contact_resistance(b, t, iz):
            return (
                1
                * pyo.units.ohm
                * pyo.units.m**2
                / b.contact_fraction
                * pyo.exp(
                    b.log_preexponential_factor
                    + b.thermal_exponent_dividend / (b.temperature[t, iz])
                )
            )

        @self.Expression(tset, iznodes)
        def voltage_drop_total(b, t, iz):
            return b.contact_resistance[t, iz] * b.current_density[t, iz]

        @self.Expression(tset, iznodes)
        def joule_heating_flux(b, t, iz):
            return b.voltage_drop_total[t, iz] * b.current_density[t, iz]

        @self.Constraint(tset, iznodes)
        def heat_flux_x_eqn(b, t, iz):
            return (
                b.heat_flux_x1[t, iz]
                == b.heat_flux_x0[t, iz] + b.joule_heating_flux[t, iz]
            )

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        fix_heat_flux_x0=True,
        temperature_guess=None,
    ):
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        if temperature_guess is not None:
            for t in self.flowsheet().time:
                for iz in self.iznodes:
                    common._set_if_unfixed(self.temperature_z[t, iz], temperature_guess)

        self.temperature_deviation_x.fix()
        if fix_heat_flux_x0:
            self.heat_flux_x0.fix()
        else:
            self.heat_flux_x1.fix()

        solver_obj = get_solver(solver, optarg)
        common._init_solve_block(self, solver_obj, solve_log)

        self.temperature_deviation_x.unfix()
        if fix_heat_flux_x0:
            self.heat_flux_x0.unfix()
        else:
            self.heat_flux_x1.unfix()

    def calculate_scaling_factors(self):
        pass

    def recursive_scaling(self):
        gsf = iscale.get_scaling_factor
        for i, c in self.heat_flux_x_eqn.items():
            if self.heat_flux_x0[i].is_reference():
                sq0 = gsf(self.heat_flux_x0[i].referent)
            else:
                sq0 = gsf(self.heat_flux_x0[i], warning=True)
            if self.heat_flux_x1[i].is_reference():
                sq1 = gsf(self.heat_flux_x1[i].referent)
            else:
                sq1 = gsf(self.heat_flux_x1[i], warning=True)
            sq = min(sq0, sq1)
            iscale.constraint_scaling_transform(c, sq, overwrite=False)
