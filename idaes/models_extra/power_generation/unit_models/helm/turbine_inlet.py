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
Steam turbine inlet stage model.  This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136.
"""
__Author__ = "John Eslick"

from pyomo.environ import Var, sqrt, value, units as pyunits
from idaes.core import declare_process_block_class
from idaes.models_extra.power_generation.unit_models.helm.turbine import (
    HelmIsentropicTurbineData,
)
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.solvers import get_solver
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale


_log = idaeslog.getLogger(__name__)


@declare_process_block_class(
    "HelmTurbineInletStage",
    doc="Inlet stage steam turbine model",
)
class HelmTurbineInletStageData(HelmIsentropicTurbineData):
    CONFIG = HelmIsentropicTurbineData.CONFIG()

    def build(self):
        super().build()

        self.flow_coeff = Var(
            self.flowsheet().time,
            initialize=1.053 / 3600.0,
            doc="Turbine flow coefficient [kg*C^0.5/Pa/s]",
            units=pyunits.kg * pyunits.K**0.5 / pyunits.Pa / pyunits.s,
        )
        self.blade_reaction = Var(initialize=0.9, doc="Blade reaction parameter")
        self.blade_velocity = Var(
            initialize=110.0,
            doc="Design blade velocity [m/s]",
            units=pyunits.m / pyunits.s,
        )
        self.eff_nozzle = Var(
            initialize=0.95,
            bounds=(0.0, 1.0),
            doc="Nozzel efficiency (typically 0.90 to 0.95)",
        )
        self.efficiency_mech = Var(initialize=1.0, doc="Turbine mechanical efficiency")

        self.eff_nozzle.fix()
        self.blade_reaction.fix()
        self.flow_coeff.fix()
        self.blade_velocity.fix()
        self.efficiency_mech.fix()
        self.efficiency_isentropic.unfix()
        self.ratioP[:] = 0.9  # make sure these have a number value
        self.deltaP[:] = 0  #   to avoid an error later in initialize

        @self.Expression(
            self.flowsheet().time,
            doc="Entering steam velocity calculation [m/s]",
        )
        def steam_entering_velocity(b, t):
            # 1.414 = 44.72/sqrt(1000) for SI if comparing to Liese (2014),
            # b.delta_enth_isentropic[t] = -(hin - hiesn), the mw converts
            # enthalpy to a mass basis
            return 1.414 * sqrt(
                (b.blade_reaction - 1)
                * b.delta_enth_isentropic[t]
                * self.eff_nozzle
                / b.control_volume.properties_in[t].mw
            )

        @self.Expression(self.flowsheet().time, doc="Efficiency expression")
        def efficiency_isentropic_expr(b, t):
            Vr = b.blade_velocity / b.steam_entering_velocity[t]
            R = b.blade_reaction
            return 2 * Vr * ((sqrt(1 - R) - Vr) + sqrt((sqrt(1 - R) - Vr) ** 2 + R))

        @self.Constraint(self.flowsheet().time, doc="Equation: Turbine inlet flow")
        def inlet_flow_constraint(b, t):
            # Some local vars to make the equation more readable
            g = b.control_volume.properties_in[t].heat_capacity_ratio
            mw = b.control_volume.properties_in[t].mw
            flow = b.control_volume.properties_in[t].flow_mol
            Tin = b.control_volume.properties_in[t].temperature
            cf = b.flow_coeff[t]
            Pin = b.control_volume.properties_in[t].pressure
            Pratio = b.ratioP[t]

            return flow**2 * mw**2 * Tin == (
                cf**2
                * Pin**2
                * g
                / (g - 1)
                * (Pratio ** (2.0 / g) - Pratio ** ((g + 1) / g))
            )

        @self.Constraint(self.flowsheet().time, doc="Equation: Efficiency")
        def efficiency_correlation(b, t):
            return b.efficiency_isentropic[t] == b.efficiency_isentropic_expr[t]

        @self.Expression(self.flowsheet().time, doc="Thermodynamic power [J/s]")
        def power_thermo(b, t):
            return b.control_volume.work[t]

        @self.Expression(self.flowsheet().time, doc="Shaft power [J/s]")
        def power_shaft(b, t):
            return b.power_thermo[t] * b.efficiency_mech

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
        calculate_cf=False,
    ):
        """
        Initialize the inlet turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves. This initializtion uses a
        flow value guess, so some reasonable flow guess should be sepecified prior
        to initializtion.

        Args:
            outlvl (int): Amount of output (0 to 3) 0 is lowest
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
            calculate_cf (bool): If True, use the flow and pressure ratio to
                calculate the flow coefficient.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # sp is what to save to make sure state after init is same as the start
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # Setup for initializtion step 1
        self.inlet_flow_constraint.deactivate()
        self.efficiency_correlation.deactivate()
        self.eff_nozzle.fix()
        self.blade_reaction.fix()
        self.flow_coeff.fix()
        self.blade_velocity.fix()
        self.inlet.fix()
        self.outlet.unfix()

        for t in self.flowsheet().time:
            self.efficiency_isentropic[t] = 0.9
        super().initialize_build(outlvl=outlvl, solver=solver, optarg=optarg)

        # Free eff_isen and activate sepcial constarints
        self.inlet_flow_constraint.activate()
        self.efficiency_correlation.activate()

        if calculate_cf:
            self.ratioP.fix()
            self.flow_coeff.unfix()

            for t in self.flowsheet().time:
                g = self.control_volume.properties_in[t].heat_capacity_ratio
                mw = self.control_volume.properties_in[t].mw
                flow = self.control_volume.properties_in[t].flow_mol
                Tin = self.control_volume.properties_in[t].temperature
                Pin = self.control_volume.properties_in[t].pressure
                Pratio = self.ratioP[t]
                self.flow_coeff[t].value = value(
                    flow
                    * mw
                    * sqrt(
                        Tin
                        / (
                            g
                            / (g - 1)
                            * (Pratio ** (2.0 / g) - Pratio ** ((g + 1) / g))
                        )
                    )
                    / Pin
                )

        # Create solver
        slvr = get_solver(solver, optarg)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = slvr.solve(self, tee=slc.tee)
        init_log.info("Initialization Complete: {}".format(idaeslog.condition(res)))
        # reload original spec
        if calculate_cf:
            cf = {}
            for t in self.flowsheet().time:
                cf[t] = value(self.flow_coeff[t])

        from_json(self, sd=istate, wts=sp)
        if calculate_cf:
            # cf was probably fixed, so will have to set the value agian here
            # if you ask for it to be calculated.
            for t in self.flowsheet().time:
                self.flow_coeff[t] = cf[t]

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t, c in self.inlet_flow_constraint.items():
            s = (
                iscale.get_scaling_factor(self.control_volume.properties_in[t].flow_mol)
                ** 2
            )
            iscale.constraint_scaling_transform(c, s, overwrite=False)
