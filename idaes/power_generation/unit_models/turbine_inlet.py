##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Steam turbine inlet stage model.  This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136.
"""
__Author__ = "John Eslick"

from pyomo.common.config import In
from pyomo.environ import Var, Expression, Constraint, sqrt, SolverFactory, value, Param
from pyomo.opt import TerminationCondition

from idaes.core import declare_process_block_class
from idaes.generic_models.unit_models.pressure_changer import (
    PressureChangerData,
    ThermodynamicAssumption,
)
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("TurbineInletStage", doc="Inlet stage steam turbine model")
class TurbineInletStageData(PressureChangerData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = False
    CONFIG.get("compressor")._default = False
    CONFIG.get("compressor")._domain = In([False])
    CONFIG.thermodynamic_assumption = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._default = ThermodynamicAssumption.isentropic
    CONFIG.get("thermodynamic_assumption")._domain = In(
        [ThermodynamicAssumption.isentropic]
    )

    def build(self):
        super(TurbineInletStageData, self).build()

        self.flow_coeff = Var(
            self.flowsheet().config.time,
            initialize=1.053 / 3600.0,
            doc="Turbine flow coefficient [kg*C^0.5/Pa/s]",
        )
        self.delta_enth_isentropic = Var(
            self.flowsheet().config.time,
            initialize=-1000,
            doc="Specific enthalpy change of isentropic process [J/mol]",
        )
        self.blade_reaction = Var(initialize=0.9, doc="Blade reaction parameter")
        self.blade_velocity = Var(initialize=110.0, doc="Design blade velocity [m/s]")
        self.eff_nozzle = Var(
            initialize=0.95,
            bounds=(0.0, 1.0),
            doc="Nozzel efficiency (typically 0.90 to 0.95)",
        )
        self.efficiency_mech = Var(initialize=0.98, doc="Turbine mechanical efficiency")
        self.flow_scale = Param(
            mutable=True,
            default=1e3,
            doc="Scaling factor for pressure flow relation should be approximatly"
            " the same order of magnitude as the expected flow.",
        )
        self.eff_nozzle.fix()
        self.blade_reaction.fix()
        self.flow_coeff.fix()
        self.blade_velocity.fix()
        self.efficiency_mech.fix()
        self.ratioP[:] = 1  # make sure these have a number value
        self.deltaP[:] = 0  #   to avoid an error later in initialize

        @self.Expression(
            self.flowsheet().config.time,
            doc="Entering steam velocity calculation [m/s]",
        )
        def steam_entering_velocity(b, t):
            # 1.414 = 44.72/sqrt(1000) for SI if comparing to Liese (2014)
            # b.delta_enth_isentropic[t] = -(hin - hiesn), the mw converts
            # enthalpy to a mass basis
            return 1.414 * sqrt(
                -(1 - b.blade_reaction)
                * b.delta_enth_isentropic[t]
                / b.control_volume.properties_in[t].mw
                * self.eff_nozzle
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Equation: Turbine inlet flow"
        )
        def inlet_flow_constraint(b, t):
            # Some local vars to make the equation more readable
            g = b.control_volume.properties_in[t].heat_capacity_ratio
            mw = b.control_volume.properties_in[t].mw
            flow = b.control_volume.properties_in[t].flow_mol
            Tin = b.control_volume.properties_in[t].temperature
            cf = b.flow_coeff[t]
            Pin = b.control_volume.properties_in[t].pressure
            Pratio = b.ratioP[t]
            return (1 / b.flow_scale ** 2) * flow ** 2 * mw ** 2 * (Tin - 273.15) == (
                1 / b.flow_scale ** 2
            ) * cf ** 2 * Pin ** 2 * (
                g / (g - 1) * (Pratio ** (2.0 / g) - Pratio ** ((g + 1) / g))
            )

        @self.Constraint(
            self.flowsheet().config.time, doc="Equation: Isentropic enthalpy change"
        )
        def isentropic_enthalpy(b, t):
            return b.work_isentropic[t] == (
                b.delta_enth_isentropic[t] * b.control_volume.properties_in[t].flow_mol
            )

        @self.Constraint(self.flowsheet().config.time, doc="Equation: Efficiency")
        def efficiency_correlation(b, t):
            Vr = b.blade_velocity / b.steam_entering_velocity[t]
            eff = b.efficiency_isentropic[t]
            R = b.blade_reaction
            return eff == 2 * Vr * (
                (sqrt(1 - R) - Vr) + sqrt((sqrt(1 - R) - Vr) ** 2 + R)
            )

        @self.Expression(self.flowsheet().config.time, doc="Thermodynamic power [J/s]")
        def power_thermo(b, t):
            return b.control_volume.work[t]

        @self.Expression(self.flowsheet().config.time, doc="Shaft power [J/s]")
        def power_shaft(b, t):
            return b.power_thermo[t] * b.efficiency_mech

    def _get_performance_contents(self, time_point=0):
        pc = super()._get_performance_contents(time_point=time_point)
        pc["vars"]["Mechanical Efficiency"] = self.efficiency_mech
        pc["vars"]["Flow Coefficient"] = self.flow_coeff[time_point]
        pc["vars"]["Isentropic Specific Enthalpy"] = self.delta_enth_isentropic[
            time_point
        ]
        pc["vars"]["Blade Reaction"] = self.blade_reaction
        pc["vars"]["Blade Velocity"] = self.blade_velocity
        pc["vars"]["Nozzel Efficiency"] = self.eff_nozzle

        pc["exprs"] = {}
        pc["exprs"]["Thermodynamic Power"] = self.power_thermo[time_point]
        pc["exprs"]["Shaft Power"] = self.power_shaft[time_point]
        pc["exprs"]["Inlet Steam Velocity"] = self.steam_entering_velocity[time_point]

        pc["params"] = {}
        pc["params"]["Flow Scaling"] = self.flow_scale

        return pc

    def initialize(
        self,
        state_args={},
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6, "max_iter": 30},
    ):
        """
        Initialize the inlet turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl (int): Amount of output (0 to 3) 0 is lowest
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)
        # Deactivate special constraints
        self.inlet_flow_constraint.deactivate()
        self.isentropic_enthalpy.deactivate()
        self.efficiency_correlation.deactivate()
        self.deltaP.unfix()
        self.ratioP.unfix()

        # Fix turbine parameters + eff_isen
        self.eff_nozzle.fix()
        self.blade_reaction.fix()
        self.flow_coeff.fix()
        self.blade_velocity.fix()

        # fix inlet and free outlet
        for t in self.flowsheet().config.time:
            for k, v in self.inlet.vars.items():
                v[t].fix()
            for k, v in self.outlet.vars.items():
                v[t].unfix()
            # If there isn't a good guess for efficeny or outlet pressure
            # provide something reasonable.
            eff = self.efficiency_isentropic[t]
            eff.fix(eff.value if value(eff) > 0.3 and value(eff) < 1.0 else 0.8)
            # for outlet pressure try outlet pressure, pressure ratio, delta P,
            # then if none of those look reasonable use a pressure ratio of 0.8
            # to calculate outlet pressure
            Pout = self.outlet.pressure[t]
            Pin = self.inlet.pressure[t]
            prdp = value((self.deltaP[t] - Pin) / Pin)
            if value(Pout / Pin) > 0.98 or value(Pout / Pin) < 0.3:
                if value(self.ratioP[t]) < 0.98 and value(self.ratioP[t]) > 0.3:
                    Pout.fix(value(Pin * self.ratioP))
                elif prdp < 0.98 and prdp > 0.3:
                    Pout.fix(value(prdp * Pin))
                else:
                    Pout.fix(value(Pin * 0.8))
            else:
                Pout.fix()
        self.deltaP[:] = value(Pout - Pin)
        self.ratioP[:] = value(Pout / Pin)

        for t in self.flowsheet().config.time:
            self.properties_isentropic[t].pressure.value = value(
                self.outlet.pressure[t]
            )
            self.properties_isentropic[t].flow_mol.value = value(self.inlet.flow_mol[t])
            self.properties_isentropic[t].enth_mol.value = value(
                self.inlet.enth_mol[t] * 0.95
            )
            self.outlet.flow_mol[t].value = value(self.inlet.flow_mol[t])
            self.outlet.enth_mol[t].value = value(self.inlet.enth_mol[t] * 0.95)

        # Make sure the initialization problem has no degrees of freedom
        # This shouldn't happen here unless there is a bug in this
        dof = degrees_of_freedom(self)
        try:
            assert dof == 0
        except:
            init_log.exception("degrees_of_freedom = {}".format(dof))
            raise

        # one bad thing about reusing this is that the log messages aren't
        # really compatible with being nested inside another initialization
        super().initialize(
            state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
        )

        # Free eff_isen and activate sepcial constarints
        self.efficiency_isentropic.unfix()
        self.outlet.pressure.unfix()
        self.inlet_flow_constraint.activate()
        self.isentropic_enthalpy.activate()
        self.efficiency_correlation.activate()

        slvr = SolverFactory(solver)
        slvr.options = optarg
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = slvr.solve(self, tee=slc.tee)

        init_log.info(
            "Initialization Complete: {}".format(idaeslog.condition(res))
        )

        # reload original spec
        from_json(self, sd=istate, wts=sp)
