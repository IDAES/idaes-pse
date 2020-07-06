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
Steam turbine stage model. This is a standard isentropic turbine. Under off-design
conditions the base efficiency and pressure ratio do not change much for the stages
between the inlet and outlet. This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136.
"""
__Author__ = "John Eslick"


from pyomo.common.config import In
from pyomo.environ import Var, Expression, SolverFactory, value
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


@declare_process_block_class("TurbineStage", doc="Basic steam turbine model")
class TurbineStageData(PressureChangerData):
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
        super(TurbineStageData, self).build()

        self.efficiency_mech = Var(initialize=0.98, doc="Turbine mechanical efficiency")
        self.efficiency_mech.fix()
        self.ratioP[:] = 0.8  # make sure these have a number value
        self.deltaP[:] = 0  #   to avoid an error later in initialize
        time_set = self.flowsheet().config.time

        self.shaft_speed = Var(time_set, doc="Shaft speed [1/s]", initialize=60.0)

        @self.Expression(time_set, doc="Specific speed [dimensionless]")
        def specific_speed(b, t):
            s = b.shaft_speed[t]  # 1/s
            v = b.control_volume.properties_out[t].flow_vol  # m3/s
            his_rate = b.work_isentropic[t]  # J/s
            m = b.control_volume.properties_out[t].flow_mass  # kg/s
            return s * v ** 0.5 * (his_rate / m) ** (-0.75)  # dimensionless

        @self.Expression(time_set, doc="Thermodynamic power [J/s]")
        def power_thermo(b, t):
            return b.control_volume.work[t]

        @self.Expression(self.flowsheet().config.time, doc="Shaft power [J/s]")
        def power_shaft(b, t):
            return b.power_thermo[t] * b.efficiency_mech

    def _get_performance_contents(self, time_point=0):
        pc = super()._get_performance_contents(time_point=time_point)
        pc["vars"]["Mechanical Efficiency"] = self.efficiency_mech

        return pc

    def initialize(
        self,
        state_args={},
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6, "max_iter": 30},
    ):
        """
        Initialize the turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl : sets output level of initialization routine
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # fix inlet and free outlet
        for t in self.flowsheet().config.time:
            for k, v in self.inlet.vars.items():
                v[t].fix()
            for k, v in self.outlet.vars.items():
                v[t].unfix()
            # If there isn't a good guess for efficiency or outlet pressure
            # provide something reasonable.
            eff = self.efficiency_isentropic[t]
            eff.fix(eff.value if value(eff) > 0.3 and value(eff) < 1.0 else 0.8)
            # for outlet pressure try outlet pressure, pressure ratio, delta P,
            # then if none of those look reasonable use a pressure ratio of 0.8
            # to calculate outlet pressure
            Pout = self.outlet.pressure[t]
            Pin = self.inlet.pressure[t]
            prdp = value((self.deltaP[t] - Pin) / Pin)
            if self.deltaP[t].fixed:
                Pout.value = value(Pin - Pout)
            if self.ratioP[t].fixed:
                Pout.value = value(self.ratioP[t] * Pin)
            if value(Pout / Pin) > 0.99 or value(Pout / Pin) < 0.1:
                if value(self.ratioP[t]) < 0.99 and value(self.ratioP[t]) > 0.1:
                    Pout.fix(value(Pin * self.ratioP[t]))
                elif prdp < 0.99 and prdp > 0.1:
                    Pout.fix(value(prdp * Pin))
                else:
                    Pout.fix(value(Pin * 0.8))
            else:
                Pout.fix()
            self.deltaP[t] = value(Pout - Pin)
            self.ratioP[t] = value(Pout / Pin)

        self.deltaP[:].unfix()
        self.ratioP[:].unfix()

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

        # one bad thing about reusing this is that the log messages aren't
        # really compatible with being nested inside another initialization
        super().initialize(
            state_args=state_args, outlvl=outlvl, solver=solver, optarg=optarg
        )

        # reload original spec
        from_json(self, sd=istate, wts=sp)
