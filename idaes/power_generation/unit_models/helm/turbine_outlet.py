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
Steam turbine outlet stage model.  This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136.
"""
__Author__ = "John Eslick"

from pyomo.common.config import In
from pyomo.environ import Var, sqrt, SolverFactory, value, Param, units as pyunits
from idaes.power_generation.unit_models.helm.turbine import HelmIsentropicTurbineData
from idaes.core import declare_process_block_class
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class(
    "HelmTurbineOutletStage",
    doc="Outlet stage steam turbine model",
)
class HelmTurbineOutletStageData(HelmIsentropicTurbineData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = HelmIsentropicTurbineData.CONFIG()

    def build(self):
        super().build()

        self.flow_coeff = Var(
            initialize=0.0333, doc="Turbine flow coefficient [kg*C^0.5/s/Pa]",
            units=pyunits.kg*pyunits.K**0.5/pyunits.s/pyunits.Pa
        )
        self.eff_dry = Var(initialize=0.87, doc="Turbine dry isentropic efficiency")
        self.design_exhaust_flow_vol = Var(
            initialize=6000.0, doc="Design exit volumetirc flowrate [m^3/s]",
            units=pyunits.m**3/pyunits.s
        )
        self.efficiency_mech = Var(initialize=1.0, doc="Turbine mechanical efficiency")
        self.efficiency_isentropic.unfix()
        self.eff_dry.fix()
        self.design_exhaust_flow_vol.fix()
        self.flow_coeff.fix()
        self.efficiency_mech.fix()

        @self.Expression(self.flowsheet().config.time, doc="Eff. fact. correlation")
        def tel(b, t):
            f = b.control_volume.properties_out[t].flow_vol / b.design_exhaust_flow_vol
            return 1e6 * (
                -0.0035 * f ** 5
                + 0.022 * f ** 4
                - 0.0542 * f ** 3
                + 0.0638 * f ** 2
                - 0.0328 * f
                + 0.0064
            )*pyunits.J/pyunits.mol

        @self.Constraint(self.flowsheet().config.time, doc="Stodola eq. choked flow")
        def stodola_equation(b, t):
            flow = b.control_volume.properties_in[t].flow_mol
            mw = b.control_volume.properties_in[t].mw
            Tin = b.control_volume.properties_in[t].temperature
            Pin = b.control_volume.properties_in[t].pressure
            Pr = b.ratioP[t]
            cf = b.flow_coeff

            return flow ** 2 * mw ** 2 * (Tin) == (
                cf ** 2 * Pin ** 2 * (1 - Pr ** 2))

        @self.Constraint(self.flowsheet().config.time, doc="Efficiency correlation")
        def efficiency_correlation(b, t):
            x = b.control_volume.properties_out[t].vapor_frac
            eff = b.efficiency_isentropic[t]
            dh_isen = b.delta_enth_isentropic[t]
            tel = b.tel[t]
            return eff == b.eff_dry * x * (1 - 0.65 * (1 - x)) * (1 + tel / dh_isen)

        @self.Expression(self.flowsheet().config.time, doc="Thermodynamic power [J/s]")
        def power_thermo(b, t):
            return b.control_volume.work[t]

        @self.Expression(self.flowsheet().config.time, doc="Shaft power [J/s]")
        def power_shaft(b, t):
            return b.power_thermo[t] * b.efficiency_mech

    def initialize(
        self,
        state_args={},
        outlvl=idaeslog.NOTSET,
        solver="ipopt",
        optarg={"tol": 1e-6, "max_iter": 30},
        calculate_cf=True,
    ):
        """
        Initialize the outlet turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl : sets output level of initialization routine
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        for t in self.flowsheet().config.time:
            if self.outlet.pressure[t].fixed:
                self.ratioP[t] = value(
                    self.outlet.pressure[t]/self.inlet.pressure[t])
                self.deltaP[t] = value(
                    self.outlet.pressure[t] - self.inlet.pressure[t])


        # Deactivate special constraints
        self.stodola_equation.deactivate()
        self.efficiency_correlation.deactivate()
        self.efficiency_isentropic.fix()
        self.deltaP.unfix()
        self.ratioP.unfix()
        self.inlet.fix()
        self.outlet.unfix()

        super().initialize(outlvl=outlvl, solver=solver, optarg=optarg)

        for t in self.flowsheet().config.time:
            mw = self.control_volume.properties_in[t].mw
            Tin = self.control_volume.properties_in[t].temperature
            Pin = self.control_volume.properties_in[t].pressure
            Pr = self.ratioP[t]
            if not calculate_cf:
                cf = self.flow_coeff
                self.inlet.flow_mol[t].fix(
                    value(cf * Pin * sqrt(1 - Pr ** 2) / mw / sqrt(Tin))
                )

        super().initialize(outlvl=outlvl, solver=solver, optarg=optarg)
        self.control_volume.properties_out[:].pressure.fix()

        # Free eff_isen and activate sepcial constarints
        self.efficiency_isentropic.unfix()
        self.outlet.pressure.fix()
        if calculate_cf:
            self.flow_coeff.unfix()
            self.inlet.flow_mol.unfix()
            self.inlet.flow_mol[0].fix()
            flow = self.control_volume.properties_in[0].flow_mol
            mw = self.control_volume.properties_in[0].mw
            Tin = self.control_volume.properties_in[0].temperature
            Pin = self.control_volume.properties_in[0].pressure
            Pr = self.ratioP[0]
            self.flow_coeff.value = value(
                flow * mw * sqrt(Tin/(1 - Pr ** 2))/Pin)

        else:
            self.inlet.flow_mol.unfix()

        self.stodola_equation.activate()
        self.efficiency_correlation.activate()
        slvr = SolverFactory(solver)
        slvr.options = optarg
        self.display()
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = slvr.solve(self, tee=slc.tee)
        init_log.info(
            "Initialization Complete (Outlet Stage): {}".format(idaeslog.condition(res))
        )

        # reload original spec
        if calculate_cf:
            cf = value(self.flow_coeff)
        from_json(self, sd=istate, wts=sp)
        if calculate_cf:
            # cf was probably fixed, so will have to set the value agian here
            # if you ask for it to be calculated.
            self.flow_coeff = cf

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        for t, c in self.stodola_equation.items():
            s = iscale.get_scaling_factor(
                self.control_volume.properties_in[t].flow_mol,
                default=1,
                warning=True)**2
            iscale.constraint_scaling_transform(c, s)
