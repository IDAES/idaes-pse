##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Steam turbine outlet stage model.  This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136, November
"""
from __future__ import division

__Author__ = "John Eslick"

import logging
_log = logging.getLogger(__name__)

from pyomo.common.config import In
from pyomo.environ import Var, Expression, Constraint, sqrt, SolverFactory, value
from pyomo.opt import TerminationCondition

from idaes.core import declare_process_block_class
from idaes.unit_models.pressure_changer import PressureChangerData
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.ui.report import degrees_of_freedom


@declare_process_block_class("TurbineOutletStage",
    doc="Outlet stage steam turbine model")
class TurbineOutletStageData(PressureChangerData):
    # Same setings as the default pressure changer, but force to expander with
    # isentroic efficiency
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = False
    CONFIG.get('compressor')._default = False
    CONFIG.get('compressor')._domain = In([False])
    CONFIG.thermodynamic_assumption = 'isentropic'
    CONFIG.get('thermodynamic_assumption')._default = 'isentropic'
    CONFIG.get('thermodynamic_assumption')._domain = In(['isentropic'])
    def build(self):
        super(TurbineOutletStageData, self).build()

        self.flow_coeff = Var(initialize=0.0333,
            doc="Turbine flow coefficient [kg*C^0.5/s/Pa]")
        self.delta_enth_isentropic = Var(self.time_ref, initialize=-100,
            doc="Specific enthalpy change of isentropic process [J/mol]")
        self.eff_dry = Var(initialize=0.87,
            doc="Turbine dry isentropic efficiency")
        self.design_exhaust_flow_vol = Var(initialize=6000.0,
            doc="Design exit volumetirc flowrate [m^3/s]")
        self.efficiency_mech = Var(initialize=0.98,
            doc="Turbine mechanical efficiency")
        self.eff_dry.fix()
        self.design_exhaust_flow_vol.fix()
        self.flow_coeff.fix()
        self.efficiency_mech.fix()
        self.ratioP[:] = 1 # make sure these have a number value
        self.deltaP[:] = 0 #   to avoid an error later in initialize

        @self.Expression(self.time_ref, doc="Efficiency factor correlation")
        def tel(b, t):
            f = b.control_volume.properties_out[t].flow_vol/b.design_exhaust_flow_vol
            return 1e6*(-0.0035*f**5 + 0.022*f**4 - 0.0542*f**3 + 0.0638*f**2 -
                        0.0328*f + 0.0064)

        @self.Constraint(self.time_ref, doc="Equation: Stodola, for choked flow")
        def stodola_equation(b, t):
            flow = b.control_volume.properties_in[t].flow_mol
            mw = b.control_volume.properties_in[t].mw
            Tin = b.control_volume.properties_in[t].temperature
            Pin = b.control_volume.properties_in[t].pressure
            Pr = b.ratioP[t]
            cf = b.flow_coeff
            return flow*mw*sqrt(Tin - 273.15) == cf*Pin*sqrt(1 - Pr**2)

        @self.Constraint( self.time_ref,
            doc="Equation: isentropic specific enthalpy change")
        def isentropic_enthalpy(b, t):
            flow = b.control_volume.properties_in[t].flow_mol
            dh_isen = b.delta_enth_isentropic[t]
            work_isen = b.work_isentropic[t]
            return work_isen == dh_isen*flow

        @self.Constraint(self.time_ref, doc="Equation: Efficiency correlation")
        def efficiency_correlation(b, t):
            x = b.control_volume.properties_out[t].vapor_frac
            eff = b.efficiency_isentropic[t]
            dh_isen = b.delta_enth_isentropic[t]
            tel = b.tel[t]
            return eff == b.eff_dry*x*(1 - 0.65*(1 - x))*(1 + tel/dh_isen)

        @self.Expression(self.time_ref, doc="Thermodynamic power [J/s]")
        def power_thermo(b, t):
            return b.control_volume.work[t]

        @self.Expression(self.time_ref, doc="Shaft power [J/s]")
        def power_shaft(b, t):
            return b.power_thermo[t]*b.efficiency_mech

    def initialize(self, state_args={}, outlvl=0, solver='ipopt',
        optarg={'tol': 1e-6}):
        """
        Initialize the outlet turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl (int): Amount of output (0 to 3) 0 is lowest
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        stee = True if outlvl >= 3 else False
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)
        # Deactivate special constraints
        self.stodola_equation.deactivate()
        self.isentropic_enthalpy.deactivate()
        self.efficiency_correlation.deactivate()
        self.deltaP[:].unfix()
        self.ratioP[:].unfix()
        # Fix turbine parameters + eff_isen
        self.eff_dry.fix()
        self.design_exhaust_flow_vol.fix()
        self.flow_coeff.fix()

        # fix inlet and free outlet
        for t in self.time_ref:
            for k, v in self.inlet[t].vars.items():
                v.fix()
            for k, v in self.outlet[t].vars.items():
                v.unfix()
            # If there isn't a good guess for efficeny or outlet pressure
            # provide something reasonable.
            eff = self.efficiency_isentropic[t]
            eff.fix(eff.value if value(eff) > 0.3 and value(eff) < 1.0 else 0.8)
            # for outlet pressure try outlet pressure, pressure ratio, delta P,
            # then if none of those look reasonable use a pressure ratio of 0.8
            # to calculate outlet pressure
            Pout = self.outlet[t].pressure
            Pin = self.inlet[t].pressure
            prdp = value((self.deltaP[t] - Pin)/Pin)
            if value(Pout/Pin) > 0.9 or value(Pout/Pin) < 0.01:
                if value(self.ratioP[t]) < 0.9 and value(self.ratioP[t]) > 0.01:
                    Pout.fix(value(Pin*self.ratioP))
                elif prdp < 0.9 and prdp > 0.01:
                    Pout.fix(value(prdp*Pin))
                else:
                    Pout.fix(value(Pin*0.3))
            else:
                Pout.fix()
        self.deltaP[:] = value(Pout - Pin)
        self.ratioP[:] = value(Pout/Pin)

        # Make sure the initialization problem has no degrees of freedom
        # This shouldn't happen here unless there is a bug in this
        dof = degrees_of_freedom(self)
        try:
            assert(dof == 0)
        except:
            _log.exception("degrees_of_freedom = {}".format(dof))
            raise

        # one bad thing about reusing this is that the log messages aren't
        # really compatible with being nested inside another initialization
        super(TurbineOutletStageData, self).initialize(state_args=state_args,
            outlvl=outlvl, solver=solver, optarg=optarg)

        # Free eff_isen and activate sepcial constarints
        self.efficiency_isentropic.unfix()
        self.outlet[:].pressure.unfix()
        self.stodola_equation.activate()
        self.isentropic_enthalpy.activate()
        self.efficiency_correlation.activate()

        slvr = SolverFactory(solver)
        slvr.options = optarg
        res = slvr.solve(self, tee=stee)

        if outlvl > 0:
            if res.solver.termination_condition == TerminationCondition.optimal:
                _log.info("{} Initialization Complete.".format(self.name))
            else:
                _log.warning(
"""{} Initialization Failed. The most likely cause of initialization failure for
the Turbine inlet stages model is that the flow coefficent is not compatable
with flow rate guess.""".format(self.name))

        # reload original spec
        from_json(self, sd=istate, wts=sp)
