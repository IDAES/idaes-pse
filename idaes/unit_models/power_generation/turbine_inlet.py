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
Steam turbine inlet stage model.  This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136, November
"""
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

@declare_process_block_class("TurbineInletStage",
    doc="Inlet stage steam turbine model")
class TurbineInletStageData(PressureChangerData):
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
        super(TurbineInletStageData, self).build()

        self.flow_coeff = Var(self.time_ref, initialize=1.053/3600.0,
            doc="Turbine flow coefficient [kg*C^0.5/Pa/s]")
        self.delta_enth_isentropic = Var(self.time_ref, initialize=-1000,
            doc="Specific enthalpy change of isentropic process [J/mol/K]")
        self.blade_reaction = Var(initialize=0.9,
            doc="Blade reaction parameter")
        self.blade_velocity = Var(initialize=110.0,
            doc="Design blade velocity [m/s]")
        self.eff_nozzle = Var(initialize=0.95, bounds=(0.0, 1.0),
            doc="Nozzel efficiency (typically 0.90 to 0.95)")
        self.eff_nozzle.fix()
        self.blade_reaction.fix()
        self.flow_coeff.fix()
        self.blade_velocity.fix()
        self.control_volume.deltaP[:] = -1000

        @self.Expression(self.time_ref,
            doc="Entering steam velocity calculation [m/s]")
        def steam_entering_velocity(b, t):
            # 1.414 = 44.72/sqrt(1000) for SI if comparing to Liese (2014)
            # b.delta_enth_isentropic[t] = -(hin - hiesn), the mw converts
            # enthalpy to a mass basis
            return 1.414*sqrt(-(1-b.blade_reaction)*b.delta_enth_isentropic[t]/
                    b.control_volume.properties_in[t].mw*self.eff_nozzle)

        @self.Constraint(self.time_ref, doc="Equation: Turbine inlet flow")
        def inlet_flow_constraint(b, t):
            # Some local vars to make the equation more readable
            g = b.control_volume.properties_in[t].heat_capacity_ratio
            mw = b.control_volume.properties_in[t].mw
            flow = b.control_volume.properties_in[t].flow_mol
            Tin = b.control_volume.properties_in[t].temperature
            cf = b.flow_coeff[t]
            Pin = b.control_volume.properties_in[t].pressure
            Pratio = b.ratioP[t]
            return (flow*mw*sqrt(Tin - 273.15) ==
                cf*Pin*sqrt(g/(g - 1)*(Pratio**(2.0/g) - Pratio**((g + 1)/g))))

        @self.Constraint(self.time_ref, doc="Equation: Isentropic enthalpy change")
        def isentropic_enthalpy(b, t):
            return b.work_isentropic[t] == (b.delta_enth_isentropic[t]*
                b.control_volume.properties_in[t].flow_mol)

        @self.Constraint(self.time_ref, doc="Equation: Efficiency")
        def efficiency_correlation(b, t):
            Vr = b.blade_velocity/b.steam_entering_velocity[t]
            eff = b.efficiency_isentropic[t]
            R = b.blade_reaction
            return eff == 2*Vr*((sqrt(1 - R) - Vr) +
                                 sqrt((sqrt(1 - R) - Vr)**2 + R))

    def initialize(self, state_args={}, outlvl=0, solver='ipopt',
        optarg={'tol': 1e-6}):
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
        stee = True if outlvl >= 3 else False
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

        # Fix turbine parameters + eff_isen
        self.eff_nozzle.fix()
        self.blade_reaction()
        self.flow_coeff.fix()
        self.blade_velocity.fix()

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
            if value(Pout/Pin) > 0.98 or value(Pout/Pin) < 0.3:
                if value(self.ratioP[t]) < 0.98 and value(self.ratioP[t]) > 0.3:
                    Pout.fix(value(Pin*self.ratioP))
                elif prdp < 0.98 and prdp > 0.3:
                    Pout.fix(value(prdp*Pin))
                else:
                    Pout.fix(value(Pin*0.8))
            else:
                Pout.fix()

        # Make sure the initialization problem has no degrees of freedom
        # This shouldn't happen here unless there is a bug in this
        dof = degrees_of_freedom(self)
        try:
            assert(dof == 0)
        except:
            _log.exception("degrees_of_freedom = {}".format(dof))
            raise

        # one bad thing about reusing this is that the log messages aren't
        # really compatable with being nested inside another initialization
        super(TurbineInletStageData, self).initialize(state_args=state_args,
            outlvl=outlvl, solver=solver, optarg=optarg)

        # Free eff_isen and activate sepcial constarints
        self.efficiency_isentropic.unfix()
        self.outlet[:].pressure.unfix()
        self.inlet_flow_constraint.activate()
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
