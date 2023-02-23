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
import pyomo.environ as pyo
from pyomo.common.config import In
from idaes.core import declare_process_block_class
from idaes.models_extra.power_generation.unit_models.balance import BalanceBlockData
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.core.solvers import get_solver
import idaes.models.properties.helmholtz.helmholtz as hltz
from idaes.models.properties.helmholtz.helmholtz import (
    HelmholtzThermoExpressions as ThermoExpr,
)
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _assert_properties(pb):
    """Assert that the properies parameter block conforms to the requirements"""
    try:
        assert isinstance(pb, hltz.HelmholtzParameterBlockData)
        assert pb.config.phase_presentation in {
            hltz.PhaseType.MIX,
            hltz.PhaseType.L,
            hltz.PhaseType.G,
        }
        assert pb.config.state_vars == hltz.StateVars.PH
    except AssertionError:
        _log.error(
            "helm.HelmPump requires a Helmholtz EOS with "
            "a single or mixed phase and pressure-enthalpy state vars."
        )
        raise


@declare_process_block_class("HelmPump")
class HelmPumpData(BalanceBlockData):
    """
    Basic isentropic 0D turbine model.  This inherits the heater block to get
    a lot of unit model boilerplate and the mass balance, enegy balance and
    pressure equations.  This model is intended to be used only with Helmholtz
    EOS property pacakges in mixed or single phase mode with P-H state vars.

    Since this inherits BalanceBlockData, and only operates in steady-state or
    pseudo-steady-state (for dynamic models) the following mass, energy and
    pressure equations are implicitly writen.

    1) Mass Balance:
        0 = flow_mol_in[t] - flow_mol_out[t]
    2) Energy Balance:
        0 = (flow_mol[t]*h_mol[t])_in - (flow_mol[t]*h_mol[t])_out + Q_in + W_in
    3) Pressure:
        0 = P_in[t] + deltaP[t] - P_out[t]
    """

    CONFIG = BalanceBlockData.CONFIG()
    # For dynamics assume pseudo-steady-state
    CONFIG.dynamic = False
    CONFIG.get("dynamic")._default = False
    CONFIG.get("dynamic")._domain = In([False])
    CONFIG.has_holdup = False
    CONFIG.get("has_holdup")._default = False
    CONFIG.get("has_holdup")._domain = In([False])
    # Rest of config to make this function like a turbine
    CONFIG.has_pressure_change = True
    CONFIG.get("has_pressure_change")._default = True
    CONFIG.get("has_pressure_change")._domain = In([True])
    CONFIG.has_work_transfer = True
    CONFIG.get("has_work_transfer")._default = True
    CONFIG.get("has_work_transfer")._domain = In([True])
    CONFIG.has_heat_transfer = False
    CONFIG.get("has_heat_transfer")._default = False
    CONFIG.get("has_heat_transfer")._domain = In([False])

    def build(self):
        """
        Add model equations to the unit model.  This is called by a default block
        construnction rule when the unit model is created.
        """
        super().build()  # Basic unit model build/read config
        config = self.config  # shorter config pointer

        # The thermodynamic expression writer object, te, writes expressions
        # including external function calls to calculate thermodynamic quantities
        # from a set of state variables.
        _assert_properties(config.property_package)
        te = ThermoExpr(blk=self, parameters=config.property_package)

        eff = self.efficiency_pump = pyo.Var(
            self.flowsheet().time, initialize=0.9, doc="Pump efficiency"
        )
        self.efficiency_isentropic = pyo.Reference(self.efficiency_pump[:])

        pratio = self.ratioP = pyo.Var(
            self.flowsheet().time,
            initialize=0.7,
            doc="Ratio of outlet to inlet pressure",
        )

        # Some shorter refernces to property blocks
        properties_in = self.control_volume.properties_in
        properties_out = self.control_volume.properties_out

        @self.Expression(self.flowsheet().time, doc="Thermodynamic work")
        def work_fluid(b, t):
            return properties_out[t].flow_vol * (self.deltaP[t])

        @self.Expression(self.flowsheet().time, doc="Work required to drive the pump.")
        def shaft_work(b, t):  # Early access to the outlet enthalpy and work
            return self.work_fluid[t] / eff[t]

        @self.Constraint(self.flowsheet().time)
        def eq_work(b, t):  # outlet enthalpy coens from energy balance
            return self.control_volume.work[t] == self.shaft_work[t]

        @self.Constraint(self.flowsheet().time)
        def eq_pressure_ratio(b, t):
            return pratio[t] * properties_in[t].pressure == properties_out[t].pressure

    def initialize_build(
        self,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        For simplicity this initialization requires you to set values for the
        efficency, inlet, and one of pressure ratio, pressure change or outlet
        pressure.
        """
        init_log = idaeslog.getInitLogger(self.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(self.name, outlvl, tag="unit")

        # Create solver
        slvr = get_solver(solver, optarg)

        # Store original specification so initialization doesn't change the model
        # This will only resore the values of varaibles that were originally fixed
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)
        # Check for alternate pressure specs
        for t in self.flowsheet().time:
            if self.outlet.pressure[t].fixed:
                self.ratioP[t] = pyo.value(
                    self.outlet.pressure[t] / self.inlet.pressure[t]
                )
            elif self.control_volume.deltaP[t].fixed:
                self.ratioP[t] = pyo.value(
                    (self.control_volume.deltaP[t] + self.inlet.pressure[t])
                    / self.inlet.pressure[t]
                )
        # Fix the variables we base the initializtion on and free the rest.
        # This requires good values to be provided for pressure, efficency,
        # and inlet conditions, but it is simple and reliable.
        self.inlet.fix()
        self.outlet.unfix()
        self.ratioP.fix()
        self.deltaP.unfix()
        self.efficiency_pump.fix()
        for t in self.flowsheet().time:
            self.outlet.pressure[t] = pyo.value(self.inlet.pressure[t] * self.ratioP[t])
            self.deltaP[t] = pyo.value(self.outlet.pressure[t] - self.inlet.pressure[t])
            self.outlet.enth_mol[t] = pyo.value(
                self.inlet.enth_mol[t] + self.shaft_work[t] / self.inlet.flow_mol[t]
            )
            self.control_volume.work[t] = pyo.value(self.shaft_work[t])
            self.outlet.flow_mol[t] = pyo.value(self.inlet.flow_mol[t])
        # Solve the model (should be already solved from above)
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = slvr.solve(self, tee=slc.tee)
        from_json(self, sd=istate, wts=sp)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        for t, c in self.eq_pressure_ratio.items():
            s = iscale.get_scaling_factor(self.control_volume.properties_in[t].pressure)
            iscale.constraint_scaling_transform(c, s, overwrite=False)
        for t, c in self.eq_work.items():
            s = iscale.get_scaling_factor(self.control_volume.work[t])
            iscale.constraint_scaling_transform(c, s, overwrite=False)
