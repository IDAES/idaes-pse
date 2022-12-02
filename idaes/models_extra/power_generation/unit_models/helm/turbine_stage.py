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
Steam turbine stage model. This is a standard isentropic turbine. Under off-design
conditions the base efficiency and pressure ratio do not change much for the stages
between the inlet and outlet. This model is based on:

Liese, (2014). "Modeling of a Steam Turbine Including Partial Arc Admission
    for Use in a Process Simulation Software Environment." Journal of Engineering
    for Gas Turbines and Power. v136.
"""
__Author__ = "John Eslick"

from pyomo.environ import Var, units as pyunits

from idaes.core import declare_process_block_class
from idaes.models_extra.power_generation.unit_models.helm.turbine import (
    HelmIsentropicTurbineData,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class("HelmTurbineStage", doc="Basic steam turbine model")
class HelmTurbineStageData(HelmIsentropicTurbineData):
    CONFIG = HelmIsentropicTurbineData.CONFIG()

    def build(self):
        super().build()

        self.efficiency_mech = Var(initialize=1.0, doc="Turbine mechanical efficiency")
        self.efficiency_mech.fix()
        time_set = self.flowsheet().time
        self.shaft_speed = Var(
            time_set, doc="Shaft speed [1/s]", initialize=60.0, units=pyunits.s**-1
        )
        self.shaft_speed.fix()

        @self.Expression(time_set, doc="Specific speed [dimensionless]")
        def specific_speed(b, t):
            s = b.shaft_speed[t]  # 1/s
            v = b.control_volume.properties_out[t].flow_vol  # m3/s
            his_rate = b.work_isentropic[t]  # J/s
            m = b.control_volume.properties_out[t].flow_mass  # kg/s
            return s * v**0.5 * (his_rate / m) ** (-0.75)  # dimensionless

        @self.Expression(time_set, doc="Thermodynamic power [J/s]")
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
    ):
        """
        Initialize the turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves.

        Args:
            outlvl : sets output level of initialization routine
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        super().initialize_build(outlvl=outlvl, solver=solver, optarg=optarg)

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
