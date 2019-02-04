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
from __future__ import division

__Author__ = "John Eslick"

import logging
_log = logging.getLogger(__name__)

from pyomo.environ import Var, Expression, Constraint, sqrt

@declare_process_block_class("TurbineInletStage")
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
    def buid(self):
        super(TurbineInletStageData, self).build()

        self.flow_coeff = Var(self.time, initialize=1.053/3600.0,
            doc="Turbine flow coefficient [kg*C^0.5/Pa/s]")
        self.delta_enth_isentropic = Var(self.time, initialize=2.2e5,
            doc="Specific enthalpy change of isentropic process [J/mol/K]")
        self.blade_reaction = Var(initialize=0.9,
            doc="Blade reaction parameter")
        self.blade_velocity = Var(initialize=110.0,
            doc="Design blade velocity [m/s]")
        self.eff_nozzle = Var(initialize=0.95, bounds=(0.0, 1.0),
            doc="Nozzel efficiency (typically 0.90 to 0.95)")
        self.eff_nozzle.fix()
        self.blade_reaction()
        self.flow_coeff.fix()
        self.blade_velocity.fix()

        @self.Expression(self.time,
            doc="Entering steam velocity calculation [m/s]")
        def steam_entering_velocity(b, t):
            # 1.414 = 44.72/sqrt(1000) for SI if comparing to Liese (2014)
            # b.delta_enth_isentropic[t] = hin - hiesn, the mw converts enthalpy
            # to a mass basis
            return (1.414*sqrt((1-b.blade_reaction)*(b.delta_enth_isentropic[t]))/
                    sqrt(b.holdup.properties_in[t].mw))

        @self.Constraint(self.time, doc="Turbine inlet flow constraint")
        def inlet_flow_constraint(b, t):
            # Some local vars to make the euqation more readable
            cp_ratio = b.control_volume.properties_in[t].heat_capacity_ratio
            mw = b.control_volume.properties_in[t].mw
            flow = b.control_volume.properties_in[t].flow_mol
            Tin = b.control_volume.properties_in[t].temperature
            cf = b.flow_coeff[t]
            Pin = b.control_volume.properties_in[t].pressure
            Pratio = b.ratioP[t]
            return (flow*mw*sqrt(Tin - 273.15) ==
                Cf*Pin*sqrt(g/(g - 1)*(Pratio**(2/g) - Pratio**((g + 1)/g))))

            @self.Constraint(self.time,
                doc="Calculation of isentropic specific enthalpy change")
                def isentropic_enthalpy(b, t):
                    return b.work_isentropic[t] == (-b.delta_enth_isentropic[t]*
                        b.holdup.properties_in[t].flow_mol)

            @self.Constraint(self.time, doc="Efficiency correaltion")
                def efficiency_correlation(b, t):
                Vr = b.blade_velocity/b.steam_entering_velocity[t]
                eff = b.efficiency_isentropic[t]
                R = b.blade_reaction
                return eff == 2*Vr*((sqrt(1 - R) - Vr) +
                                     sqrt((sqrt(1 - R) - Vr)**2 + R))
