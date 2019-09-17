##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Heat Exchanger Models.
"""

__author__ = "John Eslick"

import logging
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (Var, Param, log, Expression, Constraint, Reference,
                           PositiveReals, SolverFactory, ExternalFunction,
                           exp, log10)

from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.opt import TerminationCondition

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)


from idaes.functions import functions_lib
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.unit_models.heater import (_make_heater_config_block,
                                      _make_heater_control_volume)

_log = logging.getLogger(__name__)


class HeatExchangerFlowPattern(Enum):
    countercurrent = 1
    cocurrent = 2
    crossflow = 3


def _make_heat_exchanger_config(config):
    """
    Declare configuration options for HeatExchangerData block.
    """
    config.declare("shell", ConfigBlock(
        implicit=True,
        description="Config block for shell",
        doc="""A config block used to construct the shell control volume."""))
    config.declare("tube", ConfigBlock(
        implicit=True,
        description="Config block for tube",
        doc="""A config block used to construct the tube control volume."""))
    _make_heater_config_block(config.shell)
    _make_heater_config_block(config.tube)
    config.declare("delta_temperature_callback", ConfigValue(
        default=delta_temperature_lmtd_callback,
        description="Callback for for temperature difference calculations"))
    config.declare("flow_pattern", ConfigValue(
        default=HeatExchangerFlowPattern.countercurrent,
        domain=In(HeatExchangerFlowPattern),
        description="Heat exchanger flow pattern",
        doc="""Heat exchanger flow pattern,
**default** - HeatExchangerFlowPattern.countercurrent.
**Valid values:** {
**HeatExchangerFlowPattern.countercurrent** - countercurrent flow,
**HeatExchangerFlowPattern.cocurrent** - cocurrent flow,
**HeatExchangerFlowPattern.crossflow** - cross flow, factor times
countercurrent temperature difference.}"""))


def delta_temperature_lmtd_callback(b):
    """
    This is a callback for a temperaure difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD).  It can be supplied to "delta_temperature_callback"
    HeatExchanger configuration option.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out
    @b.Expression(b.flowsheet().config.time)
    def delta_temperature(b, t):
        return (dT1[t] - dT2[t]) / log(dT1[t] / dT2[t])


def delta_temperature_amtd_callback(b):
    """
    This is a callback for a temperaure difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using arithmetic-mean
    temperature difference (AMTD).  It can be supplied to
    "delta_temperature_callback" HeatExchanger configuration option.
    """
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out
    @b.Expression(b.flowsheet().config.time)
    def delta_temperature(b, t):
        return (dT1[t] + dT2[t]) * 0.5


def delta_temperature_underwood_callback(b):
    """
    This is a callback for a temperaure difference expression to calculate
    :math:`\Delta T` in the heat exchanger model using log-mean temperature
    difference (LMTD) approximation given by Underwood (1970).  It can be
    supplied to "delta_temperature_callback" HeatExchanger configuration option.
    This uses a cube root function that works with negative numbers returning
    the real negative root. This should always evaluate successfully.
    """
    # external function that ruturns the real root, for the cuberoot of negitive
    # numbers, so it will return without error for positive and negitive dT.
    b.cbrt = ExternalFunction(library=functions_lib(), function="cbrt")
    dT1 = b.delta_temperature_in
    dT2 = b.delta_temperature_out
    @b.Expression(b.flowsheet().config.time)
    def delta_temperature(b, t):
        return ((b.cbrt(dT1[t]) + b.cbrt(dT2[t])) / 2.0)**3


@declare_process_block_class("HeatExchanger",
                             doc="Simple 0D heat exchanger model.")
class HeatExchangerData(UnitModelBlockData):
    """
    Simple 0D heat exchange unit.
    Unit model to transfer heat from one material to another.
    """
    CONFIG = UnitModelBlockData.CONFIG()
    _make_heat_exchanger_config(CONFIG)

    def set_scaling_factor_energy(self, f):
        """
        This function sets scaling_factor_energy for both shell and tube.
        This factor multiplies the energy balance and heat transfer equations
        in the heat exchnager.  The value of this factor should be about
        1/(expected heat duty).

        Args:
            f: Energy balance scaling factor
        """
        self.shell.scaling_factor_energy.value = f
        self.tube.scaling_factor_energy.value = f

    def build(self):
        """
        Building model

        Args:
            None
        Returns:
            None
        """
        ########################################################################
        #  Call UnitModel.build to setup dynamics and configure                #
        ########################################################################
        super().build()
        config = self.config
        ########################################################################
        # Add variables                                                        #
        ########################################################################
        u = self.overall_heat_transfer_coefficient = Var(
            self.flowsheet().config.time, domain=PositiveReals,
            initialize=100.0, doc="Overall heat transfer coefficient")
        a = self.area = Var(
            domain=PositiveReals, initialize=1000.0, doc="Heat exchange area")
        self.delta_temperature_in = Var(
            self.flowsheet().config.time, initialize=10.0,
            doc="Temperature difference at the shell inlet end")
        self.delta_temperature_out = Var(
            self.flowsheet().config.time, initialize=10.0,
            doc="Temperature difference at the shell outlet end")
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            self.crossflow_factor = Var(
                self.flowsheet().config.time, initialize=1.0,
                doc="Factor to adjust coutercurrent flow heat "
                "transfer calculation for cross flow.")
            f = self.crossflow_factor
        ########################################################################
        # Add control volumes                                                  #
        ########################################################################
        _make_heater_control_volume(self, "shell", config.shell,
                                    dynamic=config.dynamic,
                                    has_holdup=config.has_holdup)
        _make_heater_control_volume(self, "tube", config.tube,
                                    dynamic=config.dynamic,
                                    has_holdup=config.has_holdup)
        # Add convienient references to heat duty.
        q = self.heat_duty = Reference(self.tube.heat)
        ########################################################################
        # Add ports                                                            #
        ########################################################################
        self.add_inlet_port(name="inlet_1", block=self.shell)
        self.add_inlet_port(name="inlet_2", block=self.tube)
        self.add_outlet_port(name="outlet_1", block=self.shell)
        self.add_outlet_port(name="outlet_2", block=self.tube)
        ########################################################################
        # Add end temperaure differnece constraints                            #
        ########################################################################
        @self.Constraint(self.flowsheet().config.time)
        def delta_temperature_in_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (b.delta_temperature_in[t] ==
                        b.shell.properties_in[t].temperature -
                        b.tube.properties_in[t].temperature)
            else:
                return (b.delta_temperature_in[t] ==
                        b.shell.properties_in[t].temperature -
                        b.tube.properties_out[t].temperature)

        @self.Constraint(self.flowsheet().config.time)
        def delta_temperature_out_equation(b, t):
            if b.config.flow_pattern == HeatExchangerFlowPattern.cocurrent:
                return (b.delta_temperature_out[t] ==
                        b.shell.properties_out[t].temperature -
                        b.tube.properties_out[t].temperature)
            else:
                return (b.delta_temperature_out[t] ==
                        b.shell.properties_out[t].temperature -
                        b.tube.properties_in[t].temperature)
        ########################################################################
        # Add a unit level energy balance                                      #
        ########################################################################
        @self.Constraint(self.flowsheet().config.time)
        def unit_heat_balance(b, t):
            return 0 == self.shell.heat[t] + self.tube.heat[t]
        ########################################################################
        # Add delta T calculations using callack function, lots of options,    #
        #   and users can provide their own if needed                          #
        ########################################################################
        config.delta_temperature_callback(self)
        ########################################################################
        # Add Heat transfer equation                                           #
        ########################################################################
        deltaT = self.delta_temperature
        scale = self.shell.scaling_factor_energy
        @self.Constraint(self.flowsheet().config.time)
        def heat_transfer_equation(b, t):
            if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
                return 0 == (f[t] * u[t] * a * deltaT[t] - q[t]) * scale
            else:
                return 0 == (u[t] * a * deltaT[t] - q[t]) * scale
        ########################################################################
        # Add symbols for LaTeX equation rendering                             #
        ########################################################################
        self.overall_heat_transfer_coefficient.latex_symbol = "U"
        self.area.latex_symbol = "A"
        self.shell.heat.latex_symbol = "Q_1"
        self.tube.heat.latex_symbol = "Q_2"
        self.delta_temperature.latex_symbol = "\\Delta T"

    def initialize(self, state_args_1=None, state_args_2=None, outlvl=0,
                   solver='ipopt', optarg={'tol': 1e-6}, duty=1000):
        """
        Heat exchanger initialization method.

        Args:
            state_args_1 : a dict of arguments to be passed to the property
                initialization for shell (see documentation of the specific
                property package) (default = {}).
            state_args_2 : a dict of arguments to be passed to the property
                initialization for tube (see documentation of the specific
                property package) (default = {}).
            outlvl : sets output level of initialisation routine
                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')
            duty : an initial guess for the amount of heat transfered
                (default = 10000)

        Returns:
            None

        """
        # Set solver options
        tee = True if outlvl >= 3 else False
        opt = SolverFactory(solver)
        opt.options = optarg
        flags1 = self.shell.initialize(outlvl=outlvl - 1,
                                       optarg=optarg,
                                       solver=solver,
                                       state_args=state_args_1)

        if outlvl > 0:
            _log.info('{} Initialization Step 1a (shell) Complete.'
                      .format(self.name))

        flags2 = self.tube.initialize(outlvl=outlvl - 1,
                                      optarg=optarg,
                                      solver=solver,
                                      state_args=state_args_2)

        if outlvl > 0:
            _log.info('{} Initialization Step 1b (tube) Complete.'
                      .format(self.name))
        # ---------------------------------------------------------------------
        # Solve unit without heat transfer equation
        self.heat_transfer_equation.deactivate()
        self.tube.heat.fix(duty)
        results = opt.solve(self, tee=tee, symbolic_solver_labels=True)
        if outlvl > 0:
            if results.solver.termination_condition == \
                    TerminationCondition.optimal:
                _log.info('{} Initialization Step 2 Complete.'
                          .format(self.name))
            else:
                _log.warning('{} Initialization Step 2 Failed.'
                             .format(self.name))
        self.tube.heat.unfix()
        self.heat_transfer_equation.activate()
        # ---------------------------------------------------------------------
        # Solve unit
        results = opt.solve(self, tee=tee, symbolic_solver_labels=True)
        if outlvl > 0:
            if results.solver.termination_condition == \
                    TerminationCondition.optimal:
                _log.info('{} Initialization Step 3 Complete.'
                          .format(self.name))
            else:
                _log.warning('{} Initialization Step 3 Failed.'
                             .format(self.name))
        # ---------------------------------------------------------------------
        # Release Inlet state
        self.shell.release_state(flags1, outlvl - 1)
        self.tube.release_state(flags2, outlvl - 1)

        if outlvl > 0:
            _log.info('{} Initialization Complete.'.format(self.name))

    def _get_performance_contents(self, time_point=0):
        var_dict = {"HX Coefficient":
                    self.overall_heat_transfer_coefficient[time_point]}
        var_dict["HX Area"] = self.area
        var_dict["Heat Duty"] = self.heat_duty[time_point]
        if self.config.flow_pattern == HeatExchangerFlowPattern.crossflow:
            var_dict = {"Crossflow Factor": self.crossflow_factor[time_point]}

        expr_dict = {}
        expr_dict["Delta T Driving"] = self.delta_temperature[time_point]
        expr_dict["Delta T In"] = self.delta_temperature_in[time_point]
        expr_dict["Delta T Out"] = self.delta_temperature_out[time_point]

        return {"vars": var_dict, "exprs": expr_dict}

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
                {"Shell Inlet": self.inlet_1,
                 "Shell Outlet": self.outlet_1,
                 "Tube Inlet": self.inlet_2,
                 "Tube Outlet": self.outlet_2},
                 time_point=time_point)

    def _get_costing(self, hx_type='U-tube', FM='stain_steel', ce_i= '2018', 
                     FL = '12ft'):
        '''
        Heat exchanger costing method
        Source: Process and Product Design Principles: Synthesis, Analysis, 
        and Evaluation
        Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
        Chapter 22. Cost Accounting and Capital Cost Estimation
        22.2 Cost Indexes and Capital Investment
        
        This method computes the purchase cost (CP) for a shell and tube heat
        exchanger (Eq. 22.43), the model computes the base cost (CB for 4 type 
        of heat exchangers, such as floating head, fixed head, U-tube, and 
        Kettle vaporizer), construction material factor (FM), pressure design 
        factor (FP), and tube lenght correction factor (FL), using CE base cost 
        index of 500. 
        
        Cp = FP*FM*FL*CB
        
        Args:
            hx_type : 'floating_head', 'fixed_head', 'U-tube'*, 'Kettle_vap'
            
            material factor (FM): 'stain_steel'*, 'carb_steel'
            
            tube length (FL): '8ft'*, '12ft', '16ft', '20ft'
            
            CE Cost Index (ce_i): '2019', '2018'*, '2017' ... '2010'
            
            * --> default options
        '''
        # ------------------------ Heat Exchanger cost ------------------------
        # heat exchanger cost
        self.base_cost = Var(initialize = 1e5,
                             bounds=(0, 1e12),
                             doc='Heat Exchanger Base Cost CB cost in $')
        self.purchase_cost = Var(initialize =1e6,
                                 bounds = (0,1e12),
                                 doc = 'Heat Exchanger Purchase Cost CP in $')
        self.FP = Var(initialize = 100,
                      doc = 'hx Pressure factor')
        self.FL = Param(mutable = True, initialize = 1.12,
                        doc='hx tube length correction factor')
        self.FM = Param(mutable = True, initialize = 3,
                        doc= 'construction material correction factor')
        self.CE_index = Param(mutable = True, initialize = 671.1,
                              doc= 'CE cost index $ year')
        self.hx_os = Param(mutable = True, initialize = 1.1,
                           doc= 'hx oversize factor 1.1 to 1.5')

        # select length correction factor 
        c_fl = {'8ft':1.25, '12ft':1.12, '16ft':1.05, '20ft':1.00}
        self.FL = c_fl[FL]

        # select construction material (stainless steel default)
        c_material = {'stain_steel':3, 'carbon_steel':2}
        self.FM = c_material[FM]
        # ToDo (9/17/2019):FM can be calculated for advanced materials Eq.22.44 
        # i.e. Carbon Steel/Cr-Mo steel, Cr-Mo steel/Cr-Mo steel, Monel/Monel, 
        # Titanium/titanium, etc.
        
        # Cost index $/year (method argument or 2018 default)
        ce_index_dic = {'2019':680, '2018':671.1, '2017':567.5, '2016':541.7,
                    '2015':556.8, '2014':576.1, '2013':567.3, '2012':584.6,
                    '2011':585.7, '2010':550.8}
        self.CE_index = ce_index_dic[ce_i]

        #--------------------------------------------------
        # base cost calculation
#       # select heat exchanger type:
        alf1 = {'floating_head':11.9052,'fixed_head':11.2927,
                'U-tube':11.3852,'Kettle_vap':12.2052}
        alf2 = {'floating_head':0.8709,'fixed_head':0.8228,
                'U-tube':0.9186,'Kettle_vap':0.8709}
        alf3 = {'floating_head':0.09005,'fixed_head':0.09861,
                'U-tube':0.09790,'Kettle_vap':0.09005}
        
        if (self.config.tube.property_package.get_metadata().
                default_units['length']) == 'm':
            area = self.area*10.7639
        elif (self.config.tube.property_package.get_metadata().
                default_units['length']) == 'ft':
            area = self.area
        else:
            raise Exception('area units not suported')

        def hx_cost_rule(self):
            return self.base_cost == self.FL * exp(alf1[hx_type]
                        - alf2[hx_type]*log10(area*self.hx_os)
                        + alf3[hx_type]*log10(area*self.hx_os)**2)
        self.base_cost_eq = Constraint(rule=hx_cost_rule)
        
        #------------------------------------------------------
        # Pressure factor calculation
        # doublecheck units (higher pressure fluid should be tube side)
        if (self.config.tube.property_package.get_metadata().
                properties['pressure']['units']) == 'Pa':
            pressure = self.tube.properties_in[0].pressure*14.69/1.01325e5
        elif (self.config.tube.property_package.get_metadata().
                properties['pressure']['units']) == 'psig':
            pressure = self.tube.properties_in[0].pressure

        #units must be in psig 
        #(ToDo(9/17/2019): check inlet pressure to select correlation)
        def hx_P_factor(self):
            # eq valid from 600 pisg to 3000 psig
            #    return self.hx_FP == 0.8510 + 0.1292*(pressure/600) 
            #                                + 0.0198*(pressure/600)**2
            # eq valid from 100 pisg to 2000 psig
            return self.FP == 0.9803 + (0.018*(pressure/100)
                                        + 0.0017*(pressure/100)**2)
        self.p_factor_eq = Constraint(rule=hx_P_factor)


        #---------------------------------------------------------
        # purchase cost equation
        def hx_CP_rule(self):
            return self.purchase_cost == (self.FP*self.FM*self.FL*
                                          (self.CE_index/500)*self.base_cost)
        self.cp_cost_eq = Constraint(rule=hx_CP_rule)
