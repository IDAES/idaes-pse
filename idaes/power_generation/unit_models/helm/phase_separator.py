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
Simplified Flash Unit Model, only for IAPWS with mixed state.
Phase separator - inlet water/steam mixture is separated into liquid and vapor
streams.

Created: August 21 2020
"""
# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from pyomo.environ import SolverFactory, value
from idaes.core.util.initialization import fix_state_vars, revert_state_vars

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa, A. Lee)"


@declare_process_block_class("HelmPhaseSeparator")
class WaterFlashData(UnitModelBlockData):
    """
Simplified Flash Unit Model Class, only for IAPWS with mixed state
    """
    CONFIG = ConfigBlock()
    CONFIG.declare("dynamic", ConfigValue(
        default=False,
        domain=In([False]),
        description="Dynamic model flag",
        doc="""Indicates whether the model is dynamic"""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=False,
        domain=In([False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """
        Args:
            None

        Returns:
            None
        """
        super(WaterFlashData, self).build()

        self.mixed_state = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            default=self.config.property_package_args)

        self.add_port("inlet", self.mixed_state)

        self.vap_state = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            default=self.config.property_package_args)

        self.liq_state = self.config.property_package.build_state_block(
            self.flowsheet().config.time,
            default=self.config.property_package_args)

        self.add_port("vap_outlet", self.vap_state)
        self.add_port("liq_outlet", self.liq_state)
        # vapor outlet state
        @self.Constraint(self.flowsheet().config.time)
        def vap_material_balance(b, t):
            return 1e-4*b.mixed_state[t].flow_mol*b.mixed_state[t].vapor_frac == \
                b.vap_state[t].flow_mol*1e-4

        @self.Constraint(self.flowsheet().config.time)
        def vap_enthalpy_balance(b, t):
            return b.mixed_state[t].enth_mol_phase["Vap"]*1e-4 == \
                b.vap_state[t].enth_mol*1e-4

        @self.Constraint(self.flowsheet().config.time)
        def vap_pressure_balance(b, t):
            return b.mixed_state[t].pressure*1e-6 == \
                b.vap_state[t].pressure*1e-6

        # liquid outlet state
        @self.Constraint(self.flowsheet().config.time)
        def liq_material_balance(b, t):
            return 1e-4*b.mixed_state[t].flow_mol*(1 - b.mixed_state[t].vapor_frac)\
                == b.liq_state[t].flow_mol*1e-4

        @self.Constraint(self.flowsheet().config.time)
        def liq_enthalpy_balance(b, t):
            return 1e-4*b.mixed_state[t].enth_mol_phase["Liq"] == \
                b.liq_state[t].enth_mol*1e-4

        @self.Constraint(self.flowsheet().config.time)
        def liq_pressure_balance(b, t):
            return b.mixed_state[t].pressure*1e-6 == \
                b.liq_state[t].pressure*1e-6

    def initialize(blk, state_args_water_steam={},
                   outlvl=0, solver='ipopt', optarg={'tol': 1e-6}):
        '''
        Drum initialization routine.

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) for the control_volume of the model to
                           provide an initial state for initialization
                           (see documentation of the specific property package)
                           (default = None).
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        opt = SolverFactory(solver)
        opt.options = optarg

        init_log.info_low("Starting initialization...")
        # fix FeedWater Inlet
        flags_fw = fix_state_vars(blk.mixed_state,
                                  state_args_water_steam)
        blk.mixed_state.initialize()
        # initialize outlet states
        for t in blk.flowsheet().config.time:
            blk.vap_state[t].flow_mol = value(blk.mixed_state[t].
                                              flow_mol*blk.mixed_state[t].
                                              vapor_frac)
            blk.vap_state[t].enth_mol = value(blk.mixed_state[t].
                                              enth_mol_phase["Vap"])
            blk.vap_state[t].pressure = value(blk.mixed_state[t].pressure)
            blk.vap_state[t].vapor_frac = 1
            blk.liq_state[t].flow_mol = value(blk.mixed_state[t].flow_mol
                                              * (1 - blk.mixed_state[t].
                                                 vapor_frac))
            blk.liq_state[t].enth_mol = value(blk.mixed_state[t].
                                              enth_mol_phase["Liq"])
            blk.liq_state[t].pressure = value(blk.mixed_state[t].pressure)
            blk.liq_state[t].vapor_frac = 0
        # unfix variables
        revert_state_vars(blk.mixed_state, flags_fw)
        init_log.info_low("Initialization Complete.")

    def set_initial_condition(self):
        pass

    def calculate_scaling_factors(self):
        pass
