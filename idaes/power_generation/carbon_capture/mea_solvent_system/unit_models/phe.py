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
IDAES PLATE HEAT EXCHANGER MODEL (PHE) USING EFFECTIVENESS-NTU METHOD

Schematic diagram

hot fluid----->|                                                   |<-- cold fluid
(even channels)|      Pass 1   |     Pass 2 to P-1 |    Pass P     |(odd channels)
               |    ___     ___|                   |    ___     ___|
               |   |   |   |   |                   |   |   |   |   |
               |   |   |   |   |                   |   |   |   |   |
               |   |   |   |   |                   |   |   |   |   |
               |   |___|   |___|                   |   |___|   |___|
               |               |                   |               |
channels         1   2   3  Np                       1   2   3  Np

P  : number of passes(same for both fluids)
if P is even: PHE is in counter-current mode
if P is odd : PHE  is in co-current mode
Np : number of channels in a pass(should be even for balanced PHE)
Hot and cold fluids flow alternatively in the channels.

Note: Total Heat transfer area of the PHE can be increased by adding more plates,
which is equilvalent to increasing the number of channels in a pass(Np)

Detailed model equations can be found in the paper :
Akula, P., Eslick, J., Bhattacharyya, D. and Miller, D.C., 2019.
Modelling and Parameter Estimation of a Plate Heat Exchanger
as Part of a Solvent-Based Post-Combustion CO2 Capture System.
In Computer Aided Chemical Engineering (Vol. 47, pp. 47-52). Elsevier.

"""

# Import Pyomo libraries
from pyomo.environ import Param, RangeSet, Constraint, Expression,\
    SolverFactory, value, Var, exp, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

__author__ = "Paul Akula"


# Set up logger
_log = idaeslog.getLogger(__name__)


class plate_param(object):
    '''
    Details of a smooth plate.
    For Chevron plates, enlargement factor may be incorporated into the plate area
    '''
    L = 1.657               # plate length
    w = 0.833               # plate width
    b = 0.0038              # gap between two plates forming a channel
    kp = 16.2               # plate thermal conductivity
    bp = 0.0006             # plate thickness
    Dp = 0.3                # plate port diameter
    Ap = L * w              # plate area
    De = 2 * b              # channel equivalent diameter


@declare_process_block_class("PHE")
class PHEData(UnitModelBlockData):
    """
    PHE Unit Model Class
    """
    CONFIG = UnitModelBlockData.CONFIG()

    # Configuration template for fluid specific  arguments
    _SideCONFIG = ConfigBlock()

    CONFIG.declare("passes", ConfigValue(
        default=4,
        domain=int,
        description="Number of passes",
        doc="""Number of passes for both the hot and cold fluid"""))

    CONFIG.declare("channel_list", ConfigValue(
        default=[12, 12, 12, 12],
        domain=list,
        description="Number of channels for each pass",
        doc="""Number of channels to be used in each pass"""))

    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.componentTotal,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
            **default** - MaterialBalanceType.useDefault.
            **Valid values:** {
        **MaterialBalanceType.useDefault -default from property package
        **MaterialBalanceType.none** - exclude material balances,
        **MaterialBalanceType.componentPhase** - use phase component balances,
        **MaterialBalanceType.componentTotal** - use total component balances,
        **MaterialBalanceType.elementTotal** - use total element balances,
        **MaterialBalanceType.total** - use total material balance.}"""))

    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
        **default** - EnergyBalanceType.enthalpyTotal
        **Valid values:** {
        **EnergyBalanceType.useDefault - default from property package
        **EnergyBalanceType.none** - exclude energy balances,
        **EnergyBalanceType.enthalpyTotal**-single enthalpy balance for material,
        **EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
        **EnergyBalanceType.energyTotal** - single energy balance for material,
        **EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))

    CONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
        **default** - MomentumBalanceType.pressureTotal.
        **Valid values:** {
        **MomentumBalanceType.none** - exclude momentum balances,
        **MomentumBalanceType.pressureTotal**-single pressure balance for material,
        **MomentumBalanceType.pressurePhase**-pressure balances for each phase,
        **MomentumBalanceType.momentumTotal**-single momentum balance for material,
        **MomentumBalanceType.momentumPhase**-momentum balances for each phase.}"""))

    CONFIG.declare("has_pressure_change", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
            constructed,
        **default** - False.
        **Valid values:** {
        **True** - include pressure change terms,
        **False** - exclude pressure change terms.}"""))

    _SideCONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use ",
        doc="""Property parameter object used to define property calculations
        **default** - useDefault.
        **Valid values:** {
        **useDefault** - use default package from parent model or flowsheet,
        **PhysicalParameterObject** - a PhysicalParameterBlock object.}"""))

    _SideCONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property package",
        doc="""A ConfigBlock with arguments to be passed to
        property block(s) and used when constructing these,
        **default** - None.
        **Valid values:** {
        see property package for documentation.}"""))

    # Create individual config blocks for hot and cold sides
    CONFIG.declare("hot_side",
                   _SideCONFIG(doc="Hot fluid config arguments"))
    CONFIG.declare("cold_side",
                   _SideCONFIG(doc="Cold fluid config arguments"))

    def build(self):
        # Call UnitModel.build to setup model
        super(PHEData, self).build()

        # Consistency check for number of passes and channels in each pass
        for i in self.config.channel_list:
            if not isinstance(i, int):
                raise ConfigurationError("number of channels ({}) must be"
                                         " an integer".format(i))

        if (self.config.passes != len(self.config.channel_list)):
            raise ConfigurationError(
                "The number of elements in the channel list: {}  "
                " does not match the number of passes ({}) given. "
                "Please provide as integers, the number of channels of each pass".format(
                    self.config.channel_list, self.config.passes))

        # ======================================================================
        # Build hot-side  Control Volume (Lean Solvent)
        self.hot_side = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.hot_side.property_package,
            "property_package_args": self.config.hot_side.property_package_args})

        self.hot_side.add_state_blocks(has_phase_equilibrium=False)

        self.hot_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=False,
            has_phase_equilibrium=False,
            has_rate_reactions=False)

        self.hot_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # ======================================================================
        # Build cold-side  Control Volume(Rich solvent)
        self.cold_side = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.cold_side.property_package,
            "property_package_args": self.config.cold_side.property_package_args})

        self.cold_side.add_state_blocks(has_phase_equilibrium=False)

        self.cold_side.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_mass_transfer=False,
            has_phase_equilibrium=False,
            has_rate_reactions=False)

        self.cold_side.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        # ======================================================================
        # Add Ports to control volumes
        # hot-side
        self.add_inlet_port(name="hot_inlet",
                            block=self.hot_side, doc='inlet Port')
        self.add_outlet_port(name="hot_outlet",
                             block=self.hot_side, doc='outlet Port')

        # cold-side
        self.add_inlet_port(name="cold_inlet",
                            block=self.cold_side, doc='inlet Port')
        self.add_outlet_port(name="cold_outlet",
                             block=self.cold_side, doc='outlet Port')
        # ======================================================================
        # Add performace equation method
        self._make_params()
        self._make_performance_method()

    def _make_params(self):
        self.P = Param(initialize=self.config.passes,
                       units=None,
                       doc="Total number of passes for hot or cold fluid")
        self.PH = RangeSet(self.P,
                           doc="Set of hot fluid passes")
        self.PC = RangeSet(self.P,
                           doc="Set of cold fluid passes(equal to PH)")
        self.plate_thermal_cond = Param(initialize=plate_param.kp,
                                        units=pyunits.W / pyunits.m / pyunits.K,
                                        doc="Plate thermal conductivity")
        self.plate_thick = Param(initialize=plate_param.bp,
                                 units=pyunits.m,
                                 doc="Plate thickness")
        self.plate_gap = Param(initialize=plate_param.b,
                               units=pyunits.m,
                               doc="Plate gap")
        self.plate_length = Param(initialize=plate_param.L,
                                  units=pyunits.m,
                                  doc="Plate length ")
        self.plate_width = Param(initialize=plate_param.w,
                                 units=pyunits.m,
                                 doc="Plate width ")
        self.plate_area = Param(initialize=plate_param.Ap,
                                units=pyunits.m**2,
                                doc="Area ")
        self.channel_dia = Param(initialize=plate_param.De,
                                 units=pyunits.m,
                                 doc=" Channel equivalent diameter")
        self.port_dia = Param(initialize=plate_param.Dp,
                              units=pyunits.m,
                              doc=" Port diameter of plate ")
        self.Np = Param(self.PH,
                        units=None,
                        doc="Number of channels in each pass", mutable=True)
        for i in self.PH:
            self.Np[i] = self.config.channel_list[i - 1]

        # heat transfer parameters
        self.param_a = Var(initialize=0.3, bounds=(0.2, 0.4), units=None,
                           doc='Nusselt parameter')
        self.param_b = Var(initialize=0.663, bounds=(0.3, 0.7), units=None,
                           doc='Nusselt parameter')
        self.param_c = Var(initialize=1 / 3.0, bounds=(1e-5, 2), units=None,
                           doc='Nusselt parameter')
        self.param_a.fix(0.4)
        self.param_b.fix(0.663)
        self.param_c.fix(0.333)

    def _make_performance_method(self):

        # Add object references - Sets
        add_object_reference(self, "component_list_ref",
                             self.config.hot_side.property_package.component_list)

        add_object_reference(
            self,
            "solvent_list",
            self.config.hot_side.property_package.component_list_solvent)

        add_object_reference(self, "cp_ref",
                             self.config.hot_side.property_package.cp_param)

        def rule_trh(blk, t):
            return (blk.hot_side.properties_out[t].temperature /
                    blk.hot_side.properties_in[t].temperature)
        self.trh = Expression(self.flowsheet().config.time, rule=rule_trh,
                              doc='Ratio of hot outlet temperature to hot'
                                  'inlet temperature')

        def rule_trc(blk, t):
            return (blk.cold_side.properties_out[t].temperature /
                    blk.cold_side.properties_in[t].temperature)
        self.trc = Expression(self.flowsheet().config.time, rule=rule_trc,
                              doc='Ratio of cold outlet temperature to cold'
                                  ' inlet temperature')

        def rule_cp_comp_hot(blk, t, j):
            return 1e3 * (
                blk.cp_ref[j, 1] +
                blk.cp_ref[j, 2] / 2 * blk.hot_side.properties_in[t].temperature *
                (blk.trh[t] + 1) +
                blk.cp_ref[j, 3] / 3 * (blk.hot_side.properties_in[t].temperature**2) *
                (blk.trh[t]**2 + blk.trh[t] + 1) +
                blk.cp_ref[j, 4] / 4 * (blk.hot_side.properties_in[t].temperature**3) *
                (blk.trh[t] + 1) * (blk.trh[t]**2 + 1) +
                blk.cp_ref[j, 5] / 5 * (blk.hot_side.properties_in[t].temperature**4) *
                (blk.trh[t]**4 + blk.trh[t]**3 + blk.trh[t]**2 + blk.trh[t] + 1))

        self.cp_comp_hot = Expression(self.flowsheet().config.time,
                                      self.solvent_list,
                                      rule=rule_cp_comp_hot,
                                      doc='Component mean specific heat capacity'
                                      ' btw inlet and outlet'
                                      ' of hot-side temperature')

        def rule_cp_hot(blk, t):
            return sum(blk.cp_comp_hot[t, j] *
                       blk.hot_side.properties_in[t].mass_frac_co2_free[j]
                       for j in blk.solvent_list)
        self.cp_hot = Expression(self.flowsheet().config.time, rule=rule_cp_hot,
                                 doc='Hot-side mean specific heat capacity on'
                                     'free CO2 basis')

        def rule_cp_comp_cold(blk, t, j):
            return 1e3 * (
                blk.cp_ref[j, 1] +
                blk.cp_ref[j, 2] / 2 * blk.cold_side.properties_in[t].temperature *
                (blk.trc[t] + 1) +
                blk.cp_ref[j, 3] / 3 * (blk.cold_side.properties_in[t].temperature**2) *
                (blk.trc[t]**2 + blk.trc[t] + 1) +
                blk.cp_ref[j, 4] / 4 * (blk.cold_side.properties_in[t].temperature**3) *
                (blk.trc[t] + 1) * (blk.trc[t]**2 + 1) +
                blk.cp_ref[j, 5] / 5 * (blk.cold_side.properties_in[t].temperature**4) *
                (blk.trc[t]**4 + blk.trc[t]**3 + blk.trc[t]**2 + blk.trc[t] + 1))

        self.cp_comp_cold = Expression(self.flowsheet().config.time,
                                       self.solvent_list,
                                       rule=rule_cp_comp_cold,
                                       doc='Component mean specific heat capacity'
                                           'btw inlet and outlet'
                                           ' of cold-side temperature')

        def rule_cp_cold(blk, t):
            return sum(blk.cp_comp_cold[t, j] *
                       blk.cold_side.properties_in[t].mass_frac_co2_free[j]
                       for j in blk.solvent_list)
        self.cp_cold = Expression(self.flowsheet().config.time, rule=rule_cp_cold,
                                  doc='Cold-side mean specific heat capacity'
                                      'on free CO2 basis')

        # Model Variables
        self.Th_in = Var(self.flowsheet().config.time,
                         self.PH, initialize=393, units=pyunits.K,
                         doc="Hot Temperature IN of pass")
        self.Th_out = Var(self.flowsheet().config.time,
                          self.PH, initialize=325, units=pyunits.K,
                          doc="Hot Temperature OUT of pass")
        self.Tc_in = Var(self.flowsheet().config.time,
                         self.PH, initialize=320, units=pyunits.K,
                         doc="Cold Temperature IN of pass")
        self.Tc_out = Var(self.flowsheet().config.time,
                          self.PH, initialize=390, units=pyunits.K,
                          doc="Cold Temperature OUT of pass")

        # ======================================================================
        # PERFORMANCE EQUATIONS

        # ---------------MASS BALANCE CONSTRAINTS------------------------------
        def rule_eq_hotside_sumx(bk, t):
            return bk.hot_side.properties_out[t].flow_mol == \
                sum(bk.hot_side.properties_out[t].flow_mol_comp[j]
                    for j in bk.component_list_ref)
        self.eq_hotside_sumx = Constraint(self.flowsheet().config.time,
                                          rule=rule_eq_hotside_sumx)

        def rule_eq_coldside_sumx(bk, t):
            return bk.cold_side.properties_out[t].flow_mol == \
                sum(bk.cold_side.properties_out[t].flow_mol_comp[j]
                    for j in bk.component_list_ref)
        self.eq_coldside_sumx = Constraint(self.flowsheet().config.time,
                                           rule=rule_eq_coldside_sumx)

        # ----------------------------------------------------------------------
        # mass flow rate in kg/s
        def rule_mh_in(blk, t):
            return blk.hot_side.properties_in[t].flow_mol *\
                blk.hot_side.properties_in[t].mw
        self.mh_in = Expression(self.flowsheet().config.time, rule=rule_mh_in,
                                doc='Hotside mass flow rate [kg/s]')

        def rule_mc_in(blk, t):
            return blk.cold_side.properties_in[t].flow_mol *\
                blk.cold_side.properties_in[t].mw
        self.mc_in = Expression(self.flowsheet().config.time, rule=rule_mc_in,
                                doc='Coldside mass flow rate [kg/s]')

        # ----------------------------------------------------------------------
        # port mass velocity[kg/m2.s]
        def rule_Gph(blk, t):
            return (4 * blk.mh_in[t] * 7) / (22 * blk.port_dia**2)

        self.Gph = Expression(self.flowsheet().config.time, rule=rule_Gph,
                              doc='Hotside port mass velocity[kg/m2.s]')

        def rule_Gpc(blk, t):
            return (4 * blk.mc_in[t] * 7) / (22 * blk.port_dia**2)

        self.Gpc = Expression(self.flowsheet().config.time, rule=rule_Gpc,
                              doc='Coldside port mass velocity[kg/m2.s]')

        # ----------------------------------------------------------------------
        # Reynold & Prandtl numbers
        def rule_Re_h(blk, t, p):
            return blk.mh_in[t] * blk.channel_dia /\
                (blk.Np[p] * blk.plate_width *
                 blk.plate_gap * blk.hot_side.properties_in[t].visc_d)
        self.Re_h = Expression(self.flowsheet().config.time, self.PH,
                               rule=rule_Re_h,
                               doc='Hotside Reynolds number')

        def rule_Re_c(blk, t, p):
            return blk.mc_in[t] * blk.channel_dia /\
                (blk.Np[p] * blk.plate_width *
                 blk.plate_gap * blk.cold_side.properties_in[t].visc_d)
        self.Re_c = Expression(self.flowsheet().config.time, self.PH,
                               rule=rule_Re_c,
                               doc='Coldside Reynolds number')

        def rule_Pr_h(blk, t):
            return blk.cp_hot[t] * blk.hot_side.properties_in[t].visc_d /\
                blk.hot_side.properties_in[t].thermal_cond
        self.Pr_h = Expression(self.flowsheet().config.time,
                               rule=rule_Pr_h,
                               doc='Hotside Prandtl number')

        def rule_Pr_c(blk, t):
            return blk.cp_cold[t] * blk.cold_side.properties_in[t].visc_d /\
                blk.cold_side.properties_in[t].thermal_cond
        self.Pr_c = Expression(self.flowsheet().config.time,
                               rule=rule_Pr_c,
                               doc='Coldside Prandtl number')
        # ----------------------------------------------------------------------
        # Film heat transfer coefficients

        def rule_hotside_transfer_coef(blk, t, p):
            return (blk.hot_side.properties_in[t].thermal_cond / blk.channel_dia *
                    blk.param_a * blk.Re_h[t, p]**blk.param_b *
                    blk.Pr_h[t]**blk.param_c)

        self.h_hot = Expression(self.flowsheet().config.time,
                                self.PH, rule=rule_hotside_transfer_coef,
                                doc='Hotside heat transfer coefficient')

        def rule_coldside_transfer_coef(blk, t, p):
            return (blk.cold_side.properties_in[t].thermal_cond / blk.channel_dia *
                    blk.param_a * blk.Re_c[t, p]**blk.param_b *
                    blk.Pr_c[t]**blk.param_c)

        self.h_cold = Expression(self.flowsheet().config.time,
                                 self.PH, rule=rule_coldside_transfer_coef,
                                 doc='Coldside heat transfer coefficient')

        # ----------------------------------------------------------------------
        # Friction factor calculation
        def rule_fric_h(blk, t):
            return 18.29 * blk.Re_h[t, 1]**(-0.652)
        self.fric_h = Expression(self.flowsheet().config.time,
                                 rule=rule_fric_h,
                                 doc='Hotside friction factor')

        def rule_fric_c(blk, t):
            return 1.441 * self.Re_c[t, 1]**(-0.206)
        self.fric_c = Expression(self.flowsheet().config.time,
                                 rule=rule_fric_c,
                                 doc='Coldside friction factor')

        # ----------------------------------------------------------------------
        # pressure drop calculation
        def rule_hotside_dP(blk, t):
            return (2 * blk.fric_h[t] * (blk.plate_length + blk.port_dia) *
                    blk.P * blk.Gph[t]**2) /\
                (blk.hot_side.properties_in[t].dens_mass *
                 blk.channel_dia) + 1.4 * blk.P * blk.Gph[t]**2 * 0.5 /\
                blk.hot_side.properties_in[t].dens_mass + \
                blk.hot_side.properties_in[t].dens_mass * \
                9.81 * (blk.plate_length + blk.port_dia)

        self.dP_h = Expression(self.flowsheet().config.time,
                               rule=rule_hotside_dP,
                               doc='Hotside pressure drop  [Pa]')

        def rule_coldside_dP(blk, t):
            return (2 * blk.fric_c[t] * (blk.plate_length + blk.port_dia) *
                    blk.P * blk.Gpc[t]**2) /\
                (blk.cold_side.properties_in[t].dens_mass * blk.channel_dia) +\
                1.4 * (blk.P * blk.Gpc[t]**2 * 0.5 /
                       blk.cold_side.properties_in[t].dens_mass) + \
                blk.cold_side.properties_in[t].dens_mass * \
                9.81 * (blk.plate_length + blk.port_dia)

        self.dP_c = Expression(self.flowsheet().config.time,
                               rule=rule_coldside_dP,
                               doc='Coldside pressure drop  [Pa]')

        def rule_eq_deltaP_hot(blk, t):
            return blk.hot_side.deltaP[t] == -blk.dP_h[t]
        self.eq_deltaP_hot = Constraint(self.flowsheet().config.time,
                                        rule=rule_eq_deltaP_hot)

        def rule_eq_deltaP_cold(blk, t):
            return blk.cold_side.deltaP[t] == -blk.dP_c[t]
        self.eq_deltaP_cold = Constraint(self.flowsheet().config.time,
                                         rule=rule_eq_deltaP_cold)

        # ----------------------------------------------------------------------
        # Overall heat transfer coefficients
        def rule_U(blk, t, p):
            return 1.0 /\
                (1.0 / blk.h_hot[t, p] + blk.plate_gap / blk.plate_thermal_cond +
                 1.0 / blk.h_cold[t, p])
        self.U = Expression(self.flowsheet().config.time, self.PH,
                            rule=rule_U,
                            doc='Overall heat transfer coefficient')
        # ----------------------------------------------------------------------
        # capacitance of hot and cold fluid

        def rule_Caph(blk, t, p):
            return blk.mh_in[t] * blk.cp_hot[t] / blk.Np[p]
        self.Caph = Expression(self.flowsheet().config.time, self.PH,
                               rule=rule_Caph,
                               doc='Hotfluid capacitance rate')

        def rule_Capc(blk, t, p):
            return blk.mc_in[t] * blk.cp_cold[t] / blk.Np[p]
        self.Capc = Expression(self.flowsheet().config.time, self.PH,
                               rule=rule_Capc,
                               doc='Coldfluid capacitance rate')

        # ----------------------------------------------------------------------
        # min n max capacitance and capacitance ratio
        def rule_Cmin(blk, t, p):
            return 0.5 * (blk.Caph[t, p] + blk.Capc[t, p] -
                          ((blk.Caph[t, p] - blk.Capc[t, p])**2 + 0.00001)**0.5)
        self.Cmin = Expression(self.flowsheet().config.time, self.PH,
                               rule=rule_Cmin,
                               doc='Minimum capacitance rate')

        def rule_Cmax(blk, t, p):
            return 0.5 * (blk.Caph[t, p] + blk.Capc[t, p] +
                          ((blk.Caph[t, p] - blk.Capc[t, p])**2 + 0.00001)**0.5)
        self.Cmax = Expression(self.flowsheet().config.time, self.PH,
                               rule=rule_Cmax,
                               doc='Maximum capacitance rate')

        def rule_CR(blk, t, p):
            return blk.Cmin[t, p] / blk.Cmax[t, p]
        self.CR = Expression(self.flowsheet().config.time, self.PH,
                             rule=rule_CR,
                             doc='Capacitance ratio')

        # ----------------------------------------------------------------------
        # Number of Transfer units for sub heat exchanger
        def rule_NTU(blk, t, p):
            return blk.U[t, p] * blk.plate_area / blk.Cmin[t, p]
        self.NTU = Expression(self.flowsheet().config.time, self.PH,
                              rule=rule_NTU,
                              doc='Number of Transfer Units')

        # ----------------------------------------------------------------------
        # effectiveness of sub-heat exchangers
        def rule_Ecf(blk, t, p):
            if blk.P.value % 2 == 0:
                return (1 - exp(-blk.NTU[t, p] * (1 - blk.CR[t, p]))) / \
                    (1 - blk.CR[t, p] *
                     exp(-blk.NTU[t, p] * (1 - blk.CR[t, p])))
            elif blk.P.value % 2 == 1:
                return (1 - exp(-blk.NTU[t, p] * (1 + blk.CR[t, p]))) / (1 + blk.CR[t, p])

        self.Ecf = Expression(self.flowsheet().config.time, self.PH,
                              rule=rule_Ecf,
                              doc='Effectiveness for sub-HX')

        # ----------------------------------------------------------------------
        # Energy balance equations for hot fluid in sub-heat exhanger
        def rule_Ebh_eq(blk, t, p):
            return blk.Th_out[t, p] == blk.Th_in[t, p] -\
                blk.Ecf[t, p] * blk.Cmin[t, p] / blk.Caph[t, p] * \
                (blk.Th_in[t, p] - blk.Tc_in[t, p])

        self.Ebh_eq = Constraint(self.flowsheet().config.time, self.PH,
                                 rule=rule_Ebh_eq,
                                 doc='Hot fluid sub-heat exchanger energy balance')

        # Hot fluid exit temperature
        def rule_Tout_hot(blk, t):
            return blk.Th_out[t, blk.P.value] ==\
                blk.hot_side.properties_out[t].temperature

        self.Tout_hot_eq = Constraint(self.flowsheet().config.time,
                                      rule=rule_Tout_hot,
                                      doc='Hot fluid exit temperature')

        # Energy balance equations for cold fluid in sub-heat exhanger
        def rule_Ebc_eq(blk, t, p):
            return blk.Tc_out[t, p] == blk.Tc_in[t, p] + \
                blk.Ecf[t, p] * blk.Cmin[t, p] / blk.Capc[t, p] * \
                (blk.Th_in[t, p] - blk.Tc_in[t, p])

        self.Ebc_eq = Constraint(self.flowsheet().config.time, self.PH,
                                 rule=rule_Ebc_eq,
                                 doc='Cold fluid sub-heat exchanger energy balance')

        # Cold fluid exit temperature
        def rule_Tout_cold(blk, t):
            return blk.Tc_out[t, 1] ==\
                blk.cold_side.properties_out[t].temperature

        self.Tout_cold_eq = Constraint(self.flowsheet().config.time,
                                       rule=rule_Tout_cold,
                                       doc='Cold fluid exit temperature')

        # ----------------------------------------------------------------------
        # Energy balance boundary conditions
        def rule_hot_BCIN(blk, t):
            return blk.Th_in[t, 1] == \
                blk.hot_side.properties_in[t].temperature

        self.hot_BCIN = Constraint(self.flowsheet().config.time,
                                   rule=rule_hot_BCIN,
                                   doc='Hot fluid inlet boundary conditions')

        def rule_cold_BCIN(blk, t):
            return blk.Tc_in[t, blk.P.value] ==\
                blk.cold_side.properties_in[t].temperature
        self.cold_BCIN = Constraint(self.flowsheet().config.time,
                                    rule=rule_cold_BCIN,
                                    doc='Cold fluid inlet boundary conditions')

        Pset = [i for i in range(1, self.P.value)]

        def rule_hot_BC(blk, t, p):
            return blk.Th_out[t, p] == blk.Th_in[t, p + 1]
        self.hot_BC = Constraint(self.flowsheet().config.time, Pset,
                                 rule=rule_hot_BC,
                                 doc='Hot fluid boundary conditions: change of pass')

        def rule_cold_BC(blk, t, p):
            return blk.Tc_out[t, p + 1] == blk.Tc_in[t, p]
        self.cold_BC = Constraint(self.flowsheet().config.time, Pset,
                                  rule=rule_cold_BC,
                                  doc='Cold fluid boundary conditions: change of pass')

        # ----------------------------------------------------------------------
        # Energy transferred
        def rule_QH(blk, t):
            return blk.mh_in[t] * blk.cp_hot[t] *\
                (blk.hot_side.properties_in[t].temperature -
                 blk.hot_side.properties_out[t].temperature)
        self.QH = Expression(self.flowsheet().config.time, rule=rule_QH,
                             doc='Heat lost by hot fluid')

        def rule_QC(blk, t):
            return blk.mc_in[t] * blk.cp_cold[t] *\
                (blk.cold_side.properties_out[t].temperature -
                 blk.cold_side.properties_in[t].temperature)
        self.QC = Expression(self.flowsheet().config.time, rule=rule_QH,
                             doc='Heat gain by cold fluid')

    def initialize(blk, hotside_state_args={}, coldside_state_args={},
                   outlvl=5, solver='ipopt', optarg={'tol': 1e-6}):
        '''
        Initialisation routine for PHE unit (default solver ipopt)

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
                 * 0 = Use default idaes.init logger setting
                 * 1 = Maximum output
                 * 2 = Include solver output
                 * 3 = Return solver state for each step in subroutines
                 * 4 = Return solver state for each step in routine
                 * 5 = Final initialization status and exceptions
                 * 6 = No output
            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating which solver to use during
                     initialization (default = 'ipopt')

        Returns:
            None
        '''
        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag='unit')
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        hotside_state_args = {
            'flow_mol': value(blk.hot_inlet.flow_mol[0]),
            'temperature': value(blk.hot_inlet.temperature[0]),
            'pressure': value(blk.hot_inlet.pressure[0]),
            'mole_frac_comp':
            {'H2O': value(blk.hot_inlet.mole_frac_comp[0, 'H2O']),
             'CO2': value(blk.hot_inlet.mole_frac_comp[0, 'CO2']),
             'MEA': value(blk.hot_inlet.mole_frac_comp[0, 'MEA'])}}

        coldside_state_args = {
            'flow_mol': value(blk.cold_inlet.flow_mol[0]),
            'temperature': value(blk.cold_inlet.temperature[0]),
            'pressure': value(blk.cold_inlet.pressure[0]),
            'mole_frac_comp':
            {'H2O': value(blk.cold_inlet.mole_frac_comp[0, 'H2O']),
             'CO2': value(blk.cold_inlet.mole_frac_comp[0, 'CO2']),
             'MEA': value(blk.cold_inlet.mole_frac_comp[0, 'MEA'])}}

        # ---------------------------------------------------------------------
        # Initialize the INLET properties
        print('==============================================================')
        print('STEP 1: PROPERTY INITIALIZATION')
        print("======================INLET Properties========================")
        blk.hot_side.properties_in.initialize(
            state_args=hotside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)
        blk.cold_side.properties_in.initialize(
            state_args=coldside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)
        # Initialize the OUTLET properties
        print("======================OUTLET Properties=======================")
        blk.hot_side.properties_out.initialize(
            state_args=hotside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False)

        blk.cold_side.properties_out.initialize(
            state_args=coldside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False)
        # ----------------------------------------------------------------------
        print('==========================STEP 2==============================')
        print('STEP 2: solve PHE model')
        # solve step 2
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("STEP 2 Complete: {}.".format(idaeslog.condition(res)))
        print('==========================END=================================')
