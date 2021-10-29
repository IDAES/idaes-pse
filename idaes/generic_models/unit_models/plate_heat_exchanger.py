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
IDAES plate heat exchanger model (PHE) using effectiveness-NTU method

Notes:
Let P = number of passes(same for both fluids)
if P is even PHE is in counter-current mode
if P is odd PHE  is in co-current mode
Hot and cold fluids flow alternatively in the channels.
Divider plates may be used to partition the passes in separate sections.
The heat transfer area  of a single plate depends on the total heat transfer
area as specified by the manufacturer and the total number of active plates.

Detailed model equations can be found in the paper :
Akula, P., Eslick, J., Bhattacharyya, D. and Miller, D.C., 2019.
Modelling and Parameter Estimation of a Plate Heat Exchanger
as Part of a Solvent-Based Post-Combustion CO2 Capture System.
In Computer Aided Chemical Engineering (Vol. 47, pp. 47-52). Elsevier.

"""

# Import Pyomo libraries
from pyomo.environ import Param, RangeSet, Constraint, Expression,\
    value, Var, exp, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        MaterialBalanceType,
                        MomentumBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.constants import Constants as CONST
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util import get_solver
import idaes.logger as idaeslog

__author__ = "Paul Akula"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PlateHeatExchanger")
class PlateHeatExchangerData(UnitModelBlockData):
    """Plate Heat Exchanger(PHE) Unit Model."""

    CONFIG = UnitModelBlockData.CONFIG()

    # Configuration template for fluid specific  arguments
    _SideCONFIG = ConfigBlock()

    CONFIG.declare("passes", ConfigValue(
        default=4,
        domain=int,
        description="Number of passes",
        doc="""Number of passes of the fluids through the heat exchanger"""))

    CONFIG.declare("channel_list", ConfigValue(
        default=[12, 12, 12, 12],
        domain=list,
        description="Number of channels for each pass",
        doc="""Number of channels to be used in each pass where a channel
               is the space between two plates with a flowing fluid"""))

    CONFIG.declare("divider_plate_number", ConfigValue(
        default=0,
        domain=int,
        description="Number of divider plates in heat exchanger",
        doc="""Divider plates are used to create separate partitions in the
        unit. Each pass can be separated by a divider plate"""))

    CONFIG.declare("port_diameter", ConfigValue(
        default=0.2045,
        domain=float,
        description="Diameter of the ports on the plate [m]",
        doc="""Diameter of the ports on the plate for fluid entry/exit
               into a channel"""))

    CONFIG.declare("plate_thermal_cond", ConfigValue(
        default=16.2,
        domain=float,
        description="Thermal conductivity [W/m.K]",
        doc="""Thermal conductivity of the plate material [W/m.K]"""))

    CONFIG.declare("total_area", ConfigValue(
        default=114.3,
        domain=float,
        description="Total heat transfer area [m2]",
        doc="""Total heat transfer area as specifed by the manufacturer"""))

    CONFIG.declare("plate_thickness", ConfigValue(
        default=0.0006,
        domain=float,
        description="Plate thickness [m]",
        doc="""Plate thickness"""))

    CONFIG.declare("plate_vertical_dist", ConfigValue(
        default=1.897,
        domain=float,
        description="Vertical distance between centers of ports [m].",
        doc="""Vertical distance between centers of ports.(Top and bottom
        ports) (approximately equals to the plate length)"""))

    CONFIG.declare("plate_horizontal_dist", ConfigValue(
        default=0.409,
        domain=float,
        description="Horizontal distance between centers of ports [m].",
        doc="""Horizontal distance between centers of ports(Left and right
        ports)"""))

    CONFIG.declare("plate_pact_length", ConfigValue(
        default=0.381,
        domain=float,
        description="Compressed plate pact length [m].",
        doc="""Compressed plate pact length.
               Length between the Head and the Follower"""))

    CONFIG.declare("surface_enlargement_factor", ConfigValue(
        default=None,
        domain=float,
        description="Surface enlargement factor",
        doc="""Surface enlargement factor is the ratio of single plate area
               (obtained from the total area) to the projected plate area"""))

    CONFIG.declare("plate_gap", ConfigValue(
        default=None,
        domain=float,
        description="Mean channel spacing or gap bewteen two plates [m]",
        doc="""The plate gap is the distance between two adjacent plates that
               forms a flow channel """))

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
        super().build()

        # Consistency check for number of passes and channels in each pass
        for i in self.config.channel_list:
            if not isinstance(i, int):
                raise ConfigurationError("number of channels ({}) must be"
                                         " an integer".format(i))

        if (self.config.passes != len(self.config.channel_list)):
            raise ConfigurationError(
                f"The number of elements in the channel list ("
                f"{self.config.channel_list}) does not match the number of "
                f"passes ({self.config.passes}) given. Please provide as "
                f"integers, the number of channels of each pass")

        # ======================================================================
        # Build hot-side  Control Volume (Lean Solvent)
        self.hot_fluid = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.hot_side.property_package,
            "property_package_args":
                self.config.hot_side.property_package_args})

        self.hot_fluid.add_state_blocks(has_phase_equilibrium=False)

        self.hot_fluid.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_mass_transfer=False,
            has_phase_equilibrium=False,
            has_rate_reactions=False)

        self.hot_fluid.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=True)

        # Energy balance is based on the effectiveness Number of Transfer units
        # (E-NTU method) and inluded as performance equations. Hence the
        # control volume energy balances are not added.

        # =====================================================================
        # Build cold-side  Control Volume(Rich solvent)
        self.cold_fluid = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.cold_side.property_package,
            "property_package_args":
                self.config.cold_side.property_package_args})

        self.cold_fluid.add_state_blocks(has_phase_equilibrium=False)

        self.cold_fluid.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_mass_transfer=False,
            has_phase_equilibrium=False,
            has_rate_reactions=False)

        self.cold_fluid.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=True)

        # =====================================================================
        # Add Ports to control volumes
        # hot-side
        self.add_inlet_port(name="hot_inlet",
                            block=self.hot_fluid, doc='inlet Port')
        self.add_outlet_port(name="hot_outlet",
                             block=self.hot_fluid, doc='outlet Port')

        # cold-side
        self.add_inlet_port(name="cold_inlet",
                            block=self.cold_fluid, doc='inlet Port')
        self.add_outlet_port(name="cold_outlet",
                             block=self.cold_fluid, doc='outlet Port')
        # =====================================================================
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

        self.plate_thermal_cond = Param(
            mutable=True,
            initialize=self.config.plate_thermal_cond,
            units=pyunits.W / pyunits.m / pyunits.K,
            doc="Plate thermal conductivity")
        self.plate_thick = Param(mutable=True,
                                 initialize=self.config.plate_thickness,
                                 units=pyunits.m,
                                 doc="Plate thickness")

        self.port_dia = Param(mutable=True,
                              initialize=self.config.port_diameter,
                              units=pyunits.m,
                              doc=" Port diameter of plate ")
        self.Np = Param(self.PH,
                        units=None,
                        doc="Number of channels in each pass", mutable=True)
        # Number of channels in each pass
        for i in self.PH:
            self.Np[i].value = self.config.channel_list[i - 1]

        # ---------------------------------------------------------------------
        # Assign plate specifications

        # effective plate length & width
        _effective_plate_length = self.config.plate_vertical_dist*pyunits.m - \
            self.config.port_diameter*pyunits.m
        _effective_plate_width = self.config.plate_horizontal_dist*pyunits.m + \
            self.config.port_diameter*pyunits.m

        self.plate_length = Expression(expr=_effective_plate_length)
        self.plate_width = Expression(expr=_effective_plate_width)

        # Area of single plate
        _total_active_plate_number = 2 * sum(self.config.channel_list) - 1 -\
            self.config.divider_plate_number

        self.plate_area = Expression(expr=self.config.total_area*pyunits.m**2 /
                                     _total_active_plate_number,
                                     doc="Heat transfer area of single plate")

        # Plate gap
        if self.config.plate_gap is None:
            _total_plate_number = 2 * sum(self.config.channel_list) + 1 +\
                self.config.divider_plate_number
            _plate_pitch = self.config.plate_pact_length*pyunits.m / _total_plate_number

            _plate_gap = _plate_pitch - self.config.plate_thickness*pyunits.m
        else:
            _plate_gap = self.config.plate_gap*pyunits.m

        self.plate_gap = Expression(expr=_plate_gap)

        # Surface enlargement factor
        if self.config.surface_enlargement_factor is None:
            _projected_plate_area = (
                _effective_plate_length * _effective_plate_width)
            _surface_enlargement_factor = (
                self.plate_area / _projected_plate_area)
        else:
            _surface_enlargement_factor = \
                self.config.surface_enlargement_factor

        self.surface_enlargement_factor = Expression(
            expr=_surface_enlargement_factor)

        # Channel equivalent diameter
        self.channel_dia = Expression(expr=2 * self.plate_gap /
                                      _surface_enlargement_factor,
                                      doc=" Channel equivalent diameter")

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

        solvent_list = \
            self.config.hot_side.property_package.solvent_set

        def rule_mass_frac_solvent_hot(blk, t, j):
            return blk.hot_fluid.properties_in[t].mass_frac_phase_comp['Liq', j] /\
                sum(blk.hot_fluid.properties_in[t].mass_frac_phase_comp['Liq', k]
                    for k in solvent_list)

        self.mass_frac_solvent_hot = Expression(
            self.flowsheet().time,
            solvent_list,
            rule=rule_mass_frac_solvent_hot,
            doc='mass fraction of solvent on free solute basis')

        def rule_cp_mass_comp_hot(blk, t, j):
            return 0.5 * (
                blk.hot_fluid.properties_in[t].cp_mol_phase_comp['Liq', j] /
                blk.hot_fluid.properties_in[t].mw_comp[j] +
                blk.hot_fluid.properties_out[t].cp_mol_phase_comp['Liq', j] /
                blk.hot_fluid.properties_out[t].mw_comp[j])

        self.cp_mass_comp_hot = Expression(
            self.flowsheet().time,
            solvent_list,
            rule=rule_cp_mass_comp_hot,
            doc='Component mean specific heat capacity between inlet and '
            'outlet of hot-side temperature')

        def rule_cp_hot(blk, t):
            return sum(blk.cp_mass_comp_hot[t, j] *
                       blk.mass_frac_solvent_hot[t, j]
                       for j in solvent_list)
        self.cp_hot = Expression(self.flowsheet().time, rule=rule_cp_hot,
                                 doc='Hot-side mean specific heat capacity on'
                                     'free solute basis')

        def rule_mass_frac_solvent_cold(blk, t, j):
            return blk.cold_fluid.properties_in[t].mass_frac_phase_comp['Liq', j] /\
                sum(blk.cold_fluid.properties_in[t].mass_frac_phase_comp['Liq', k]
                    for k in solvent_list)

        self.mass_frac_solvent_cold = Expression(
            self.flowsheet().time,
            solvent_list,
            rule=rule_mass_frac_solvent_cold,
            doc='mass fraction of solvent on free solute basis')

        def rule_cp_mass_comp_cold(blk, t, j):
            return 0.5 * (
                blk.cold_fluid.properties_in[t].cp_mol_phase_comp['Liq', j] /
                blk.cold_fluid.properties_in[t].mw_comp[j] +
                blk.cold_fluid.properties_out[t].cp_mol_phase_comp['Liq', j] /
                blk.cold_fluid.properties_out[t].mw_comp[j])

        self.cp_mass_comp_cold = Expression(
            self.flowsheet().time,
            solvent_list,
            rule=rule_cp_mass_comp_cold,
            doc='Component mean specific heat capacity between inlet and '
            'outlet of hot-side temperature')

        def rule_cp_cold(blk, t):
            return sum(blk.cp_mass_comp_cold[t, j] *
                       blk.mass_frac_solvent_cold[t, j]
                       for j in solvent_list)
        self.cp_cold = Expression(self.flowsheet().time, rule=rule_cp_cold,
                                  doc='Cold-side mean specific heat capacity on'
                                      'free solute basis')

        # Model Variables
        self.Th_in = Var(self.flowsheet().time,
                         self.PH, initialize=393, units=pyunits.K,
                         doc="Hot Temperature IN of pass")
        self.Th_out = Var(self.flowsheet().time,
                          self.PH, initialize=325, units=pyunits.K,
                          doc="Hot Temperature OUT of pass")
        self.Tc_in = Var(self.flowsheet().time,
                         self.PH, initialize=320, units=pyunits.K,
                         doc="Cold Temperature IN of pass")
        self.Tc_out = Var(self.flowsheet().time,
                          self.PH, initialize=390, units=pyunits.K,
                          doc="Cold Temperature OUT of pass")

        # =====================================================================
        # PERFORMANCE EQUATIONS
        # mass flow rate in kg/s
        def rule_mh_in(blk, t):
            return blk.hot_fluid.properties_in[t].flow_mol *\
                blk.hot_fluid.properties_in[t].mw
        self.mh_in = Expression(self.flowsheet().time, rule=rule_mh_in,
                                doc='Hotside mass flow rate [kg/s]')

        def rule_mc_in(blk, t):
            return blk.cold_fluid.properties_in[t].flow_mol *\
                blk.cold_fluid.properties_in[t].mw
        self.mc_in = Expression(self.flowsheet().time, rule=rule_mc_in,
                                doc='Coldside mass flow rate [kg/s]')

        # ---------------------------------------------------------------------
        # port mass velocity[kg/m2.s]
        def rule_Gph(blk, t):
            return (4 * blk.mh_in[t] * 7) / (22 * blk.port_dia**2)

        self.Gph = Expression(self.flowsheet().time, rule=rule_Gph,
                              doc='Hotside port mass velocity[kg/m2.s]')

        def rule_Gpc(blk, t):
            return (4 * blk.mc_in[t] * 7) / (22 * blk.port_dia**2)

        self.Gpc = Expression(self.flowsheet().time, rule=rule_Gpc,
                              doc='Coldside port mass velocity[kg/m2.s]')

        # ---------------------------------------------------------------------
        # Reynold & Prandtl numbers
        def rule_Re_h(blk, t, p):
            return blk.mh_in[t] * blk.channel_dia /\
                (blk.Np[p] * blk.plate_width *
                 blk.plate_gap * blk.hot_fluid.properties_in[t].visc_d_phase["Liq"])
        self.Re_h = Expression(self.flowsheet().time, self.PH,
                               rule=rule_Re_h,
                               doc='Hotside Reynolds number')

        def rule_Re_c(blk, t, p):
            return blk.mc_in[t] * blk.channel_dia /\
                (blk.Np[p] * blk.plate_width *
                 blk.plate_gap * blk.cold_fluid.properties_in[t].visc_d_phase["Liq"])
        self.Re_c = Expression(self.flowsheet().time, self.PH,
                               rule=rule_Re_c,
                               doc='Coldside Reynolds number')

        def rule_Pr_h(blk, t):
            return blk.cp_hot[t] * blk.hot_fluid.properties_in[t].visc_d_phase["Liq"] /\
                blk.hot_fluid.properties_in[t].therm_cond_phase["Liq"]
        self.Pr_h = Expression(self.flowsheet().time,
                               rule=rule_Pr_h,
                               doc='Hotside Prandtl number')

        def rule_Pr_c(blk, t):
            return blk.cp_cold[t] * blk.cold_fluid.properties_in[t].visc_d_phase["Liq"] /\
                blk.cold_fluid.properties_in[t].therm_cond_phase["Liq"]
        self.Pr_c = Expression(self.flowsheet().time,
                               rule=rule_Pr_c,
                               doc='Coldside Prandtl number')
        # ---------------------------------------------------------------------
        # Film heat transfer coefficients

        def rule_hotside_transfer_coef(blk, t, p):
            return (blk.hot_fluid.properties_in[t].therm_cond_phase["Liq"] /
                    blk.channel_dia *
                    blk.param_a * blk.Re_h[t, p]**blk.param_b *
                    blk.Pr_h[t]**blk.param_c)

        self.h_hot = Expression(self.flowsheet().time,
                                self.PH, rule=rule_hotside_transfer_coef,
                                doc='Hotside heat transfer coefficient')

        def rule_coldside_transfer_coef(blk, t, p):
            return (blk.cold_fluid.properties_in[t].therm_cond_phase["Liq"] /
                    blk.channel_dia *
                    blk.param_a * blk.Re_c[t, p]**blk.param_b *
                    blk.Pr_c[t]**blk.param_c)

        self.h_cold = Expression(self.flowsheet().time,
                                 self.PH, rule=rule_coldside_transfer_coef,
                                 doc='Coldside heat transfer coefficient')

        # ---------------------------------------------------------------------
        # Friction factor calculation
        def rule_fric_h(blk, t):
            return 18.29 * blk.Re_h[t, 1]**(-0.652)
        self.fric_h = Expression(self.flowsheet().time,
                                 rule=rule_fric_h,
                                 doc='Hotside friction factor')

        def rule_fric_c(blk, t):
            return 1.441 * self.Re_c[t, 1]**(-0.206)
        self.fric_c = Expression(self.flowsheet().time,
                                 rule=rule_fric_c,
                                 doc='Coldside friction factor')

        # ---------------------------------------------------------------------
        # pressure drop calculation
        def rule_hotside_dP(blk, t):
            return (2 * blk.fric_h[t] * (blk.plate_length + blk.port_dia) *
                    blk.P * blk.Gph[t]**2) /\
                (blk.hot_fluid.properties_in[t].dens_mass *
                 blk.channel_dia) + 1.4 * blk.P * blk.Gph[t]**2 * 0.5 /\
                blk.hot_fluid.properties_in[t].dens_mass + \
                blk.hot_fluid.properties_in[t].dens_mass* \
                CONST.acceleration_gravity * (blk.plate_length + blk.port_dia)

        self.dP_h = Expression(self.flowsheet().time,
                               rule=rule_hotside_dP,
                               doc='Hotside pressure drop  [Pa]')

        def rule_coldside_dP(blk, t):
            return (2 * blk.fric_c[t] * (blk.plate_length + blk.port_dia) *
                    blk.P * blk.Gpc[t]**2) /\
                (blk.cold_fluid.properties_in[t].dens_mass * blk.channel_dia) +\
                1.4 * (blk.P * blk.Gpc[t]**2 * 0.5 /
                       blk.cold_fluid.properties_in[t].dens_mass) + \
                blk.cold_fluid.properties_in[t].dens_mass * \
                CONST.acceleration_gravity * (blk.plate_length + blk.port_dia)

        self.dP_c = Expression(self.flowsheet().time,
                               rule=rule_coldside_dP,
                               doc='Coldside pressure drop  [Pa]')

        def rule_eq_deltaP_hot(blk, t):
            return blk.hot_fluid.deltaP[t] == -blk.dP_h[t]
        self.eq_deltaP_hot = Constraint(self.flowsheet().time,
                                        rule=rule_eq_deltaP_hot)

        def rule_eq_deltaP_cold(blk, t):
            return blk.cold_fluid.deltaP[t] == -blk.dP_c[t]
        self.eq_deltaP_cold = Constraint(self.flowsheet().time,
                                         rule=rule_eq_deltaP_cold)

        # ---------------------------------------------------------------------
        # Overall heat transfer coefficients
        def rule_U(blk, t, p):
            return 1.0 /\
                (1.0 / blk.h_hot[t, p] +
                 blk.plate_gap / blk.plate_thermal_cond +
                 1.0 / blk.h_cold[t, p])
        self.U = Expression(self.flowsheet().time, self.PH,
                            rule=rule_U,
                            doc='Overall heat transfer coefficient')
        # ---------------------------------------------------------------------
        # capacitance of hot and cold fluid

        def rule_Caph(blk, t, p):
            return blk.mh_in[t] * blk.cp_hot[t] / blk.Np[p]
        self.Caph = Expression(self.flowsheet().time, self.PH,
                               rule=rule_Caph,
                               doc='Hotfluid capacitance rate')

        def rule_Capc(blk, t, p):
            return blk.mc_in[t] * blk.cp_cold[t] / blk.Np[p]
        self.Capc = Expression(self.flowsheet().time, self.PH,
                               rule=rule_Capc,
                               doc='Coldfluid capacitance rate')

        # ---------------------------------------------------------------------
        # min n max capacitance and capacitance ratio
        def rule_Cmin(blk, t, p):
            return 0.5 * (blk.Caph[t, p] + blk.Capc[t, p] -
                          ((blk.Caph[t, p] - blk.Capc[t, p])**2 +
                           0.00001)**0.5)
        self.Cmin = Expression(self.flowsheet().time, self.PH,
                               rule=rule_Cmin,
                               doc='Minimum capacitance rate')

        def rule_Cmax(blk, t, p):
            return 0.5 * (blk.Caph[t, p] + blk.Capc[t, p] +
                          ((blk.Caph[t, p] - blk.Capc[t, p])**2 +
                           0.00001)**0.5)
        self.Cmax = Expression(self.flowsheet().time, self.PH,
                               rule=rule_Cmax,
                               doc='Maximum capacitance rate')

        def rule_CR(blk, t, p):
            return blk.Cmin[t, p] / blk.Cmax[t, p]
        self.CR = Expression(self.flowsheet().time, self.PH,
                             rule=rule_CR,
                             doc='Capacitance ratio')

        # ---------------------------------------------------------------------
        # Number of Transfer units for sub heat exchanger
        def rule_NTU(blk, t, p):
            return blk.U[t, p] * blk.plate_area / blk.Cmin[t, p]
        self.NTU = Expression(self.flowsheet().time, self.PH,
                              rule=rule_NTU,
                              doc='Number of Transfer Units')

        # ---------------------------------------------------------------------
        # effectiveness of sub-heat exchangers
        def rule_Ecf(blk, t, p):
            if blk.P.value % 2 == 0:
                return (1 - exp(-blk.NTU[t, p] * (1 - blk.CR[t, p]))) / \
                    (1 - blk.CR[t, p] *
                     exp(-blk.NTU[t, p] * (1 - blk.CR[t, p])))
            elif blk.P.value % 2 == 1:
                return ((1 - exp(-blk.NTU[t, p] * (1 + blk.CR[t, p]))) /
                        (1 + blk.CR[t, p]))

        self.Ecf = Expression(self.flowsheet().time, self.PH,
                              rule=rule_Ecf,
                              doc='Effectiveness for sub-HX')

        # ---------------------------------------------------------------------
        # Energy balance equations for hot fluid in sub-heat exhanger
        def rule_Ebh_eq(blk, t, p):
            return blk.Th_out[t, p] == blk.Th_in[t, p] -\
                blk.Ecf[t, p] * blk.Cmin[t, p] / blk.Caph[t, p] * \
                (blk.Th_in[t, p] - blk.Tc_in[t, p])

        self.Ebh_eq = Constraint(
            self.flowsheet().time,
            self.PH,
            rule=rule_Ebh_eq,
            doc='Hot fluid sub-heat exchanger energy balance')

        # Hot fluid exit temperature
        def rule_Tout_hot(blk, t):
            return blk.Th_out[t, blk.P.value] ==\
                blk.hot_fluid.properties_out[t].temperature

        self.Tout_hot_eq = Constraint(self.flowsheet().time,
                                      rule=rule_Tout_hot,
                                      doc='Hot fluid exit temperature')

        # Energy balance equations for cold fluid in sub-heat exhanger
        def rule_Ebc_eq(blk, t, p):
            return blk.Tc_out[t, p] == blk.Tc_in[t, p] + \
                blk.Ecf[t, p] * blk.Cmin[t, p] / blk.Capc[t, p] * \
                (blk.Th_in[t, p] - blk.Tc_in[t, p])

        self.Ebc_eq = Constraint(
            self.flowsheet().time,
            self.PH,
            rule=rule_Ebc_eq,
            doc='Cold fluid sub-heat exchanger energy balance')

        # Cold fluid exit temperature
        def rule_Tout_cold(blk, t):
            return blk.Tc_out[t, 1] ==\
                blk.cold_fluid.properties_out[t].temperature

        self.Tout_cold_eq = Constraint(self.flowsheet().time,
                                       rule=rule_Tout_cold,
                                       doc='Cold fluid exit temperature')

        # ---------------------------------------------------------------------
        # Energy balance boundary conditions
        def rule_hot_BCIN(blk, t):
            return blk.Th_in[t, 1] == \
                blk.hot_fluid.properties_in[t].temperature

        self.hot_BCIN = Constraint(self.flowsheet().time,
                                   rule=rule_hot_BCIN,
                                   doc='Hot fluid inlet boundary conditions')

        def rule_cold_BCIN(blk, t):
            return blk.Tc_in[t, blk.P.value] ==\
                blk.cold_fluid.properties_in[t].temperature
        self.cold_BCIN = Constraint(self.flowsheet().time,
                                    rule=rule_cold_BCIN,
                                    doc='Cold fluid inlet boundary conditions')

        Pset = [i for i in range(1, self.P.value)]

        def rule_hot_BC(blk, t, p):
            return blk.Th_out[t, p] == blk.Th_in[t, p + 1]
        self.hot_BC = Constraint(
            self.flowsheet().time,
            Pset,
            rule=rule_hot_BC,
            doc='Hot fluid boundary conditions: change of pass')

        def rule_cold_BC(blk, t, p):
            return blk.Tc_out[t, p + 1] == blk.Tc_in[t, p]
        self.cold_BC = Constraint(
            self.flowsheet().time, Pset,
            rule=rule_cold_BC,
            doc='Cold fluid boundary conditions: change of pass')

        # ---------------------------------------------------------------------
        # Energy transferred
        def rule_QH(blk, t):
            return blk.mh_in[t] * blk.cp_hot[t] *\
                (blk.hot_fluid.properties_in[t].temperature -
                 blk.hot_fluid.properties_out[t].temperature)
        self.QH = Expression(self.flowsheet().time, rule=rule_QH,
                             doc='Heat lost by hot fluid')

        def rule_QC(blk, t):
            return blk.mc_in[t] * blk.cp_cold[t] *\
                (blk.cold_fluid.properties_out[t].temperature -
                 blk.cold_fluid.properties_in[t].temperature)
        self.QC = Expression(self.flowsheet().time, rule=rule_QH,
                             doc='Heat gain by cold fluid')

    def initialize(blk, hotside_state_args=None, coldside_state_args=None,
                   outlvl=idaeslog.NOTSET, solver=None, optarg=None):
        '''
        Initialisation routine for PHE unit (default solver ipopt)

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                           package(s) to provide an initial state for
                           initialization (see documentation of the specific
                           property package) (default = {}).
            outlvl : sets output level of initialization routine
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during
                     initialization (default = None)
        Returns:
            None
        '''
        # Set solver options
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag='unit')
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Create solver
        opt = get_solver(solver, optarg)

        hot_fluid_mole_frac_comp_dict = dict()
        cold_fluid_mole_frac_comp_dict = dict()

        for i in blk.config.hot_side.property_package.apparent_species_set:
            hot_fluid_mole_frac_comp_dict[i] =\
            value(blk.hot_inlet.mole_frac_comp[0, i])

        for i in blk.config.cold_side.property_package.apparent_species_set:
            cold_fluid_mole_frac_comp_dict[i] =\
            value(blk.cold_inlet.mole_frac_comp[0, i])

        hotside_state_args = {
            'flow_mol': value(blk.hot_inlet.flow_mol[0]),
            'temperature': value(blk.hot_inlet.temperature[0]),
            'pressure': value(blk.hot_inlet.pressure[0]),
            'mole_frac_comp': hot_fluid_mole_frac_comp_dict}

        coldside_state_args = {
            'flow_mol': value(blk.cold_inlet.flow_mol[0]),
            'temperature': value(blk.cold_inlet.temperature[0]),
            'pressure': value(blk.cold_inlet.pressure[0]),
            'mole_frac_comp': cold_fluid_mole_frac_comp_dict}

        # ---------------------------------------------------------------------
        # Initialize the Inlet properties
        init_log.info('Step 1: Property Initialization')
        init_log.info_high("Inlet Properties initialization")
        blk.hot_fluid.properties_in.initialize(
            state_args=hotside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)
        blk.cold_fluid.properties_in.initialize(
            state_args=coldside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)

        # Initialize the Outlet properties
        init_log.info_high("Outlet Properties initialization")
        blk.hot_fluid.properties_out.initialize(
            state_args=hotside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False)

        blk.cold_fluid.properties_out.initialize(
            state_args=coldside_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=False)
        # ----------------------------------------------------------------------
        init_log.info('Step 2: PHE Initialization')
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "Step 2 Complete: {}.".format(idaeslog.condition(res)))
        init_log.info('Initialization Completed')

    def _get_stream_table_contents(self, time_point=0):
        return create_stream_table_dataframe(
            {
                "Hot Inlet": self.hot_inlet,
                "Hot Outlet": self.hot_outlet,
                "Cold Inlet": self.cold_inlet,
                "Cold Outlet": self.cold_outlet,
            },
            time_point=time_point,
        )


if __name__ == '__main__':

    from pyomo.environ import ConcreteModel, value
    from idaes.core import FlowsheetBlock
    import matplotlib.pyplot as plt

    from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
    from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent \
        import configuration as aqueous_mea

    from idaes.core.util.model_statistics import (degrees_of_freedom,
                                                  number_variables,
                                                  number_total_constraints,
                                                  number_unused_variables)

    solver = get_solver()
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Set up property package
    m.fs.hotside_properties = GenericParameterBlock(default=aqueous_mea)
    m.fs.coldside_properties = GenericParameterBlock(default=aqueous_mea)

    # create instance of plate heat exchanger  on flowsheet
    m.fs.unit = PlateHeatExchanger(default={'passes': 4,
                                            'channel_list': [12, 12, 12, 12],
                                            'divider_plate_number': 2,
                                            "hot_side": {
                                                "property_package": m.fs.hotside_properties
                                            },
                                            "cold_side": {
                                                "property_package": m.fs.coldside_properties
                                            }})
    # hot fluid
    m.fs.unit.hot_inlet.flow_mol[0].fix(60.54879)
    m.fs.unit.hot_inlet.temperature[0].fix(392.23)
    m.fs.unit.hot_inlet.pressure[0].fix(202650)
    m.fs.unit.hot_inlet.mole_frac_comp[0, "CO2"].fix(0.0158)
    m.fs.unit.hot_inlet.mole_frac_comp[0, "H2O"].fix(0.8747)
    m.fs.unit.hot_inlet.mole_frac_comp[0, "MEA"].fix(0.1095)

    # cold fluid
    m.fs.unit.cold_inlet.flow_mol[0].fix(63.01910)
    m.fs.unit.cold_inlet.temperature[0].fix(326.36)
    m.fs.unit.cold_inlet.pressure[0].fix(202650)
    m.fs.unit.cold_inlet.mole_frac_comp[0, "CO2"].fix(0.0414)
    m.fs.unit.cold_inlet.mole_frac_comp[0, "H2O"].fix(0.8509)
    m.fs.unit.cold_inlet.mole_frac_comp[0, "MEA"].fix(0.1077)

    print('dof = {}'.format(degrees_of_freedom(m.fs.unit)))
    m.fs.unit.initialize()
    solver.solve(m.fs.unit)
    print('PC_out = {}'.format(value(m.fs.unit.cold_outlet.pressure[0])))
    print('PH_out = {}'.format(value(m.fs.unit.hot_outlet.pressure[0])))
    m.fs.unit.report()

    '''
    #molar flowrate
    hot_molar_flowrate =  [60.54879,102.07830,60.05750,58.86128,104.30124,105.49475,27.53904,
                           60.88287,60.88904,60.04379,77.69829,76.66811,59.68288]
    cold_molar_flowrate = [63.01910,104.99350,62.53341,60.40000,106.18207,107.38606, 28.19399,
                           63.60044,63.16421,61.14541,81.36657,79.89472,60.39896]

    # plant inlet temperature
    hot_Temp_IN  = [392.23,389.57,393.78,382.42,376.32,392.69,389.69,392.10,392.36,392.30,391.19,390.95,392.26]
    cold_Temp_IN = [326.36,332.26,329.12,318.82,319.58,330.54,321.42,327.72,327.47,326.72,328.59,325.44,326.11]

    # plant exit temperature
    NCCC_hot_Temp_OUT=[330.42,336.70,331.44,323.41,324.57,334.63,324.83,331.25,331.08,330.43,332.96,329.96,329.75]
    NCCC_cold_Temp_OUT= [384.9111111,383.2111111,383.6944444,376.1722222,370.5,384.9555556,382.0388889,
                        384.6388889,384.8055556,384.5888889,384.3277778,383.5333333,384.5277778]

    # hot-side (lean solvent) inlet mole fraction
    xh_H2O=[0.8747,0.8569,0.8739,0.8505,0.8573,0.8808,0.8582,0.8757,0.8758,0.8702,0.8646,0.8590,0.8676]
    xh_MEA=[0.1095,0.1148,0.1138,0.1110,0.1020,0.1033,0.1144,0.1070,0.1071,0.1115,0.1106,0.1152,0.1134]
    xh_CO2=[0.0158,0.0284,0.0123,0.0385,0.0407,0.0160,0.0274,0.0172,0.0171,0.0183,0.0248,0.0258,0.0190]

    # cold-side(rich solvent) inlet mole fraction
    xc_H2O=[0.8509,0.8426,0.8546,0.8378,0.8501,0.8676,0.8315,0.8554,0.8586,0.8490,0.8468,0.8408,0.8455]
    xc_MEA=[0.1077,0.1137,0.1123,0.1104,0.1019,0.1038,0.1143,0.1050,0.1055,0.1110,0.1079,0.1127,0.1141]
    xc_CO2=[0.0414,0.0438,.0331,0.0518,0.0480,0.0286,0.0541,0.0397,0.0360,0.0400,0.0453,0.0465,0.0404]

    # create list to save results
    PHE_THOUT=[]
    PHE_TCOUT=[]

    # fix the inputs, solve and save results
    for i in range(1):
        # hot fluid
        m.fs.unit.hot_inlet.flow_mol[0].fix(hot_molar_flowrate[i])
        m.fs.unit.hot_inlet.temperature[0].fix(hot_Temp_IN[i])
        m.fs.unit.hot_inlet.mole_frac_comp[0,"CO2"].fix(xh_CO2[i])
        m.fs.unit.hot_inlet.mole_frac_comp[0,"H2O"].fix(xh_H2O[i])
        m.fs.unit.hot_inlet.mole_frac_comp[0,"MEA"].fix(xh_MEA[i])

        #cold fluid
        m.fs.unit.cold_inlet.flow_mol[0].fix(cold_molar_flowrate[i])
        m.fs.unit.cold_inlet.temperature[0].fix(cold_Temp_IN[i])
        m.fs.unit.cold_inlet.mole_frac_comp[0,"CO2"].fix(xc_CO2[i])
        m.fs.unit.cold_inlet.mole_frac_comp[0,"H2O"].fix(xc_H2O[i])
        m.fs.unit.cold_inlet.mole_frac_comp[0,"MEA"].fix(xc_MEA[i])

        solver.solve(m.fs.unit,tee=False)
        PHE_THOUT.append(value(m.fs.unit.hot_outlet.temperature[0]))
        PHE_TCOUT.append(value(m.fs.unit.cold_outlet.temperature[0]))

    m.fs.unit.report()



    fontsize = 16
    labelsize = 16
    markersize=12

    x= [i for i in range(1,len(NCCC_cold_Temp_OUT)+1)]

    plt.figure(figsize=(8,6))

    plt.plot(x,NCCC_cold_Temp_OUT,
            color='g',
            linestyle='',
            label='Data: Rich solvent',
            mfc="None",
            marker='s',
            markersize=markersize)
    plt.plot(x,PHE_TCOUT,
            color='g',
            linestyle='',
            label='Model: Rich solvent',
            mfc="None",
            marker='x',
            markersize=markersize)
    plt.plot(x,NCCC_hot_Temp_OUT,
            color='r',
            linestyle='',
            label='Data: Lean solvent',
            mfc="None",
            marker='o',
            markersize=markersize)
    plt.plot(x,PHE_THOUT,
            color='r',
            linestyle='',
            label='Model: Lean solvent',
            mfc="None",
            marker='+',
            markersize=markersize)

    plt.ylim(278,400)
    plt.ylabel('Temperature (K)',fontsize=fontsize,fontweight='bold')
    plt.xlabel('NCCC Case No.',fontsize=fontsize,fontweight='bold')
    plt.legend(loc='lower right',fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    plt.tight_layout()
    plt.show()
    '''
    xgy = 1
