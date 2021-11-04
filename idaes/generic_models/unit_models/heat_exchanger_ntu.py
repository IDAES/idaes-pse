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

Assumptions:
    * No phase equilibrium or reactions occur within unit

Reference
Detailed model equations can be found in the paper :
[1] Akula, P., Eslick, J., Bhattacharyya, D. and Miller, D.C., 2019.
Modelling and Parameter Estimation of a Plate Heat Exchanger
as Part of a Solvent-Based Post-Combustion CO2 Capture System.
In Computer Aided Chemical Engineering (Vol. 47, pp. 47-52). Elsevier.

"""

# Import Pyomo libraries
from pyomo.environ import Constraint, Reference, units as pyunits
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
from idaes.core.util.constants import Constants
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.core.util.math import smooth_min, smooth_max
from idaes.core.util import get_solver
import idaes.logger as idaeslog

__author__ = "Paul Akula, Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PlateHeatExchanger")
class PlateHeatExchangerData(UnitModelBlockData):
    """Plate Heat Exchanger(PHE) Unit Model."""

    CONFIG = UnitModelBlockData.CONFIG()

    # Configuration template for fluid specific  arguments
    _SideCONFIG = ConfigBlock()

    # CONFIG.declare("passes", ConfigValue(
    #     default=4,
    #     domain=int,
    #     description="Number of passes",
    #     doc="""Number of passes of the fluids through the heat exchanger"""))

    # CONFIG.declare("channel_list", ConfigValue(
    #     default=[12, 12, 12, 12],
    #     domain=list,
    #     description="Number of channels for each pass",
    #     doc="""Number of channels to be used in each pass where a channel
    #            is the space between two plates with a flowing fluid"""))

    # CONFIG.declare("divider_plate_number", ConfigValue(
    #     default=0,
    #     domain=int,
    #     description="Number of divider plates in heat exchanger",
    #     doc="""Divider plates are used to create separate partitions in the
    #     unit. Each pass can be separated by a divider plate"""))

    # CONFIG.declare("port_diameter", ConfigValue(
    #     default=(0.2045, pyunits.m),
    #     description="Diameter of the ports on the plate [m]",
    #     doc="""Diameter of the ports on the plate for fluid entry/exit
    #            into a channel"""))

    # CONFIG.declare("plate_thermal_cond", ConfigValue(
    #     default=(16.2, pyunits.W / pyunits.m / pyunits.K),
    #     description="Thermal conductivity [W/m.K]",
    #     doc="""Thermal conductivity of the plate material [W/m.K]"""))

    # CONFIG.declare("total_area", ConfigValue(
    #     default=(114.3, pyunits.m**2),
    #     description="Total heat transfer area [m2]",
    #     doc="""Total heat transfer area as specifed by the manufacturer"""))

    # CONFIG.declare("plate_thickness", ConfigValue(
    #     default=(0.0006, pyunits.m),
    #     description="Plate thickness [m]",
    #     doc="""Plate thickness"""))

    # CONFIG.declare("plate_vertical_dist", ConfigValue(
    #     default=(1.897, pyunits.m),
    #     description="Vertical distance between centers of ports [m].",
    #     doc="""Vertical distance between centers of ports.(Top and bottom
    #     ports) (approximately equals to the plate length)"""))

    # CONFIG.declare("plate_horizontal_dist", ConfigValue(
    #     default=(0.409, pyunits.m),
    #     description="Horizontal distance between centers of ports [m].",
    #     doc="""Horizontal distance between centers of ports(Left and right
    #     ports)"""))

    # CONFIG.declare("plate_pact_length", ConfigValue(
    #     default=(0.381, pyunits.m),
    #     description="Compressed plate pact length [m].",
    #     doc="""Compressed plate pact length.
    #            Length between the Head and the Follower"""))

    # CONFIG.declare("surface_enlargement_factor", ConfigValue(
    #     default=None,
    #     description="Surface enlargement factor",
    #     doc="""Surface enlargement factor is the ratio of single plate area
    #            (obtained from the total area) to the projected plate area"""))

    # CONFIG.declare("plate_gap", ConfigValue(
    #     default=None,
    #     description="Mean channel spacing or gap bewteen two plates [m]",
    #     doc="""The plate gap is the distance between two adjacent plates that
    #            forms a flow channel """))

    _SideCONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.useDefault,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    _SideCONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.useDefault,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    _SideCONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))

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

        # # Consistency check for number of passes and channels in each pass
        # for i in self.config.channel_list:
        #     if not isinstance(i, int):
        #         raise ConfigurationError("number of channels ({}) must be"
        #                                  " an integer".format(i))

        # if (self.config.passes != len(self.config.channel_list)):
        #     raise ConfigurationError(
        #         f"The number of elements in the channel list ("
        #         f"{self.config.channel_list}) does not match the number of "
        #         f"passes ({self.config.passes}) given. Please provide as "
        #         f"integers, the number of channels of each pass")

        # ---------------------------------------------------------------------
        # Build hot-side  Control Volume (Lean Solvent)
        self.hot_side = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.hot_side.property_package,
            "property_package_args":
                self.config.hot_side.property_package_args})

        # TODO : Add support for phase equilibrium?
        self.hot_side.add_state_blocks(has_phase_equilibrium=False)

        self.hot_side.add_material_balances(
            balance_type=self.config.hot_side.material_balance_type,
            has_phase_equilibrium=False)

        self.hot_side.add_energy_balances(
            balance_type=self.config.hot_side.energy_balance_type,
            has_heat_transfer=True)

        # TODO : Make pressure drop optional?
        self.hot_side.add_momentum_balances(
            balance_type=self.config.hot_side.momentum_balance_type,
            has_pressure_change=True)

        # ---------------------------------------------------------------------
        # Build cold-side  Control Volume(Rich solvent)
        self.cold_side = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.cold_side.property_package,
            "property_package_args":
                self.config.cold_side.property_package_args})

        self.cold_side.add_state_blocks(has_phase_equilibrium=False)

        self.cold_side.add_material_balances(
            balance_type=self.config.cold_side.material_balance_type,
            has_phase_equilibrium=False)

        self.cold_side.add_energy_balances(
            balance_type=self.config.cold_side.energy_balance_type,
            has_heat_transfer=True)

        # TODO : Make pressure drop optional?
        self.cold_side.add_momentum_balances(
            balance_type=self.config.cold_side.momentum_balance_type,
            has_pressure_change=True)

        # ---------------------------------------------------------------------
        # Add Ports to control volumes
        self.add_inlet_port(name="hot_inlet",
                            block=self.hot_side,
                            doc='Hot side inlet port')
        self.add_outlet_port(name="hot_outlet",
                             block=self.hot_side,
                             doc='Hot side outlet port')

        self.add_inlet_port(name="cold_inlet",
                            block=self.cold_side,
                            doc='Cold side inlet port')
        self.add_outlet_port(name="cold_outlet",
                             block=self.cold_side,
                             doc='Cold side outlet port')

        # ---------------------------------------------------------------------
        # Add unit level References
        # Set references to balance terms at unit level
        self.heat_duty = Reference(self.cold_side.heat[:])

        # ---------------------------------------------------------------------
        # Add performance equations
        hunits = self.config.hot_side.property_package.get_metadata(
            ).get_derived_units

        # Overall energy balance
        def rule_energy_balance(blk, t):
            return blk.hot_side.heat[t] == pyunits.convert(
                blk.cold_side.heat[t], to_units=hunits("power"))
        self.energy_balance_constraint = Constraint(
            self.flowsheet().time, rule=rule_energy_balance)


    # def _make_params(self):
    #     self.pass_num = Param(initialize=self.config.passes,
    #                           units=pyunits.dimensionless,
    #                           doc="Total number of passes for hot or cold fluid")

    #     # Design for equal number of passes for both fluids
    #     self.PH = RangeSet(self.pass_num,
    #                        doc="Set of hot/cold fluid passes")

    #     self.plate_thermal_cond = Param(
    #         mutable=True,
    #         units=pyunits.W / pyunits.m / pyunits.K,
    #         doc="Plate thermal conductivity")
    #     set_param_value(self, "plate_thermal_cond",
    #                     pyunits.W / pyunits.m / pyunits.K)

    #     self.plate_thickness = Param(mutable=True,
    #                                  units=pyunits.m,
    #                                  doc="Plate thickness")
    #     set_param_value(self, "plate_thickness", pyunits.m)

    #     self.port_diameter = Param(mutable=True,
    #                                units=pyunits.m,
    #                                doc=" Port diameter of plate ")
    #     set_param_value(self, "port_diameter", pyunits.m)

    #     self.divider_plate_number = Param(mutable=True,
    #                                       units=pyunits.dimensionless,
    #                                       doc="Number of divider plates in heat exchanger ")
    #     set_param_value(self, "divider_plate_number", pyunits.dimensionless)

    #     self.plate_pact_length = Param(mutable=True,
    #                                    units=pyunits.m,
    #                                    doc="Compressed plate pact length")
    #     set_param_value(self, "plate_pact_length", pyunits.m)

    #     self.Np = Param(self.PH,
    #                     units=pyunits.dimensionless,
    #                     doc="Number of channels in each pass", mutable=True)
    #     # Number of channels in each pass
    #     for i in self.PH:
    #         self.Np[i].value = self.config.channel_list[i - 1]

    #     # ---------------------------------------------------------------------
    #     # Assign plate specifications

    #     self.plate_vertical_dist = Param(
    #         mutable=True,
    #         units=pyunits.m,
    #         doc="Vertical distance between centers of ports")
    #     set_param_value(self, "plate_vertical_dist", pyunits.m)

    #     self.plate_horizontal_dist = Param(
    #         mutable=True,
    #         units=pyunits.m,
    #         doc="Horizontal distance between centers of ports")
    #     set_param_value(self, "plate_horizontal_dist", pyunits.m)

    #     self.plate_length = Expression(expr=self.plate_vertical_dist -
    #                                    self.port_diameter)
    #     self.plate_width = Expression(expr=self.plate_horizontal_dist +
    #                                   self.port_diameter)
    #     # heat transer area
    #     self.total_area = Var(doc='Total heat transfer area',
    #                           units=pyunits.m**2)
    #     set_param_value(self, "total_area", pyunits.m**2)
    #     self.total_area.fix()

    #     # total number of active plates
    #     _total_active_plate_number = 2 * sum(self.config.channel_list) - 1 -\
    #         self.divider_plate_number

    #     self.total_active_plate_number = Param(initialize=_total_active_plate_number,
    #                                            units=pyunits.dimensionless,
    #                                            doc="Total number of active plates")

    #     # Area of single plate
    #     self.plate_area = Var(initialize=1,
    #                           units=pyunits.m**2,
    #                           doc='Heat transfer area of single plate')

    #     # Equation for design of new plate_area/total_area/No.of active plates
    #     self.plate_area_eq = Constraint(expr=self.total_area ==
    #                                     self.plate_area * self.total_active_plate_number)

    #     # Plate gap
    #     if self.config.plate_gap is None:
    #         _total_plate_number = 2 * sum(self.config.channel_list) + 1 +\
    #             self.divider_plate_number
    #         _plate_pitch = self.plate_pact_length / _total_plate_number

    #         _plate_gap = _plate_pitch - self.plate_thickness
    #         self.plate_gap = Expression(expr=_plate_gap)
    #     else:
    #         self.plate_gap = Param(mutable=True,
    #                                units=pyunits.m,
    #                                doc="Gap bewteen two plates ")
    #         set_param_value(self, "plate_gap", pyunits.m)

    #     # Surface enlargement factor
    #     # projected plate area = plate length *plate width
    #     if self.config.surface_enlargement_factor is None:
    #         self.surface_enlargement_factor = Expression(
    #             expr=self.plate_area / (self.plate_length * self.plate_width))
    #     else:
    #         self.surface_enlargement_factor = Param(mutable=True,
    #                                                 units=pyunits.dimensionless,
    #                                                 doc="Surface enlargement factor ")
    #         set_param_value(self, "surface_enlargement_factor",
    #                         pyunits.dimensionless)

    #     # Channel equivalent diameter
    #     self.channel_diameter = Expression(expr=2 * self.plate_gap /
    #                                        self.surface_enlargement_factor,
    #                                        doc=" Channel equivalent diameter")

    #     # heat transfer parameters
    #     self.nusselt_param_a = Var(initialize=0.3,
    #                                bounds=(0.2, 0.4),
    #                                units=pyunits.dimensionless,
    #                                doc='Nusselt parameter')
    #     self.nusselt_param_b = Var(initialize=0.663,
    #                                bounds=(0.3, 0.7),
    #                                units=pyunits.dimensionless,
    #                                doc='Nusselt parameter')
    #     self.nusselt_param_c = Var(initialize=1 / 3.0,
    #                                bounds=(1e-5, 2),
    #                                units=pyunits.dimensionless,
    #                                doc='Nusselt parameter')

    #     self.nusselt_param_a.fix(0.4)
    #     self.nusselt_param_b.fix(0.663)
    #     self.nusselt_param_c.fix(0.333)

    #     # friction factor parameters
    #     self.frict_factor_param_a = Var(initialize=0.0,
    #                                     units=pyunits.dimensionless,
    #                                     doc='Friction factor parameter')
    #     self.frict_factor_param_b = Var(initialize=18.29,
    #                                     units=pyunits.dimensionless,
    #                                     doc='Friction factor parameter')
    #     self.frict_factor_param_c = Var(initialize=-0.652,
    #                                     units=pyunits.dimensionless,
    #                                     doc='Friction factor parameter')

    #     self.frict_factor_param_a.fix(0.0)
    #     self.frict_factor_param_b.fix(18.29)
    #     self.frict_factor_param_c.fix(-0.652)

    #     # espilon parameter for smooth Cmin and Cmax functions
    #     self.epsilon_param = Param(initialize= 1e-2,
    #                                units=pyunits.W/pyunits.K,
    #                                doc="Epsilon parameter for smooth Cmin and Cmax")

    # def _make_performance_method(self):

    #     solvent_list = \
    #         self.config.hot_side.property_package.solvent_set

    #     def rule_mass_frac_solvent_hot(blk, t, j):
    #         return blk.hot_fluid.properties_in[t].mass_frac_phase_comp['Liq', j] /\
    #             sum(blk.hot_fluid.properties_in[t].mass_frac_phase_comp['Liq', k]
    #                 for k in solvent_list)

    #     self.mass_frac_solvent_hot = Expression(
    #         self.flowsheet().time,
    #         solvent_list,
    #         rule=rule_mass_frac_solvent_hot,
    #         doc='Mass fraction of hot solvent on free solute basis')

    #     def rule_cp_mass_comp_hot(blk, t, j):
    #         return 0.5 * (
    #             blk.hot_fluid.properties_in[t].cp_mol_phase_comp['Liq', j] /
    #             blk.hot_fluid.properties_in[t].mw_comp[j] +
    #             blk.hot_fluid.properties_out[t].cp_mol_phase_comp['Liq', j] /
    #             blk.hot_fluid.properties_out[t].mw_comp[j])

    #     self.cp_mass_comp_hot = Expression(
    #         self.flowsheet().time,
    #         solvent_list,
    #         rule=rule_cp_mass_comp_hot,
    #         doc='Component mean specific heat capacity between inlet and '
    #         'outlet of hot-side temperature')

    #     def rule_cp_hot(blk, t):
    #         return sum(blk.cp_mass_comp_hot[t, j] *
    #                    blk.mass_frac_solvent_hot[t, j]
    #                    for j in solvent_list)
    #     self.cp_hot = Expression(self.flowsheet().time, rule=rule_cp_hot,
    #                              doc='Hot-side mean specific heat capacity on'
    #                                  'free solute basis')

    #     def rule_mass_frac_solvent_cold(blk, t, j):
    #         return blk.cold_fluid.properties_in[t].mass_frac_phase_comp['Liq', j] /\
    #             sum(blk.cold_fluid.properties_in[t].mass_frac_phase_comp['Liq', k]
    #                 for k in solvent_list)

    #     self.mass_frac_solvent_cold = Expression(
    #         self.flowsheet().time,
    #         solvent_list,
    #         rule=rule_mass_frac_solvent_cold,
    #         doc='Mass fraction of cold solvent on free solute basis')

    #     def rule_cp_mass_comp_cold(blk, t, j):
    #         return 0.5 * (
    #             blk.cold_fluid.properties_in[t].cp_mol_phase_comp['Liq', j] /
    #             blk.cold_fluid.properties_in[t].mw_comp[j] +
    #             blk.cold_fluid.properties_out[t].cp_mol_phase_comp['Liq', j] /
    #             blk.cold_fluid.properties_out[t].mw_comp[j])

    #     self.cp_mass_comp_cold = Expression(
    #         self.flowsheet().time,
    #         solvent_list,
    #         rule=rule_cp_mass_comp_cold,
    #         doc='Component mean specific heat capacity between inlet and '
    #         'outlet of hot-side temperature')

    #     def rule_cp_cold(blk, t):
    #         return sum(blk.cp_mass_comp_cold[t, j] *
    #                    blk.mass_frac_solvent_cold[t, j]
    #                    for j in solvent_list)
    #     self.cp_cold = Expression(self.flowsheet().time, rule=rule_cp_cold,
    #                               doc='Cold-side mean specific heat capacity on'
    #                                   'free solute basis')

    #     # Model Variables
    #     self.Th_in = Var(self.flowsheet().time,
    #                      self.PH, initialize=393, units=pyunits.K,
    #                      doc="Hot Temperature IN of pass")
    #     self.Th_out = Var(self.flowsheet().time,
    #                       self.PH, initialize=325, units=pyunits.K,
    #                       doc="Hot Temperature OUT of pass")
    #     self.Tc_in = Var(self.flowsheet().time,
    #                      self.PH, initialize=320, units=pyunits.K,
    #                      doc="Cold Temperature IN of pass")
    #     self.Tc_out = Var(self.flowsheet().time,
    #                       self.PH, initialize=390, units=pyunits.K,
    #                       doc="Cold Temperature OUT of pass")

    #     # =====================================================================
    #     # PERFORMANCE EQUATIONS

    #     # G:mass flow velocity, p :port , h:hotfluid , c :coldfluid
    #     # Gph is hotfluid port mass velocity[kg/m2.s]
    #     def rule_Gph(blk, t):
    #         return (4 * blk.hot_fluid.properties_in[t].flow_mass * 7) /\
    #                (22 * blk.port_diameter**2)

    #     Gph = self.port_mass_velocity_h = \
    #         Expression(self.flowsheet().time,
    #                    rule=rule_Gph,
    #                    doc='Hotside port mass velocity[kg/m2.s]')

    #     # Gpc is coldfluid port mass velocity[kg/m2.s]
    #     def rule_Gpc(blk, t):
    #         return (4 * blk.cold_fluid.properties_in[t].flow_mass * 7) /\
    #                (22 * blk.port_diameter**2)

    #     Gpc = self.port_mass_velocity_c =\
    #         Expression(self.flowsheet().time,
    #                    rule=rule_Gpc,
    #                    doc='Coldside port mass velocity[kg/m2.s]')

    #     # ---------------------------------------------------------------------
    #     # Gch is hotfluid channel mass velocity[kg/m2.s]
    #     def rule_Gch(blk, t):
    #         return blk.hot_fluid.properties_in[t].flow_mass / \
    #             (blk.Np[1] * blk.plate_width * blk.plate_gap)

    #     Gch = self.channel_mass_velocity_h = \
    #         Expression(self.flowsheet().time,
    #                    rule=rule_Gch,
    #                    doc='Hotside channel mass velocity[kg/m2.s]')

    #     # Gcc is coldfluid channel mass velocity[kg/m2.s]
    #     def rule_Gcc(blk, t):
    #         return blk.cold_fluid.properties_in[t].flow_mass / \
    #             (blk.Np[1] * blk.plate_width * blk.plate_gap)

    #     Gcc = self.channel_mass_velocity_c =\
    #         Expression(self.flowsheet().time,
    #                    rule=rule_Gcc,
    #                    doc='Coldside channel mass velocity[kg/m2.s]')

    #     # ---------------------------------------------------------------------
    #     # Reynold & Prandtl numbers
    #     def rule_Re_h(blk, t, p):
    #         return blk.hot_fluid.properties_in[t].flow_mass * blk.channel_diameter /\
    #             (blk.Np[p] * blk.plate_width *
    #              blk.plate_gap * blk.hot_fluid.properties_in[t].visc_d_phase["Liq"])


    #     Re_h = self.Reynold_number_h = \
    #         Expression(self.flowsheet().time,
    #                    self.PH,
    #                    rule=rule_Re_h,
    #                    doc='Hotside Reynolds number')

    #     def rule_Re_c(blk, t, p):
    #         return blk.cold_fluid.properties_in[t].flow_mass * blk.channel_diameter /\
    #             (blk.Np[p] * blk.plate_width *
    #              blk.plate_gap * blk.cold_fluid.properties_in[t].visc_d_phase["Liq"])


    #     Re_c = self.Reynold_number_c = \
    #         Expression(self.flowsheet().time,
    #                    self.PH,
    #                    rule=rule_Re_c,
    #                    doc='Coldside Reynolds number')

    #     def rule_Pr_h(blk, t):
    #         return blk.cp_hot[t] * blk.hot_fluid.properties_in[t].visc_d_phase["Liq"] /\
    #             blk.hot_fluid.properties_in[t].therm_cond_phase["Liq"]

    #     Pr_h = self.Prandtl_number_h = \
    #         Expression(self.flowsheet().time,
    #                    rule=rule_Pr_h,
    #                    doc='Hotside Prandtl number')

    #     def rule_Pr_c(blk, t):
    #         return blk.cp_cold[t] * blk.cold_fluid.properties_in[t].visc_d_phase["Liq"] /\
    #             blk.cold_fluid.properties_in[t].therm_cond_phase["Liq"]

    #     Pr_c = self.Prandtl_number_c = \
    #         Expression(self.flowsheet().time,
    #                    rule=rule_Pr_c,
    #                    doc='Coldside Prandtl number')
    #     # ---------------------------------------------------------------------
    #     # Film heat transfer coefficients

    #     def rule_hotside_transfer_coef(blk, t, p):
    #         return (blk.hot_fluid.properties_in[t].therm_cond_phase["Liq"] /
    #                 blk.channel_diameter *
    #                 blk.nusselt_param_a * Re_h[t, p]**blk.nusselt_param_b *
    #                 Pr_h[t]**blk.nusselt_param_c)

    #     h_hot = self.heat_trans_coef_hot =\
    #         Expression(self.flowsheet().time,
    #                    self.PH,
    #                    rule=rule_hotside_transfer_coef,
    #                    doc='Hotside heat transfer coefficient')

    #     def rule_coldside_transfer_coef(blk, t, p):
    #         return (blk.cold_fluid.properties_in[t].therm_cond_phase["Liq"] /
    #                 blk.channel_diameter *
    #                 blk.nusselt_param_a * Re_c[t, p]**blk.nusselt_param_b *
    #                 Pr_c[t]**blk.nusselt_param_c)

    #     h_cold = self.heat_trans_coef_cold =\
    #         Expression(self.flowsheet().time,
    #                    self.PH,
    #                    rule=rule_coldside_transfer_coef,
    #                    doc='Coldside heat transfer coefficient')

    #     # ---------------------------------------------------------------------
    #     # Friction factor calculation
    #     def rule_fric_h(blk, t):
    #         return blk.frict_factor_param_a + blk.frict_factor_param_b *\
    #             Re_h[t, 1]**(blk.frict_factor_param_c)

    #     self.fric_h = Expression(self.flowsheet().time,
    #                              rule=rule_fric_h,
    #                              doc='Hotside friction factor')

    #     def rule_fric_c(blk, t):
    #         return blk.frict_factor_param_a + blk.frict_factor_param_b *\
    #             Re_c[t, 1]**(blk.frict_factor_param_c)

    #     self.fric_c = Expression(self.flowsheet().time,
    #                              rule=rule_fric_c,
    #                              doc='Coldside friction factor')

    #     # ---------------------------------------------------------------------
    #     # pressure drop calculation
    #     def rule_hotside_dP(blk, t):
    #         return (2 * blk.fric_h[t] * (blk.plate_length + blk.port_diameter) *
    #                 blk.pass_num * Gch[t]**2) /\
    #             (blk.hot_fluid.properties_in[t].dens_mass *
    #              blk.channel_diameter) + 0.7 * blk.pass_num * Gph[t]**2 /\
    #             blk.hot_fluid.properties_in[t].dens_mass + \
    #             blk.hot_fluid.properties_in[t].dens_mass * \
    #             CONST.acceleration_gravity * \
    #             (blk.plate_length + blk.port_diameter)

    #     self.deltaP_h = Expression(self.flowsheet().time,
    #                                rule=rule_hotside_dP,
    #                                doc='Hotside pressure drop  [Pa]')

    #     def rule_coldside_dP(blk, t):
    #         return (2 * blk.fric_c[t] * (blk.plate_length + blk.port_diameter) *
    #                 blk.pass_num * Gcc[t]**2) /\
    #             (blk.cold_fluid.properties_in[t].dens_mass * blk.channel_diameter) +\
    #             0.7 * (blk.pass_num * Gpc[t]**2 /
    #                    blk.cold_fluid.properties_in[t].dens_mass) + \
    #             blk.cold_fluid.properties_in[t].dens_mass * \
    #             CONST.acceleration_gravity * \
    #             (blk.plate_length + blk.port_diameter)

    #     self.deltaP_c = Expression(self.flowsheet().time,
    #                                rule=rule_coldside_dP,
    #                                doc='Coldside pressure drop  [Pa]')

    #     def rule_eq_deltaP_hot(blk, t):
    #         return blk.hot_fluid.deltaP[t] == -blk.deltaP_h[t]
    #     self.eq_deltaP_hot = Constraint(self.flowsheet().time,
    #                                     rule=rule_eq_deltaP_hot)

    #     def rule_eq_deltaP_cold(blk, t):
    #         return blk.cold_fluid.deltaP[t] == -blk.deltaP_c[t]
    #     self.eq_deltaP_cold = Constraint(self.flowsheet().time,
    #                                      rule=rule_eq_deltaP_cold)

    #     # ---------------------------------------------------------------------
    #     # Overall heat transfer coefficients
    #     def rule_U(blk, t, p):
    #         return 1.0 /\
    #             (1.0 / h_hot[t, p] +
    #              blk.plate_gap / blk.plate_thermal_cond +
    #              1.0 / h_cold[t, p])
    #     self.U = Expression(self.flowsheet().time, self.PH,
    #                         rule=rule_U,
    #                         doc='Overall heat transfer coefficient')
    #     # ---------------------------------------------------------------------
    #     # capacitance of hot and cold fluid

    #     def rule_Caph(blk, t, p):
    #         return blk.hot_fluid.properties_in[t].flow_mass * blk.cp_hot[t] / blk.Np[p]
    #     self.Caph = Expression(self.flowsheet().time, self.PH,
    #                            rule=rule_Caph,
    #                            doc='Hotfluid capacitance rate')

    #     def rule_Capc(blk, t, p):
    #         return blk.cold_fluid.properties_in[t].flow_mass * blk.cp_cold[t] / blk.Np[p]
    #     self.Capc = Expression(self.flowsheet().time, self.PH,
    #                            rule=rule_Capc,
    #                            doc='Coldfluid capacitance rate')

    #     # ---------------------------------------------------------------------
    #     # min n max capacitance and capacitance ratio
    #     def rule_Cmin(blk, t, p):
    #         return smooth_min(blk.Caph[t, p], blk.Capc[t, p], eps=blk.epsilon_param)

    #     self.Cmin = Expression(self.flowsheet().time, self.PH,
    #                            rule=rule_Cmin,
    #                            doc='Minimum capacitance rate')

    #     # smooth core/
    #     def rule_Cmax(blk, t, p):
    #         return smooth_max(blk.Caph[t, p], blk.Capc[t, p], eps=blk.epsilon_param)

    #     self.Cmax = Expression(self.flowsheet().time, self.PH,
    #                            rule=rule_Cmax,
    #                            doc='Maximum capacitance rate')

    #     def rule_CR(blk, t, p):
    #         return blk.Cmin[t, p] / blk.Cmax[t, p]
    #     self.CR = Expression(self.flowsheet().time, self.PH,
    #                          rule=rule_CR,
    #                          doc='Capacitance ratio')

    #     # ---------------------------------------------------------------------
    #     # Number of Transfer units for sub heat exchanger
    #     def rule_NTU(blk, t, p):
    #         return blk.U[t, p] * blk.plate_area / blk.Cmin[t, p]
    #     self.NTU = Expression(self.flowsheet().time, self.PH,
    #                           rule=rule_NTU,
    #                           doc='Number of Transfer Units')

    #     # ---------------------------------------------------------------------
    #     # effectiveness of sub-heat exchangers
    #     def rule_Ecf(blk, t, p):
    #         if blk.pass_num.value % 2 == 0:
    #             return (1 - exp(-blk.NTU[t, p] * (1 - blk.CR[t, p]))) / \
    #                 (1 - blk.CR[t, p] *
    #                  exp(-blk.NTU[t, p] * (1 - blk.CR[t, p])))
    #         elif blk.pass_num.value % 2 == 1:
    #             return ((1 - exp(-blk.NTU[t, p] * (1 + blk.CR[t, p]))) /
    #                     (1 + blk.CR[t, p]))

    #     self.Ecf = Expression(self.flowsheet().time, self.PH,
    #                           rule=rule_Ecf,
    #                           doc='Effectiveness for sub-HX')

    #     # ---------------------------------------------------------------------
    #     # Energy balance equations for hot fluid in sub-heat exhanger
    #     def rule_hotside_energy_balance_per_pass(blk, t, p):
    #         return blk.Th_out[t, p] == blk.Th_in[t, p] -\
    #             blk.Ecf[t, p] * blk.Cmin[t, p] / blk.Caph[t, p] * \
    #             (blk.Th_in[t, p] - blk.Tc_in[t, p])

    #     self.Ebh_eq = Constraint(
    #         self.flowsheet().time,
    #         self.PH,
    #         rule=rule_hotside_energy_balance_per_pass,
    #         doc='Hot fluid sub-heat exchanger energy balance')

    #     # Hot fluid exit temperature
    #     def rule_Tout_hot(blk, t):
    #         return blk.Th_out[t, blk.pass_num.value] ==\
    #             blk.hot_fluid.properties_out[t].temperature

    #     self.Tout_hot_eq = Constraint(self.flowsheet().time,
    #                                   rule=rule_Tout_hot,
    #                                   doc='Hot fluid exit temperature')

    #     # Energy balance equations for cold fluid in sub-heat exhanger
    #     def rule_coldside_energy_balance_per_pass(blk, t, p):
    #         return blk.Tc_out[t, p] == blk.Tc_in[t, p] + \
    #             blk.Ecf[t, p] * blk.Cmin[t, p] / blk.Capc[t, p] * \
    #             (blk.Th_in[t, p] - blk.Tc_in[t, p])

    #     self.Ebc_eq = Constraint(
    #         self.flowsheet().time,
    #         self.PH,
    #         rule=rule_coldside_energy_balance_per_pass,
    #         doc='Cold fluid sub-heat exchanger energy balance')

    #     # Cold fluid exit temperature
    #     def rule_Tout_cold(blk, t):
    #         return blk.Tc_out[t, 1] ==\
    #             blk.cold_fluid.properties_out[t].temperature

    #     self.Tout_cold_eq = Constraint(self.flowsheet().time,
    #                                    rule=rule_Tout_cold,
    #                                    doc='Cold fluid exit temperature')

    #     # ---------------------------------------------------------------------
    #     # Energy balance boundary conditions
    #     def rule_hot_BCIN(blk, t):
    #         return blk.Th_in[t, 1] == \
    #             blk.hot_fluid.properties_in[t].temperature

    #     self.hot_BCIN = Constraint(self.flowsheet().time,
    #                                rule=rule_hot_BCIN,
    #                                doc='Hot fluid inlet boundary conditions')

    #     def rule_cold_BCIN(blk, t):
    #         return blk.Tc_in[t, blk.pass_num.value] ==\
    #             blk.cold_fluid.properties_in[t].temperature
    #     self.cold_BCIN = Constraint(self.flowsheet().time,
    #                                 rule=rule_cold_BCIN,
    #                                 doc='Cold fluid inlet boundary conditions')

    #     # list of  passes : 1,...., p-1
    #     Pset = [i for i in range(1, self.pass_num.value)]

    #     def rule_hot_BC(blk, t, p):
    #         return blk.Th_out[t, p] == blk.Th_in[t, p + 1]
    #     self.hot_BC = Constraint(
    #         self.flowsheet().time,
    #         Pset,
    #         rule=rule_hot_BC,
    #         doc='Hot fluid boundary conditions: change of pass')

    #     def rule_cold_BC(blk, t, p):
    #         return blk.Tc_out[t, p + 1] == blk.Tc_in[t, p]
    #     self.cold_BC = Constraint(
    #         self.flowsheet().time, Pset,
    #         rule=rule_cold_BC,
    #         doc='Cold fluid boundary conditions: change of pass')

    #     # ---------------------------------------------------------------------
    #     # Energy transferred from Hotside to Coldside
    #     def rule_heat_lost(blk, t):
    #         return blk.hot_fluid.properties_in[t].flow_mass * blk.cp_hot[t] *\
    #             (blk.hot_fluid.properties_in[t].temperature -
    #              blk.hot_fluid.properties_out[t].temperature)
    #     self.heat_lost = Expression(self.flowsheet().time, rule=rule_heat_lost,
    #                                 doc='Heat lost by hot fluid')

    #     def rule_heat_gain(blk, t):
    #         return blk.cold_fluid.properties_in[t].flow_mass * blk.cp_cold[t] *\
    #             (blk.cold_fluid.properties_out[t].temperature -
    #              blk.cold_fluid.properties_in[t].temperature)
    #     self.heat_gain = Expression(self.flowsheet().time, rule=rule_heat_gain,
    #                                 doc='Heat gain by cold fluid')

    # def initialize(blk, hot_side_state_args=None, cold_side_state_args=None,
    #                outlvl=idaeslog.NOTSET, solver=None, optarg=None):
    #     '''
    #     Initialisation routine for PHE unit (default solver ipopt)

    #     Keyword Arguments:
    #         state_args : a dict of arguments to be passed to the property
    #                        package(s) to provide an initial state for
    #                        initialization (see documentation of the specific
    #                        property package) (default = {}).
    #         outlvl : sets output level of initialization routine
    #         optarg : solver options dictionary object (default=None, use
    #                  default solver options)
    #         solver : str indicating which solver to use during
    #                  initialization (default = None)
    #     Returns:
    #         None
    #     '''
    #     # Set solver options
    #     init_log = idaeslog.getInitLogger(blk.name, outlvl, tag='unit')
    #     solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

    #     # Create solver
    #     opt = get_solver(solver, optarg)

    #     hot_fluid_mole_frac_comp_dict = dict()
    #     cold_fluid_mole_frac_comp_dict = dict()

    #     for i in blk.config.hot_side.property_package.apparent_species_set:
    #         hot_fluid_mole_frac_comp_dict[i] =\
    #             value(blk.hot_inlet.mole_frac_comp[0, i])

    #     for i in blk.config.cold_side.property_package.apparent_species_set:
    #         cold_fluid_mole_frac_comp_dict[i] =\
    #             value(blk.cold_inlet.mole_frac_comp[0, i])

    #     hotside_state_args = {
    #         'flow_mol': value(blk.hot_inlet.flow_mol[0]),
    #         'temperature': value(blk.hot_inlet.temperature[0]),
    #         'pressure': value(blk.hot_inlet.pressure[0]),
    #         'mole_frac_comp': hot_fluid_mole_frac_comp_dict}

    #     coldside_state_args = {
    #         'flow_mol': value(blk.cold_inlet.flow_mol[0]),
    #         'temperature': value(blk.cold_inlet.temperature[0]),
    #         'pressure': value(blk.cold_inlet.pressure[0]),
    #         'mole_frac_comp': cold_fluid_mole_frac_comp_dict}

    #     # ---------------------------------------------------------------------
    #     # Initialize the Inlet properties
    #     init_log.info('Step 1: Property Initialization')
    #     init_log.info_high("Inlet Properties initialization")
    #     blk.hot_fluid.properties_in.initialize(
    #         state_args=hotside_state_args,
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         hold_state=True)
    #     blk.cold_fluid.properties_in.initialize(
    #         state_args=coldside_state_args,
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         hold_state=True)

    #     # Initialize the Outlet properties
    #     init_log.info_high("Outlet Properties initialization")
    #     blk.hot_fluid.properties_out.initialize(
    #         state_args=hotside_state_args,
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         hold_state=False)

    #     blk.cold_fluid.properties_out.initialize(
    #         state_args=coldside_state_args,
    #         outlvl=outlvl,
    #         optarg=optarg,
    #         solver=solver,
    #         hold_state=False)
    #     # ----------------------------------------------------------------------
    #     init_log.info('Step 2: PHE Initialization')
    #     with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
    #         res = opt.solve(blk, tee=slc.tee)
    #     init_log.info_high(
    #         "Step 2 Complete: {}.".format(idaeslog.condition(res)))
    #     init_log.info('Initialization Completed')

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
