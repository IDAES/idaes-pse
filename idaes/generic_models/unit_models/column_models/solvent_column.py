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
General purpose gas-solvent contactor column model.

"""

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           Expression,
                           NonNegativeReals,
                           Param,
                           Reals,
                           units as pyunits,
                           Var)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES Libraries
from idaes.core import (ControlVolume1DBlock,
                        UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection)
from idaes.core.util.constants import Constants
from idaes.core.util import get_solver, scaling as iscale
from idaes.core.util.config import is_physical_parameter_block
from idaes.generic_models.unit_models.heat_exchanger_1D import \
    HeatExchangerFlowPattern as FlowPattern
from idaes.core.control_volume1d import DistributedVars
import idaes.logger as idaeslog


__author__ = "Andrew Lee, Paul Akula, John Eslick, Anuja Deshpande"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PackedColumn")
class PackedColumnData(UnitModelBlockData):
    """
    Standard Continous Differential Contactor (CDC) Model Class.

    """

    # Configuration template for unit level arguments applicable to both phases
    CONFIG = UnitModelBlockData.CONFIG()

    # Configuration template for phase specific  arguments
    _PhaseCONFIG = ConfigBlock()

    CONFIG.declare("finite_elements", ConfigValue(
        default=10,
        domain=int,
        description="Number of finite elements length domain",
        doc="""Number of finite elements to use when discretizing length
domain (default=20)"""))
    CONFIG.declare("length_domain_set", ConfigValue(
        default=[0.0, 1.0],
        domain=list,
        description="Number of finite elements length domain",
        doc="""length_domain_set - (optional) list of point to use to
initialize a new ContinuousSet if length_domain is not
provided (default = [0.0, 1.0])"""))
    CONFIG.declare("flow_type", ConfigValue(
        default=FlowPattern.countercurrent,
        domain=In(FlowPattern),
        description="Flow configuration of PackedColumn",
        doc="""PackedColumn flow pattern,
**default** - FlowPattern.countercurrent.
**Valid values:** {
**FlowPattern.countercurrent** - countercurrent flow,
**FlowPattern.cocurrent** - cocurrent flow}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed, **default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))

    # Populate the phase side template to default values
    _PhaseCONFIG.declare("property_package", ConfigValue(
        default=None,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object"""))
    _PhaseCONFIG.declare("property_package_args", ConfigValue(
        default={},
        description="Arguments for constructing vapor property package",
        doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)"""))

    # Create individual config blocks for vapor(gas) and liquid sides
    CONFIG.declare("vapor_side",
                   _PhaseCONFIG(doc="vapor side config arguments"))

    CONFIG.declare("liquid_side",
                   _PhaseCONFIG(doc="liquid side config arguments"))

    # -------------------------------------------------------------------------

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to build default attributes
        super().build()

        # TODO : Property package validation
        # Same reference state

        # ---------------------------------------------------------------------
        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, Liquid flows from 1 to 0
        if self.config.flow_type == FlowPattern.countercurrent:
            set_direction_vapor = FlowDirection.forward
            set_direction_liquid = FlowDirection.backward
        else:
            raise NotImplementedError(
                "{} Packed Column class has implemented only counter-current "
                "flow pattern. Please contact the "
                "developer of the unit model you are using.".format(self.name))

        # ---------------------------------------------------------------------
        # TODO : Get units from properties
        # Geometric variables and constraints
        self.length_column = Var(domain=Reals,
                                 initialize=4.9,
                                 units=pyunits.m,
                                 doc='Column length')

        self.diameter_column = Var(domain=Reals,
                                   initialize=0.1,
                                   units=pyunits.m,
                                   doc='Column diameter')

        self.area_column = Var(domain=Reals,
                               initialize=0.5,
                               units=pyunits.m**2,
                               doc='Column cross-sectional area')

        @self.Constraint(doc="Column cross-sectional area")
        def column_cross_section_area(blk):
            return blk.area_column == (
                Constants.pi*0.25*blk.diameter_column**2)

        # Packing  parameters
        self.voidage_packing = Param(initialize=0.97,
                                     units=pyunits.dimensionless,
                                     mutable=True,
                                     doc="Packing void space m3/m3")

        self.packing_specific_surface_area = Param(
            initialize=250,
            units=pyunits.m**2/pyunits.m**3,
            mutable=True,
            doc="Packing specific surface area")

        self.diameter_hydraulic = Expression(
            expr=4*self.voidage_packing/self.packing_specific_surface_area,
            doc="Hydraulic diameter")

        # ---------------------------------------------------------------------
        # Build control volumnes
        # Vapor phase
        # TODO: Add enthalpy transfer terms to energy balances
        self.vapor_phase = ControlVolume1DBlock(default={
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "BACKWARD",
            "finite_elements": self.config.finite_elements,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
            "property_package": self.config.vapor_side.property_package,
            "property_package_args":
                self.config.vapor_side.property_package_args})

        self.vapor_phase.add_geometry(
            flow_direction=set_direction_vapor,
            length_domain_set=self.config.length_domain_set,
            length_var=self.length_column)

        self.vapor_phase.add_state_blocks(
            information_flow=set_direction_vapor,
            has_phase_equilibrium=False)

        # TODO: Overrides for default balance types
        self.vapor_phase.add_material_balances(
            balance_type=MaterialBalanceType.useDefault,
            has_phase_equilibrium=False,
            has_mass_transfer=True)

        self.vapor_phase.add_energy_balances(
            balance_type=EnergyBalanceType.useDefault,
            has_heat_transfer=True,
            has_enthalpy_transfer=True)

        self.vapor_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=self.config.has_pressure_change)

        self.vapor_phase.apply_transformation()

        # Liquid phase
        self.liquid_phase = ControlVolume1DBlock(default={
            "transformation_method": "dae.finite_difference",
            "transformation_scheme": "FORWARD",
            "finite_elements": self.config.finite_elements,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
            "property_package": self.config.liquid_side.property_package,
            "property_package_args":
                self.config.liquid_side.property_package_args})

        self.liquid_phase.add_geometry(
            flow_direction=set_direction_liquid,
            length_domain_set=self.config.length_domain_set,
            length_var=self.length_column)

        self.liquid_phase.add_state_blocks(
            information_flow=set_direction_liquid,
            has_phase_equilibrium=False)

        self.liquid_phase.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_phase_equilibrium=False,
            has_mass_transfer=True)

        self.liquid_phase.add_energy_balances(
            balance_type=EnergyBalanceType.enthalpyTotal,
            has_heat_transfer=True,
            has_enthalpy_transfer=True)

        # No liquid phase momentum balance due to assumption of mechanical
        # equilibrium

        self.liquid_phase.apply_transformation()

        # ---------------------------------------------------------------------
        # Add Ports for vapor side
        self.add_inlet_port(name="vapor_inlet", block=self.vapor_phase)
        self.add_outlet_port(name="vapor_outlet", block=self.vapor_phase)

        # Add Ports for liquid side
        self.add_inlet_port(name="liquid_inlet", block=self.liquid_phase)
        self.add_outlet_port(name="liquid_outlet", block=self.liquid_phase)

        # ---------------------------------------------------------------------
        # Vapor and liquid holdup
        # TODO: Link to holdup model
        self.liquid_holdup_fraction = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            initialize=0.03,
            domain=NonNegativeReals,
            bounds=(0, 1),
            units=pyunits.dimensionless,
            doc="Liquid phase volumetric holdup fraction")

        def rule_vapor_holdup(blk, t, x):
            return blk.voidage_packing - blk.liquid_holdup_fraction[t, x]
        self.vapor_holdup_fraction = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_vapor_holdup,
            doc="Vapor phase volumetric holdup fraction")

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="Cross-sectional area occupied by vapor phase")
        def vapor_side_area(blk, t, x):
            return blk.vapor_phase.area[t, x] == (
                blk.area_column * blk.vapor_holdup_fraction[t, x])

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="Cross-sectional area occupied by liquid phase")
        def liquid_side_area(blk, t, x):
            return blk.liquid_phase.area[t, x] == (
                blk.area_column*blk.liquid_holdup_fraction[t, x])

        self.area_interfacial = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=200,
            domain=NonNegativeReals,
            units=pyunits.m**2/pyunits.m**3,
            doc='Specific interfacial area between vapor and liquid phases')

        # ---------------------------------------------------------------------
        # Mechanical equilibrium
        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="Mechanical equilibrium between phases")
        def mechanincal_equilibrium(blk, t, x):
            if x == blk.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return (blk.liquid_phase.properties[t, x].pressure ==
                        blk.vapor_phase.properties[t, x].pressure)

        # ---------------------------------------------------------------------
        # Define component lists for all species and those in equilibrium
        vap_comp = self.config.vapor_side.property_package.component_list
        liq_comp = self.config.liquid_side.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp

        # ---------------------------------------------------------------------
        # Enthalpy transfer due to mass transfer
        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Enthalpy transfer from vapor phase due to mass transfer")
        def vapor_enthalpy_transfer_equation(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            return blk.vapor_phase.enthalpy_transfer[t, x] == sum(
                blk.vapor_phase.mass_transfer_term[t, x, "Vap", j] *
                blk.vapor_phase.properties[t, x].enth_mol_phase_comp["Vap", j]
                for j in equilibrium_comp)

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Equating enthalpy transfer between phases")
        def enthalpy_transfer_balance(blk, t, x):
            if x == blk.vapor_phase.length_domain.first():
                return Constraint.Skip
            xl = blk.liquid_phase.length_domain.prev(x)
            return blk.vapor_phase.enthalpy_transfer[t, x] == \
                -blk.liquid_phase.enthalpy_transfer[t, xl]

        # ---------------------------------------------------------------------
        # Mass transfer
        self.pressure_equil_comp = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            domain=NonNegativeReals,
            initialize=500,
            units=pyunits.Pa,
            doc='Equilibruim pressure of components at the interface')

        self.mass_transfer_coeff_vap_comp = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            domain=NonNegativeReals,
            initialize=10000,
            units=pyunits.mol/pyunits.s/pyunits.m**2/pyunits.Pa,
            doc='Vapor phase mass transfer coefficient')

        def rule_mass_transfer_liq(blk, t, x, j):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            elif j not in equilibrium_comp:
                return blk.liquid_phase.mass_transfer_term[t, x, "Liq", j] == \
                    0.0
            else:
                xv = self.liquid_phase.length_domain.next(x)
                return blk.liquid_phase.mass_transfer_term[t, x, "Liq", j] == (
                    blk.mass_transfer_coeff_vap_comp[t, xv, j] *
                    blk.area_interfacial[t, x] * blk.area_column *
                    (blk.vapor_phase.properties[t, x].fug_phase_comp[
                        "Vap", j] -
                     blk.pressure_equil_comp[t, x, j]))
        self.mass_transfer_liq = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            self.liquid_phase.properties.component_list,
            rule=rule_mass_transfer_liq,
            doc="Mass transfer to liquid phase calculation")

        def rule_mass_transfer_vap(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            elif j not in equilibrium_comp:
                return blk.vapor_phase.mass_transfer_term[t, x, "Vap", j] == \
                    0.0
            else:
                xl = blk.liquid_phase.length_domain.prev(x)
                return blk.vapor_phase.mass_transfer_term[t, x, "Vap", j] == (
                    -blk.liquid_phase.mass_transfer_term[t, xl, "Liq", j])
        self.mass_transfer_vap = Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            self.vapor_phase.properties.component_list,
            rule=rule_mass_transfer_vap,
            doc="Mass transfer to vapor phase calculation")

        # Add equilibrium pressure constraint
        def equil_pressure(blk, t, x, j):
            return blk.pressure_equil_comp[t, x, j] == (
                blk.liquid_phase.properties[t, x].fug_phase_comp["Liq", j])
        self.liquid_fugacity_constraint = Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            equilibrium_comp,
            rule=equil_pressure)

        # ---------------------------------------------------------------------
        # Heat transfer
        self.heat_transfer_coeff = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            domain=NonNegativeReals,
            initialize=5000,
            units=pyunits.J/pyunits.K/pyunits.s/pyunits.m,
            doc='Overall heat transfer coefficient')

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="Heat transfer calculation")
        def vapor_phase_heat_transfer(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                xl = self.liquid_phase.length_domain.prev(x)
                return blk.vapor_phase.heat[t, x] == (
                    blk.heat_transfer_coeff[t, x] *
                    (blk.liquid_phase.properties[t, xl].temperature -
                     blk.vapor_phase.properties[t, x].temperature))

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="Overall energy balance")
        def column_energy_balance(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                xl = self.liquid_phase.length_domain.prev(x)
                return blk.vapor_phase.heat[t, x] == (
                    -blk.liquid_phase.heat[t, xl])

    # -------------------------------------------------------------------------
    # Model scaling
    def calculate_scaling_factors(self):
        # TODO : Improve this
        super().calculate_scaling_factors()

    # -------------------------------------------------------------------------
    # Model initialization routine
    def initialize(blk,
                   vapor_state_args=None,
                   liquid_state_args=None,
                   state_vars_fixed=False,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None):
        """
        Column initialization.

        Arguments:
            vapor_state_args : a dict of arguments to be passed to the vapor
                property package to provide an initial state for initialization
                (see documentation of the specific property package)
                (default = None).
            liquid_state_args : a dict of arguments to be passed to the liquid
                property package to provide an initial state for initialization
                (see documentation of the specific property package)
                (default = None).
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during initialization
                    (default = None, use IDAES default solver)

        """

        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        opt = get_solver(solver, optarg)

        # ---------------------------------------------------------------------
        # Initialize vapor phase control volume block
        vflags = blk.vapor_phase.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=vapor_state_args,
        )

        init_log.info_high('Initialization Step 1 Complete.')

        # ---------------------------------------------------------------------
        # Initialize liquid phase control volume block
        lflags = blk.liquid_phase.initialize(
                outlvl=outlvl,
                optarg=optarg,
                solver=solver,
                state_args=liquid_state_args,
        )

        init_log.info_high('Initialization Step 2 Complete.')

        # ---------------------------------------------------------------------
        # Deactivate heat & mass transfer and solve for intermediate variables
        blk.mass_transfer_liq.deactivate()
        blk.mass_transfer_vap.deactivate()
        blk.vapor_phase_heat_transfer.deactivate()
        blk.column_energy_balance.deactivate()
        blk.vapor_enthalpy_transfer_equation.deactivate()
        blk.enthalpy_transfer_balance.deactivate()

        blk.vapor_phase.mass_transfer_term.fix(0)
        blk.liquid_phase.mass_transfer_term.fix(0)
        blk.vapor_phase.heat.fix(0)
        blk.liquid_phase.heat.fix(0)
        blk.vapor_phase.enthalpy_transfer.fix(0)
        blk.liquid_phase.enthalpy_transfer.fix(0)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 3 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Reactivate mass transfer and solve
        blk.mass_transfer_liq.activate()
        blk.mass_transfer_vap.activate()

        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 4 {}.".format(idaeslog.condition(results))
            )

        # ---------------------------------------------------------------------
        # Reactivate enthalpy transfer and solve
        blk.vapor_enthalpy_transfer_equation.activate()
        blk.enthalpy_transfer_balance.activate()

        blk.vapor_phase.enthalpy_transfer.unfix()
        blk.liquid_phase.enthalpy_transfer.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 5 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Reactivate heat transfer and solve
        blk.vapor_phase_heat_transfer.activate()
        blk.column_energy_balance.activate()

        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            results = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Initialization Step 6 {}.".format(idaeslog.condition(results))
        )

        # ---------------------------------------------------------------------
        # Release states
        blk.vapor_phase.release_state(vflags, outlvl)
        blk.liquid_phase.release_state(lflags, outlvl)

        init_log.info('Initialization Complete: {}'
                      .format(idaeslog.condition(results)))
