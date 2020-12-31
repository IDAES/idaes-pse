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
IDAES Continous Differential Contactor (CDC) Model, MEA GEN1 PackedColumn

Schematic Diagram:

        Clean Gas  <----+     +<---- Lean Solvent
                        |     |
                       +----  -+ --- L,  z=1.0
                       |      /|
                       |     / |
                       |    /  |
                       |  ABS  |
                       |  /    |
                       | /     |
                       |/      |
                       +-------+ --- 0, z=0.0
                        |     |
        Flue gas ------>+     +---->Rich Solvent

    z: dimensionless distance in axial direction
    L: column length
    L*z: distance in axial direction
    The column height can be changed without changing the discretiation
    by changing L.


Detailed model equations can be found in the supplimentary information
of the paper :
Akula, Paul; Eslick, John; Bhattacharyya, Debangsu; Miller, David
"Model Development, Validation, and Part-Load Optimization of a
MEA-Based Post-Combustion CO2 Capture Process
Under Part-Load and Variable Capture Operation,
Industrial & Engineering Chemistry Research,2020. (submitted)
"""
# Import Python libraries and third-party
import os
import numpy as np
import matplotlib.pyplot as plt

# Import Pyomo libraries
from pyomo.environ import (Constraint, Expression, Param, Reals, NonNegativeReals,
                           Set, value, Var, exp, SolverFactory, SolverStatus,
                           units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES Libraries
from idaes.core.util.constants import Constants as CONST
from idaes.core import (ControlVolume1DBlock, UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import add_object_reference
from idaes.core.control_volume1d import DistributedVars
import idaes.logger as idaeslog
from idaes.core.util import to_json, from_json, StoreSpec

__author__ = "Paul Akula, John Eslick"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PackedColumn")
class PackedColumnData(UnitModelBlockData):
    """Standard Continous Differential Contactor (CDC)  Model Class."""

    # Configuration template for unit level arguments applicable to both phases
    CONFIG = ConfigBlock()

    # Configuration template for phase specific  arguments
    _PhaseCONFIG = ConfigBlock()

    # Unit level config arguments applicable to all phases
    CONFIG.declare("dynamic", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Dynamic model flag",
        doc="""Indicates whether this model will be dynamic or not,
            **default** = useDefault.
            **Valid values:** {
            **useDefault** - get flag from parent (default = False),
            **True** - set as a dynamic model,
            **False** - set as a steady-state model.}"""))
    CONFIG.declare("has_holdup", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
            Must be True if dynamic = True,
            **default** - False.
            **Valid values:** {
            **True** - construct holdup terms,
            **False** - do not construct holdup terms}"""))
    CONFIG.declare("area_definition", ConfigValue(
        default=DistributedVars.variant,
        domain=In(DistributedVars),
        description="Argument for defining form of area variable",
        doc="""Argument defining whether area variable should be spatially
        variant or not. **default** - DistributedVars.uniform.
        **Valid values:** {
        DistributedVars.uniform - area does not vary across spatial domian,
        DistributedVars.variant - area can vary over the domain and is indexed
        by time and space.}"""))

    CONFIG.declare("finite_elements", ConfigValue(
        default=10,
        domain=int,
        description="Number of finite elements length domain",
        doc="""Number of finite elements to use when discretizing length
            domain (default=10)"""))
    CONFIG.declare("length_domain_set", ConfigValue(
        default=[0.0, 1.0],
        domain=list,
        description="Number of finite elements length domain",
        doc="""length_domain_set - (optional) list of point to use to
            initialize a new ContinuousSet if length_domain is not
            provided (default = [0.0, 1.0]).
            domain (default = [0.0, 1.0])"""))
    CONFIG.declare("transformation_method", ConfigValue(
        default="dae.finite_difference",
        description="Method to use for DAE transformation",
        doc="""Method to use to transform domain. Must be a method recognised
            by the Pyomo TransformationFactory,
            **default** - "dae.finite_difference".
            **Valid values:** {
            **"dae.finite_difference"**
                   - Use a finite difference transformation method,
            **"dae.collocation"**
                   - use a collocation transformation method}"""))
    CONFIG.declare("collocation_points", ConfigValue(
        default=3,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
            discretizing length domain (default=3)"""))
    CONFIG.declare("flow_type", ConfigValue(
        default="counter_current",
        domain=In(['counter_current']),
        description="Flow configuration of PackedColumn",
        doc="""Flow configuration of PackedColumn
                - counter_current: gas side flows from 0 to 1
                                   liquid side flows from 1 to 0"""))
    _PhaseCONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.componentTotal,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
        **default** - MaterialBalanceType.componentTotal.
        **Valid values:** {
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
        **default** - EnergyBalanceType.enthalpyTotal.
        **Valid values:** {
        **EnergyBalanceType.none** - exclude energy balances,
        **EnergyBalanceType.enthalpyTotal**
                 - single enthalpy balance for material,
        **EnergyBalanceType.enthalpyPhase**
                 - enthalpy balances for each phase,
        **EnergyBalanceType.energyTotal**
                 - single energy balance for material,
        **EnergyBalanceType.energyPhase**
                 - energy balances for each phase.}"""))
    CONFIG.declare("process_type", ConfigValue(
        default='Absorber',
        domain=In(['Absorber', 'Stripper', 'Regenerator']),
        description="Flag indicating the type of  process",
        doc="""Flag indicating either absoprtion or stripping process.
             Absorber process has O2 and N2 in the vapor phase
                """))

    # Populate the phase side template to default values
    _PhaseCONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.none,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
            **default** - MomentumBalanceType.pressureTotal.
            **Valid values:** {
            **MomentumBalanceType.none** - exclude momentum balances,
            **MomentumBalanceType.pressureTotal**
                    - single pressure balance for material,
            **MomentumBalanceType.pressurePhase**
                    - pressure balances for each phase,
            **MomentumBalanceType.momentumTotal**
                    - single momentum balance for material,
            **MomentumBalanceType.momentumPhase**
                    - momentum balances for each phase.}"""))
    _PhaseCONFIG.declare("has_mass_transfer", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Mass transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
            **default** - False.
            **Valid values:** {
            **True** - include mass transfer terms}"""))
    _PhaseCONFIG.declare("has_heat_transfer", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
            **default** - False.
            **Valid values:** {
            **True** - include heat transfer terms}"""))
    _PhaseCONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
            constructed,
            **default** - False.
            **Valid values:** {
            **True** - include pressure change terms,
            **False** - exclude pressure change terms.}"""))
    _PhaseCONFIG.declare("pressure_drop_type", ConfigValue(
        default=None,
        domain=In([None, "Billet_Schultes_correlation",
                   "Stichlmair_Fair_Bravo_correlation",
                   "GPDC-Kister"]),
        description="Construction flag for type of pressure drop",
        doc="""Indicates what type of pressure drop correlation should be used,
            **default** - None
            **Valid values:** {
            **None** - set pressure drop to zero,
            **"Stichlmair_Fair_Bravo_correlation"**
                      - Use the Stichlmair_Fair_Bravo_correlation model
            **"GPDC-Kister"**
                - Use the Generalized Pressure Drop Correlation of Kister 2007
            **"Billet_Schultes_correlation"**
                     - Use the Billet_Schultes_correlation model}"""))
    _PhaseCONFIG.declare("has_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Phase equilibrium term construction flag",
        doc="""Argument to enable phase equilibrium on the gas side.
            - True - include phase equilibrium term
            - False - do not include phase equilirium term"""))

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

    _PhaseCONFIG.declare("transformation_scheme", ConfigValue(
        default="BACKWARD",
        description="Scheme to use for DAE transformation",
        doc="""Scheme to use when transformating domain. See Pyomo
            documentation for supported schemes,
            **default** - "BACKWARD".
            **Valid values:** {
            **"BACKWARD"**
                 - Use a BACKWARD finite difference transformation method,
            **"FORWARD""**
                 - Use a FORWARD finite difference transformation method,
            **"LAGRANGE-RADAU""**
                 - use a collocation transformation method}"""))
    # Create individual config blocks for vapor(gas) and liquid sides
    CONFIG.declare("vapor_side",
                   _PhaseCONFIG(doc="vapor side config arguments"))
    CONFIG.declare("liquid_side",
                   _PhaseCONFIG(doc="liquid side config arguments"))

    # Set vapor side momentum balance to pressureTotal as default
    CONFIG.vapor_side.momentum_balance_type = MomentumBalanceType.pressureTotal

    #==========================================================================

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to build default attributes
        super(PackedColumnData, self).build()

    #==========================================================================
        """ Set argument values for vapor and liquid sides"""
        # Consistency check for process type
        if (self.config.process_type
                not in ['Absorber', 'Stripper', 'Regenerator']):
            raise ConfigurationError("{} got an invalid value for process_type "
                                     " process type must be"
                                     " 'Absorber', 'Stripper', or 'Regenerator'"
                                     .format(self.name))

        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, Liquid flows from 1 to 0
        if self.config.flow_type == "counter_current":
            set_direction_vapor = FlowDirection.forward
            set_direction_liquid = FlowDirection.backward

        # Set  material balance to componentTotal as default
        self.config.liquid_side.material_balance_type = \
            MaterialBalanceType.componentTotal
        self.config.vapor_side.material_balance_type = \
            MaterialBalanceType.componentTotal

        # set as rate-based
        self.config.liquid_side.has_mass_transfer = True
        self.config.vapor_side.has_mass_transfer = True

        # Set  energy_balance_type to enthalpyTotal as default
        self.config.energy_balance_type == EnergyBalanceType.enthalpyTotal

        # Set heat transfer terms for rate-based model
        self.config.vapor_side.has_heat_transfer = True
        self.config.liquid_side.has_heat_transfer = True

    # ==========================================================================
        """ Build Control volume 1D for vapor phase and
            populate vapor control volume"""

        self.vapor_phase = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.vapor_side.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": self.config.area_definition,
            "property_package": self.config.vapor_side.property_package,
            "property_package_args":
                self.config.vapor_side.property_package_args})

        self.vapor_phase.add_geometry(flow_direction=set_direction_vapor,
                                      length_domain_set=self.config.length_domain_set)

        self.vapor_phase.add_state_blocks(
            information_flow=set_direction_vapor,
            has_phase_equilibrium=False)

        self.vapor_phase.add_material_balances(
            balance_type=self.config.vapor_side.material_balance_type,
            has_phase_equilibrium=self.config.vapor_side.has_phase_equilibrium,
            has_mass_transfer=self.config.vapor_side.has_mass_transfer)

        self.vapor_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.vapor_side.has_heat_transfer)

        self.vapor_phase.add_momentum_balances(
            balance_type=self.config.vapor_side.momentum_balance_type,
            has_pressure_change=self.config.vapor_side.has_pressure_change)

        self.vapor_phase.apply_transformation()

    # ==========================================================================
        """ Build Control volume 1D for liquid phase and
            populate liquid control volume"""
        self.liquid_phase = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme": self.config.liquid_side.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": self.config.area_definition,
            "property_package": self.config.liquid_side.property_package,
            "property_package_args":
                self.config.liquid_side.property_package_args})

        self.liquid_phase.add_geometry(flow_direction=set_direction_liquid,
                                       length_domain_set=self.config.
                                       length_domain_set)

        self.liquid_phase.add_state_blocks(
            information_flow=set_direction_liquid,
            has_phase_equilibrium=False)

        self.liquid_phase.add_material_balances(
            balance_type=self.config.liquid_side.material_balance_type,
            has_phase_equilibrium=self.config.liquid_side.has_phase_equilibrium,
            has_mass_transfer=self.config.liquid_side.has_mass_transfer)

        self.liquid_phase.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.liquid_side.has_heat_transfer)

        self.liquid_phase.apply_transformation()

        # Add Ports for vapor side
        self.add_inlet_port(name="vapor_inlet", block=self.vapor_phase)
        self.add_outlet_port(name="vapor_outlet", block=self.vapor_phase)

        # Add Ports for liquid side
        self.add_inlet_port(name="liquid_inlet", block=self.liquid_phase)
        self.add_outlet_port(name="liquid_outlet", block=self.liquid_phase)

    # ==========================================================================
        """ Add performace equation method"""
        self._make_performance()

    def _make_performance(self):
        """
        Constraints for unit model.

        Args:
            None

        Returns:
            None
        """

        # ======================================================================
        # Add object references - Sets
        add_object_reference(self,
                             "vap_comp",
                             self.config.vapor_side.property_package.component_list)
        add_object_reference(self,
                             "liq_comp",
                             self.config.liquid_side.property_package.component_list)
        add_object_reference(self,
                             "vapor_phase_list_ref",
                             self.config.vapor_side.property_package.phase_list)
        add_object_reference(self,
                             "liquid_phase_list_ref",
                             self.config.liquid_side.property_package.phase_list)

        # Add object reference - time
        add_object_reference(self,
                             "t",
                             self.flowsheet().config.time)

        # Add object references - Transport parameters
        add_object_reference(self,
                             "eps_ref",
                             self.config.vapor_side.property_package.eps_p)
        add_object_reference(self,
                             "a_ref",
                             self.config.vapor_side.property_package.a)

        add_object_reference(self,
                             "Cl_ref",
                             self.config.vapor_side.property_package.Cl)
        add_object_reference(self,
                             "Cv_ref",
                             self.config.vapor_side.property_package.Cv)
        add_object_reference(self,
                             "dh_ref",
                             self.config.vapor_side.property_package.dia_hydraulic)

        # Add object references - others
        R_ref = CONST.gas_constant

        # Unit Model Parameters/sets
        self.zi = Param(self.vapor_phase.length_domain, mutable=True,
                        doc='''integer indexing parameter required for transfer
                             across boundaries of a given volume element''')
        # store the integer  indexing
        for i, x in enumerate(self.vapor_phase.length_domain, 1):
            self.zi[x] = i

        self.homotopy_par_m = Param(initialize=0, mutable=True, units=None,
                                    doc='''continuation parameter to turn on mass
                                    transfer terms gradually''')
        self.homotopy_par_h = Param(initialize=0, mutable=True, units=None, doc='''continuation parameter to turn on heat
                                    transfer terms gradually''')

        # diffusing components
        self.dcomp = Set(initialize=['CO2', 'H2O'])

        # Unit Model Variables
        # Geometry
        self.dia_col = Var(domain=Reals, initialize=0.1, units=pyunits.m,
                           doc='Column diameter')
        self.area_col = Var(domain=Reals, initialize=0.5, units=pyunits.m**2,
                            doc='Column cross-sectional area')
        self.length_col = Var(domain=Reals, initialize=4.9, units=pyunits.m,
                              doc='Column length')

        # Vapor Inlet Conditions
        self.vap_in_flow = Var(self.t, domain=Reals, units=pyunits.m / pyunits.s,
                               doc='Inlet vapor molar flow rate')
        self.vap_in_temperature = Var(self.t, domain=Reals, units=pyunits.K,
                                      doc='Temperature  at vapor inlet')
        self.vap_in_mole_frac = Var(self.t, self.vap_comp, domain=Reals, units=None,
                                    doc='Mole fractions of  '
                                    'species at vapor inlet')

        # Liquid Inlet Conditions
        self.liq_in_flow = Var(self.t, domain=Reals, units=pyunits.mol / pyunits.s,
                               doc='Inlet liquid molar flow rate')
        self.liq_in_temperature = Var(self.t, domain=Reals, units=pyunits.K,
                                      doc='Temperature  at liquid inlet')
        self.liq_in_mole_frac = Var(self.t, self.liq_comp, domain=Reals, units=None,
                                    doc='Mole fraction of '
                                    'species at solvent inlet')

        # Hydrodynamics
        self.bot_pressure = Var(self.t, domain=Reals, initialize=101325, units=pyunits.Pa,
                                doc='Column Bottom pressure')
        self.vel_vap = Var(self.t,
                           self.vapor_phase.length_domain,
                           domain=NonNegativeReals, initialize=2, units=pyunits.m / pyunits.s,
                           doc='Vapor superficial velocity')
        self.vel_liq = Var(self.t,
                           self.liquid_phase.length_domain, units=pyunits.m / pyunits.s,
                           domain=NonNegativeReals, initialize=0.01,
                           doc='Liquid superficial velocity')
        # mass and heat transfer terms
        # mass transfer
        self.pressure_equil = Var(self.t,
                                  self.vapor_phase.length_domain,
                                  self.dcomp,
                                  domain=NonNegativeReals, initialize=500, units=pyunits.Pa,
                                  doc='''Equilibruim pressure of diffusing
                                      components at the interface ''')
        self.N_v = Var(self.t,
                       self.liquid_phase.length_domain,
                       self.dcomp,
                       domain=Reals, initialize=0.0, units=pyunits.mol / (pyunits.s * pyunits.m**3),
                       doc='''Moles of diffusing species transfered
                                     into liquid ''')
        self.E = Var(self.t,
                     self.liquid_phase.length_domain, units=None,
                     domain=NonNegativeReals, initialize=160,
                     doc='Enhancement factor')
        self.yi_MEA = Var(self.t,
                          self.liquid_phase.length_domain,
                          domain=NonNegativeReals, initialize=0.5, units=None,
                          doc='''Dimensionless concentration of MEA
                                    at interface ''')
        self.yeq_CO2 = Var(self.t,
                           self.liquid_phase.length_domain,
                           domain=NonNegativeReals, initialize=0.5, units=None,
                           doc='''Dimensionless concentration of CO2
                                      in equilibruim with the bulk''')

        # heat transfer
        self.Q_v = Var(self.t,
                       self.vapor_phase.length_domain,
                       domain=Reals, initialize=0.0, units=pyunits.J / (pyunits.s * pyunits.m**3),
                       doc='Heat transfer rate in vapor phase')
        self.Q_l = Var(self.t,
                       self.vapor_phase.length_domain,
                       domain=Reals, initialize=0.0, units=pyunits.J / (pyunits.s * pyunits.m**3),
                       doc='Heat transfer rate in liquid phase')

        # =====================================================================
        # Add performance equations

        # inter-facial Area model:
        # reference: Tsai correlation,regressed by Chinen et al. 2018
        def rule_interfacial_area(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                return blk.a_ref * 0.6486 * (
                    blk.liquid_phase.properties[t, x].dens_mass /
                    blk.liquid_phase.properties[t, x].surf_tens *
                    (blk.vel_liq[t, x])**(4.0 / 3.0))**0.12

        self.area_interfacial = Expression(self.t,
                                           self.vapor_phase.length_domain,
                                           rule=rule_interfacial_area,
                                           doc='inter-facial area [m2/m3]')

        # liquid holdup model
        # reference: Tsai correlation,regressed by Chinen et al. 2018
        def rule_holdup_liq(blk, t, x):
            return 24.2355 * (blk.vel_liq[t, x] *
                              (blk.liquid_phase.properties[t, x].visc_d /
                               blk.liquid_phase.properties[t, x].dens_mass)**(1.0 / 3.0))**0.6471

        self.holdup_liq = Expression(self.t,
                                     self.liquid_phase.length_domain,
                                     rule=rule_holdup_liq,
                                     doc='volumetric liquid holdup [-]')

        # vapor holdup model
        # reference: Tsai correlation,regressed by Chinen et al. 2018
        def rule_holdup_vap(blk, t, x):
            return blk.eps_ref - blk.holdup_liq[t, x]

        self.holdup_vap = Expression(self.t,
                                     self.vapor_phase.length_domain,
                                     rule=rule_holdup_vap,
                                     doc='volumetric vapor holdup [-]')

        # ---------------------------------------------------------------------
        # Geometry contraints

        # Column area
        @self.Constraint(doc="Column cross-sectional area [m2]")
        def column_cross_section_area(blk):
            return blk.area_col == (CONST.pi * 0.25 * (blk.dia_col)**2)

        # Area of control volume : vapor side and liquid side
        control_volume_area_definition = ''' column_area * phase_holdup.
        The void fraction of the vapor phase (volumetric vapor holdup) and that
        of the liquid phase(volumetric liquid holdup) are
        lumped into the definition of the cross-sectional area of the
        vapor-side and liquid-side control volume respectively. Hence, the
        cross-sectional area of the control volume changes with time and space.
        '''

        @self.Constraint(self.t,
                         self.vapor_phase.length_domain,
                         doc=control_volume_area_definition)
        def vapor_side_area(bk, t, x):
            return bk.vapor_phase.area[t, x] == bk.area_col * bk.holdup_vap[t, x]

        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc=control_volume_area_definition)
        def liquid_side_area(bk, t, x):
            return bk.liquid_phase.area[t, x] == bk.area_col * bk.holdup_liq[t, x]

        # Length of control volume : vapor side and liquid side
        @self.Constraint(doc="Vapor side length")
        def vapor_side_length(blk):
            return blk.vapor_phase.length == blk.length_col

        @self.Constraint(doc="Liquid side length")
        def liquid_side_length(blk):
            return blk.liquid_phase.length == blk.length_col

        # ---------------------------------------------------------------------
        # Hydrodynamic contraints
        # Vapor superficial velocity
        @self.Constraint(self.t,
                         self.vapor_phase.length_domain,
                         doc="Vapor superficial velocity")
        def eq_vel_vap(blk, t, x):
            return blk.vel_vap[t, x] * blk.area_col * \
                blk.vapor_phase.properties[t, x].conc_mol == \
                blk.vapor_phase.properties[t, x].flow_mol

        # Liquid superficial velocity
        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc="Liquid superficial velocity")
        def eq_vel_liq(blk, t, x):
            return blk.vel_liq[t, x] * blk.area_col * \
                blk.liquid_phase.properties[t, x].conc_mol == \
                blk.liquid_phase.properties[t, x].flow_mol

        # pressure drop calculation
        if (self.config.vapor_side.has_pressure_change and
            self.config.vapor_side.pressure_drop_type ==
                "Stichlmair_Fair_Bravo_correlation"):
            raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "pressure drop. Please contact the "
                "developer of the property_package you are using."
                .format(self.name))

        if (self.config.vapor_side.has_pressure_change and
                self.config.vapor_side.pressure_drop_type == "GPDC-Kister"):
            raise NotImplementedError(
                "{} control volume class has not implemented a method for "
                "pressure drop. Please contact the "
                "developer of the property_package you are using."
                .format(self.name))

        # ---------------------------------------------------------------------
        # sum  of mole fraction to unity contraints:use component flows
        # Sum of component  flows  in vapor phase
        @self.Constraint(self.t,
                         self.vapor_phase.length_domain,
                         doc="Sum of vapor component flows")
        def eq_sumy(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return blk.vapor_phase.properties[t, x].flow_mol == \
                    sum(blk.vapor_phase.properties[t, x].flow_mol_comp[j]
                        for j in blk.vap_comp)

        # Sum of component  flows  in liquid phase
        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc="Sum of liquid component flows")
        def eq_sumx(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return blk.liquid_phase.properties[t, x].flow_mol == \
                    sum(blk.liquid_phase.properties[t, x].flow_mol_comp[j]
                        for j in blk.liq_comp)

        # ---------------------------------------------------------------------
        # mass and heat transfer coefficients

        # vapor  mass transfer coefficients for diffusing components
        def rule_mass_transfer_coeff_vap(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                return 1 /\
                    (R_ref * blk.vapor_phase.properties[t, x].temperature) *\
                    blk.Cv_ref / (blk.holdup_vap[t, x])**0.5 *\
                    (blk.a_ref / blk.dh_ref)**0.5 *\
                    (blk.vapor_phase.properties[t, x].diffus[j])**(2 / 3) *\
                    (blk.vapor_phase.properties[t, x].visc_d /
                        blk.vapor_phase.properties[t, x].dens_mass)**(1 / 3) *\
                    ((blk.vel_vap[t, x] * blk.vapor_phase.properties[t, x].dens_mass) /
                        (blk.a_ref * blk.vapor_phase.properties[t, x].visc_d))**(3 / 4)

        self.k_v = Expression(self.t,
                              self.vapor_phase.length_domain,
                              self.dcomp,
                              rule=rule_mass_transfer_coeff_vap,
                              doc=' Vapor mass transfer coefficient [mol/m2.s.Pa]')

        # mass transfer coefficients of CO2 in liquid phase
        def rule_mass_transfer_coeff_CO2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return blk.Cl_ref * 12**(1 / 6) * (blk.vel_liq[t, x] *
                                                   blk.liquid_phase.properties[t, x].diffus['CO2'] /
                                                   (blk.dh_ref * blk.holdup_liq[t, x]))**0.5

        self.k_l_CO2 = Expression(self.t,
                                  self.liquid_phase.length_domain,
                                  rule=rule_mass_transfer_coeff_CO2,
                                  doc='''CO2 mass transfer coefficient in solvent
                                     [m/s]''')
        # mass tranfer terms

        def rule_phi(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                zb = self.vapor_phase.length_domain[self.zi[x].value - 1]
                return blk.E[t, zb] * blk.k_l_CO2[t, zb] / blk.k_v[t, x, 'CO2']

        self.phi = Expression(self.t,
                              self.vapor_phase.length_domain,
                              rule=rule_phi,
                              doc='''CO2 Equilibruim partial pressure
                                   intermediate  term''')

        # Equilibruim partial pressure of diffusing components at interface
        @self.Constraint(self.t,
                         self.vapor_phase.length_domain,
                         self.dcomp,
                         doc='''Equilibruim partial pressure of diffusing
                                components at interface''')
        def pressure_at_interface(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.pressure_equil[t, x, j] == 0.0
            else:
                zb = self.vapor_phase.length_domain[self.zi[x].value - 1]
                if j == 'H2O':
                    return blk.pressure_equil[t, x, j] == (
                        blk.liquid_phase.properties[t, zb].vol_mol *
                        blk.liquid_phase.properties[t, zb].conc_mol_comp_true[j] *
                        blk.liquid_phase.properties[t, zb].pressure_sat[j])
                elif j == 'CO2':
                    return blk.pressure_equil[t, x, j] == (
                        (blk.vapor_phase.properties[t, x].mole_frac_comp[j] *
                         blk.vapor_phase.properties[t, x].pressure + blk.phi[t, x] *
                         blk.liquid_phase.properties[t, zb].conc_mol_comp_true[j]) /
                        (1 + blk.phi[t, x] /
                         blk.liquid_phase.properties[t, zb].henry_N2O_analogy))

        # mass transfer of  diffusing components
        def rule_mass_transfer(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.N_v[t, x, j] == 0.0
            else:
                return blk.N_v[t, x, j] == (blk.k_v[t, x, j] *
                                            blk.area_interfacial[t, x] * blk.area_col *
                                            (blk.vapor_phase.properties[t, x].mole_frac_comp[j] *
                                             blk.vapor_phase.properties[t, x].pressure -
                                             blk.pressure_equil[t, x, j])) * blk.homotopy_par_m

        self.mass_transfer = Constraint(self.t,
                                        self.vapor_phase.length_domain,
                                        self.dcomp, rule=rule_mass_transfer,
                                        doc="mass transfer to liquid")

        # mass tranfer term handle
        # liquid side
        if self.config.liquid_side.has_mass_transfer:
            @self.Constraint(self.t,
                             self.liquid_phase.length_domain,
                             self.liquid_phase_list_ref,
                             self.liq_comp,
                             doc="mass transfer to liquid")
            def liquid_phase_mass_transfer_handle(blk, t, x, p, j):
                if x == self.liquid_phase.length_domain.last():
                    return blk.liquid_phase.mass_transfer_term[t, x, p, j] == 0.0
                else:
                    zf = self.vapor_phase.length_domain[self.zi[x].value + 1]
                    if j == 'MEA':
                        return blk.liquid_phase.mass_transfer_term[t, x, p, j] == \
                            0.0
                    else:
                        return blk.liquid_phase.mass_transfer_term[t, x, p, j] == \
                            blk.N_v[t, zf, j]
        # vapor side
        if self.config.vapor_side.has_mass_transfer:
            @self.Constraint(self.t,
                             self.vapor_phase.length_domain,
                             self.vapor_phase_list_ref,
                             self.vap_comp,
                             doc="mass transfer from vapor")
            def vapor_phase_mass_transfer_handle(blk, t, x, p, j):
                if x == self.vapor_phase.length_domain.first():
                    return blk.vapor_phase.mass_transfer_term[t, x, p, j] == 0.0
                else:
                    if j in ['N2', 'O2']:
                        return blk.vapor_phase.mass_transfer_term[t, x, p, j] == \
                            0.0
                    else:
                        return blk.vapor_phase.mass_transfer_term[t, x, p, j] == \
                            -blk.N_v[t, x, j]

        # Vapor-liquid heat transfer coefficient, Chilton Colburn  analogy
        def rule_heat_transfer_coeff(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                return blk.k_v[t, x, 'CO2'] *\
                    blk.vapor_phase.properties[t, x].pressure *\
                    blk.vapor_phase.properties[t, x].cp_mol_mean *\
                    (blk.vapor_phase.properties[t, x].therm_cond /
                     (blk.vapor_phase.properties[t, x].conc_mol *
                      blk.vapor_phase.properties[t, x].cp_mol_mean *
                      blk.vapor_phase.properties[t, x].diffus['CO2']))**(2 / 3)

        self.h_v = Expression(self.t,
                              self.vapor_phase.length_domain,
                              rule=rule_heat_transfer_coeff,
                              doc='''vap-liq heat transfer coefficient
                                     [J/m2.s.K]''')

        # Vapor-liquid heat transfer coefficient modified by Ackmann factor
        def rule_heat_transfer_coeff_Ack(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                Ackmann_factor =\
                    (blk.vapor_phase.properties[t, x].cp_mol_comp_mean['CO2'] *
                     blk.N_v[t, x, 'CO2'] +
                     blk.vapor_phase.properties[t, x].cp_mol_comp_mean['H2O'] *
                     blk.N_v[t, x, 'H2O'])
                return Ackmann_factor /\
                    (1 - exp(-Ackmann_factor /
                             (blk.h_v[t, x] * blk.area_interfacial[t, x] * blk.area_col)))
        self.h_v_Ack = Expression(self.t,
                                  self.vapor_phase.length_domain,
                                  rule=rule_heat_transfer_coeff_Ack,
                                  doc='''vap-liq heat transfer coefficient corrected
                                     by Ackmann factor [J/m3.s.K]''')

        # heat transfer vapor  side
        @self.Constraint(self.t,
                         self.vapor_phase.length_domain,
                         doc="heat transfer - vapor side [J/s.m]")
        def vapor_phase_heat_transfer(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return blk.Q_v[t, x] == 0
            else:
                zb = self.vapor_phase.length_domain[value(self.zi[x]) - 1]
                return blk.Q_v[t, x] == blk.h_v_Ack[t, x] * \
                    (blk.liquid_phase.properties[t, zb].temperature -
                     blk.vapor_phase.properties[t, x].temperature) * \
                    blk.homotopy_par_h

        # heat transfer liquid side
        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc="heat transfer - liquid side [J/s.m]")
        def liquid_phase_heat_transfer(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.Q_l[t, x] == 0
            else:
                zf = self.vapor_phase.length_domain[value(self.zi[x]) + 1]
                return blk.Q_l[t, x] == blk.Q_v[t, zf] + \
                    (blk.liquid_phase.properties[t, x].habs * blk.N_v[t, zf, 'CO2'] -
                     blk.liquid_phase.properties[t, x].hvap * blk.N_v[t, zf, 'H2O']) *\
                    blk.homotopy_par_h

        # heat transfer handle
        if self.config.energy_balance_type != EnergyBalanceType.none:
            # vapor  heat transfer handle
            @self.Constraint(self.t,
                             self.vapor_phase.length_domain,
                             doc="vapor - heat transfer handle")
            def vapor_phase_heat_transfer_handle(blk, t, x):
                return blk.vapor_phase.heat[t, x] == blk.Q_v[t, x]

            # liquid  heat transfer handle
            @self.Constraint(self.t,
                             self.liquid_phase.length_domain,
                             doc="liquid - heat transfer handle")
            def liquid_phase_heat_transfer_handle(blk, t, x):
                return blk.liquid_phase.heat[t, x] == -blk.Q_l[t, x]

        # Enhancement factor model
        # reference: Jozsef Gaspar,Philip Loldrup Fosbol, (2015)
        # self.yi_MEA[z] is equivalent to sqrt(yi_MEA) in the document

        def rule_conc_mol_comp_interface_CO2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                zf = self.liquid_phase.length_domain[self.zi[x].value + 1]
                return blk.pressure_equil[t, zf, 'CO2'] /\
                    blk.liquid_phase.properties[t, x].henry_N2O_analogy

        self.conc_mol_comp_CO2_eq = Expression(self.t,
                                               self.liquid_phase.length_domain,
                                               rule=rule_conc_mol_comp_interface_CO2,
                                               doc='''Concentration of CO2
                                                at the interface [mol/m3]''')

        def rule_Hatta(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return (blk.liquid_phase.properties[t, x].k2_rxn *
                        blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEA'] *
                        blk.liquid_phase.properties[t, x].diffus['CO2'])**0.5 /\
                    blk.k_l_CO2[t, x]

        self.Hatta = Expression(self.t,
                                self.liquid_phase.length_domain,
                                rule=rule_Hatta,
                                doc='Hatta number')

        def rule_yb_CO2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return blk.liquid_phase.properties[t, x].conc_mol_comp_true['CO2'] /\
                    blk.conc_mol_comp_CO2_eq[t, x]

        self.yb_CO2 = Expression(self.t,
                                 self.liquid_phase.length_domain,
                                 rule=rule_yb_CO2,
                                 doc='''Dimensionless concentration of CO2,
                                      Driving force term where
                                      Absortion implies yb_CO2 < 1 and
                                      Desorption impies yb_CO2 > 1 ''')

        def rule_instantaneous_E(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return 1 + (blk.liquid_phase.properties[t, x].diffus['MEA'] *
                            blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEA']) /\
                    (2 * blk.liquid_phase.properties[t, x].diffus['CO2'] *
                        blk.conc_mol_comp_CO2_eq[t, x])

        self.instant_E = Expression(self.t,
                                    self.liquid_phase.length_domain,
                                    rule=rule_instantaneous_E,
                                    doc='Instantaneous Enhancement factor')

        def rule_yi_MEACOO(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return 1 + \
                    (blk.liquid_phase.properties[t, x].diffus['MEA'] *
                     blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEA']) * \
                    (1 - blk.yi_MEA[t, x] * blk.yi_MEA[t, x]) / \
                    (2 * blk.liquid_phase.properties[t, x].diffus['MEACOO-'] *
                        blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEACOO-'])

        self.yi_MEACOO = Expression(self.t,
                                    self.liquid_phase.length_domain,
                                    rule=rule_yi_MEACOO,
                                    doc='Dimensionless concentration of MEACOO-')

        def rule_yi_MEAH(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return 1 + \
                    (blk.liquid_phase.properties[t, x].diffus['MEA'] *
                     blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEA']) * \
                    (1 - blk.yi_MEA[t, x] * blk.yi_MEA[t, x]) / \
                    (2 * blk.liquid_phase.properties[t, x].diffus['MEA+'] *
                        blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEA+'])

        self.yi_MEAH = Expression(self.t,
                                  self.liquid_phase.length_domain,
                                  rule=rule_yi_MEAH,
                                  doc='Dimensionless concentration of MEA+')

        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc='''dimensionless concentration of CO2
                                at equilibruim with the bulk ''')
        def yeq_CO2_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.yeq_CO2[t, x] == 0.0
            else:
                return blk.yeq_CO2[t, x] * blk.yi_MEA[t, x]**4 == \
                    blk.yb_CO2[t, x] * blk.yi_MEAH[t, x] * blk.yi_MEACOO[t, x]

        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc='''Enhancement factor model Eqn 1 ''')
        def E1_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.E[t, x] == 1
            else:
                return (blk.E[t, x] - 1) * (1 - blk.yb_CO2[t, x]) == \
                    (blk.instant_E[t, x] - 1) * (1 - blk.yi_MEA[t, x]**2)

        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc='''Enhancement factor model Eqn 2 ''')
        def E2_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.yi_MEA[t, x] == 0
            else:
                return blk.E[t, x] * (1 - blk.yb_CO2[t, x]) == \
                    blk.Hatta[t, x] * blk.yi_MEA[t, x] * \
                    (1 - blk.yeq_CO2[t, x])

        @self.Constraint(self.t,
                         self.liquid_phase.length_domain,
                         doc='Enhancement factor lower bound ')
        def E3_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return 1 - blk.E[t, x] <= 0.0

        # ---------------------------------------------------------------------
        # Boundary Conditions
        # Enthalpy balance

        @self.Constraint(self.t,
                         doc='Boundary condition for vapor temperature @ inlet ')
        def vap_BC_temperature(blk, t):
            return blk.vap_in_temperature[t] ==\
                blk.vapor_phase.properties[t, 0].temperature

        @self.Constraint(self.t,
                         doc='Boundary condition for liquid temperature @ inlet ')
        def liq_BC_temperature(blk, t):
            return blk.liq_in_temperature[t] ==\
                blk.liquid_phase.properties[t, 1].temperature

        # Mass balance-mole_frac
        @self.Constraint(self.t, self.liq_comp,
                         doc='Boundary condition for liquid mole fraction @ inlet ')
        def liq_BC_mole_frac(blk, t, j):
            return blk.liq_in_mole_frac[t, j] ==\
                blk.liquid_phase.properties[t, 1].mole_frac_comp[j]

        @self.Constraint(self.t, self.vap_comp,
                         doc='Boundary condition for vapor mole fraction @ inlet ')
        def vap_BC_mole_frac(blk, t, j):
            return blk.vap_in_mole_frac[t, j] ==\
                blk.vapor_phase.properties[t, 0].mole_frac_comp[j]

        # Mass balance-flow-mol
        @self.Constraint(self.t,
                         doc='Boundary condition for vapor flow @ inlet ')
        def vap_BC_flow_mol(blk, t):
            return blk.vap_in_flow[t] ==\
                blk.vapor_phase.properties[t, 0].flow_mol

        @self.Constraint(self.t,
                         doc='Boundary condition for liquid flow @ inlet ')
        def liq_BC_flow_mol(blk, t):
            return blk.liq_in_flow[t] ==\
                blk.liquid_phase.properties[t, 1].flow_mol

        # Pressure balance
        @self.Constraint(self.t,
                         doc='Boundary condition for  pressure @ inlet ')
        def vap_BC_pressure(blk, t):
            return blk.bot_pressure[t] ==\
                blk.vapor_phase.properties[t, 0].pressure

    #==========================================================================
    # Model initialization routine
    def initialize(blk, method='default', *args, **kwargs):
        """
        Initialize the mea column.

        Args:
            method: Method of initialization to use in {**'default'** -
                bootstrap initilization, **json** - reload a saved model state
                from a json file
            *args: Positional args fed to the selected initialization function,
            **kwargs: keyword arguments fed to selected initialization function}
        """
        if method == 'default':
            cwd = os.getcwd()
            initialization_folder = os.path.join(cwd, r"initialized_files")
            if not os.path.exists(initialization_folder):
                os.makedirs(initialization_folder)

            # get names for json files
            if blk.config.dynamic:
                s1 = 'dyn'
            else:
                s1 = 'ss'
            if blk.config.process_type == 'Absorber':
                s2 = 'abs'
            else:
                s2 = 'reg'

            fname_init =\
                os.path.join(initialization_folder,
                             "{}_{}.json".format(s1, s2))

            if os.path.exists(fname_init):
                from_json(blk, fname=fname_init,
                          wts=StoreSpec(ignore_missing=True))
                if not blk.config.dynamic:
                    blk.make_steady_state_column_profile()
                if blk.config.dynamic:
                    blk.make_dynamic_column_profile()
            else:
                blk._default_init(*args, **kwargs)
                # save result
                fname =\
                    os.path.join(initialization_folder,
                                 "{}_{}.json".format(s1, s2))
                to_json(blk, fname=fname)
        else:
            raise(Exception("Unknown initialization method {}".format(method)))

    def _default_init(blk, vapor_phase_state_args={}, liquid_phase_state_args={},
                      state_vars_fixed=False,
                      homotopy_steps_m=[0, 0.1, 0.2,
                                        0.3, 0.4, 0.5, 0.6, 0.8, 1],
                      homotopy_steps_h=[0, 0.1, 0.2, 0.3,
                                        0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                      outlvl=idaeslog.NOTSET, solver='ipopt', optarg={'tol': 1e-6}):
        """
        Initialisation routine for MEA PackedColumn unit (default solver ipopt).

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            homotopy_steps_m : List of continuations steps [0 -1]
                               No mass transfer --->mass transfer
            homotopy_steps_h : List of continuations steps [0 -1]
                               No heat transfer ---> heat transfer
            outlvl : sets output level of initialisation routine

                     * 0 = no output (default)
                     * 1 = return solver state for each step in routine
                     * 2 = return solver state for each step in subroutines
                     * 3 = include solver output infomation (tee=True)

            optarg : solver options dictionary object (default={'tol': 1e-6})
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
#
#        Returns:
#            None
#        """
        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        opt = SolverFactory(solver)
        opt.options = optarg

        dynamic_constraints = [
            "vapor_side_area",
            "liquid_side_area",
            "eq_vel_vap",
            "eq_vel_liq",
            "eq_sumy",
            "eq_sumx",
            "pressure_at_interface",
            "mass_transfer",
            "liquid_phase_mass_transfer_handle",
            "vapor_phase_mass_transfer_handle",
            "vapor_phase_heat_transfer",
            "liquid_phase_heat_transfer",
            "vapor_phase_heat_transfer_handle",
            "liquid_phase_heat_transfer_handle",
            "yeq_CO2_eqn",
            "E1_eqn"
            "E2_eqn",
            #"E3_eqn",
            "vap_BC_temperature",
            "liq_BC_temperature",
            "liq_BC_mole_frac",
            "vap_BC_mole_frac",
            "vap_BC_flow_mol",
            "liq_BC_flow_mol",
            "vap_BC_pressure",
            "material_balances",
            "material_flow_linking_constraints",
            "material_holdup_calculation",
            "enthalpy_flow_linking_constraint",
            "energy_holdup_calculation",
            "material_flow_dx_disc_eq",
            "enthalpy_flow_dx_disc_eq",
            "pressure_dx_disc_eq",
            "enthalpy_balances",
            "material_accumulation_disc_eq",
            "energy_accumulation_disc_eq",
            "pressure_balance"]

        steady_state_constraints = [
            "vapor_side_area",
            "liquid_side_area",
            "eq_vel_vap",
            "eq_vel_liq",
            "eq_sumy",
            "eq_sumx",
            "pressure_at_interface",
            "mass_transfer",
            "liquid_phase_mass_transfer_handle",
            "vapor_phase_mass_transfer_handle",
            "vapor_phase_heat_transfer",
            "liquid_phase_heat_transfer",
            "vapor_phase_heat_transfer_handle",
            "liquid_phase_heat_transfer_handle",
            "yeq_CO2_eqn",
            "E1_eqn"
            "E2_eqn",
            #"E3_eqn",
            "vap_BC_temperature",
            "liq_BC_temperature",
            "liq_BC_mole_frac",
            "vap_BC_mole_frac",
            "vap_BC_flow_mol",
            "liq_BC_flow_mol",
            "vap_BC_pressure",
            "material_balances",
            "material_flow_linking_constraints",
            "enthalpy_flow_linking_constraint",
            "material_flow_dx_disc_eq",
            "enthalpy_flow_dx_disc_eq",
            "pressure_dx_disc_eq",
            "enthalpy_balances",
            "pressure_balance"]

        # ---------------------------------------------------------------------
        # Deactivate unit model level constraints (asides geometry constraints)
        for c in blk.component_objects(Constraint, descend_into=True):
            if blk.config.dynamic:
                if c.local_name in dynamic_constraints:
                    c.deactivate()
            elif not blk.config.dynamic:
                if c.local_name in steady_state_constraints:
                    c.deactivate()

        # For some reasons E1_eqn & E2_eqn  are still active
        # So deactivating again.
        blk.E1_eqn.deactivate()
        blk.E2_eqn.deactivate()
        # Fix some variables
        # Hydrodynamics - velocity
        blk.vel_liq[:, :].fix()
        blk.vel_vap[:, :].fix()
        # interface pressure
        blk.pressure_equil[:, :, :].fix()
        # flux
        blk.N_v[:, :, :].fix(0.0)
        blk.vapor_phase.mass_transfer_term.fix(0.0)
        blk.liquid_phase.mass_transfer_term.fix(0.0)
        # Enhancement factor model
        blk.E[:, :].fix()
        blk.yi_MEA[:, :].fix()
        blk.yeq_CO2[:, :].fix()
        # heat transfer
        blk.Q_v[:, :].fix(0.0)
        blk.Q_l[:, :].fix(0.0)
        blk.vapor_phase.heat[:, :].fix(0.0)
        blk.liquid_phase.heat[:, :].fix(0.0)
        # area
        blk.vapor_phase.area.fix(1)
        blk.liquid_phase.area.fix(1)

        # other variables
        # Pressure_dx
        blk.vapor_phase.pressure_dx[:, :].fix(0.0)
        # vapor side flow terms
        blk.vapor_phase._enthalpy_flow[:, :, :].fix(1.0)
        blk.vapor_phase.enthalpy_flow_dx[:, :, :].fix(0.0)
        blk.vapor_phase._flow_terms[:, :, :, :].fix(1.0)
        blk.vapor_phase.material_flow_dx[:, :, :, :].fix(0.0)
        # liquid side flow terms
        blk.liquid_phase._enthalpy_flow[:, :, :].fix(1.0)
        blk.liquid_phase.enthalpy_flow_dx[:, :, :].fix(0.0)
        blk.liquid_phase._flow_terms[:, :, :, :].fix(1.0)
        blk.liquid_phase.material_flow_dx[:, :, :, :].fix(0.0)
        # accumulation terms
        # fix accumulation terms to zero and holdup to 1
        if blk.config.dynamic:
            # liquid
            blk.liquid_phase.energy_holdup[:, :, :].fix(1.0)
            blk.liquid_phase.energy_accumulation[:, :, :].fix(0.0)
            blk.liquid_phase.material_holdup[:, :, :, :].fix(1.0)
            blk.liquid_phase.material_accumulation[:, :, :, :].fix(0.0)
            # vapor
            blk.vapor_phase.energy_holdup[:, :, :].fix(1.0)
            blk.vapor_phase.energy_accumulation[:, :, :].fix(0.0)
            blk.vapor_phase.material_holdup[:, :, :, :].fix(1.0)
            blk.vapor_phase.material_accumulation[:, :, :, :].fix(0.0)

        # ---------------------------------------------------------------------
        # get values for state variables for initialization
        if blk.config.process_type == 'Absorber':
            vapor_phase_state_args = {
                'flow_mol': blk.vap_in_flow[0].value,
                'temperature': blk.vap_in_temperature[0].value,
                'pressure': blk.bot_pressure[0].value,
                'mole_frac_comp':
                {'H2O': blk.vap_in_mole_frac[0, 'H2O'].value,
                 'CO2': blk.vap_in_mole_frac[0, 'CO2'].value,
                 'N2': blk.vap_in_mole_frac[0, 'N2'].value,
                 'O2': blk.vap_in_mole_frac[0, 'O2'].value}}
        elif blk.config.process_type in ['Stripper', 'Regenerator']:
            vapor_phase_state_args = {
                'flow_mol': blk.vap_in_flow[0].value,
                'temperature': blk.vap_in_temperature[0].value,
                'pressure': blk.bot_pressure[0].value,
                'mole_frac_comp':
                {'H2O': blk.vap_in_mole_frac[0, 'H2O'].value,
                 'CO2': blk.vap_in_mole_frac[0, 'CO2'].value}}

        liquid_phase_state_args = {
            'flow_mol': blk.liq_in_flow[0].value,
            'temperature': blk.liq_in_temperature[0].value,
            'pressure': blk.bot_pressure[0].value,
            'mole_frac_comp':
            {'H2O': blk.liq_in_mole_frac[0, 'H2O'].value,
             'CO2': blk.liq_in_mole_frac[0, 'CO2'].value,
             'MEA': blk.liq_in_mole_frac[0, 'MEA'].value}}

        print('==============================================================')
        print("STEP 1: Property Package initialization")
        print('==============================================================')

        # Initialize vapor_phase block
        blk.vapor_phase.properties.initialize(
            state_args=vapor_phase_state_args,
            state_vars_fixed=False,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)
        # Initialize liquid_phase properties block
        blk.liquid_phase.properties.initialize(
            state_args=liquid_phase_state_args,
            state_vars_fixed=False,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)

        print('==============================================================')
        print("STEP 2: Steady-State ISOTHERMAL MASS BALANCE")
        print('        -->no mass transfer ')
        print('        -->no heat transfer')

        # unfix flow variable terms
        # vapor side
        blk.vapor_phase.properties[:, :].flow_mol.unfix()
        blk.vapor_phase.properties[:, :].mole_frac_comp[:].unfix()
        blk.vapor_phase._flow_terms[:, :, :, :].unfix()
        blk.vapor_phase.material_flow_dx[:, :, :, :].unfix()
        # liquid-side
        blk.liquid_phase.properties[:, :].flow_mol.unfix()
        blk.liquid_phase.properties[:, :].mole_frac_comp[:].unfix()
        blk.liquid_phase._flow_terms[:, :, :, :].unfix()
        blk.liquid_phase.material_flow_dx[:, :, :, :].unfix()

        # activate mass balance related equations
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in [
                "eq_sumy",
                "eq_sumx",
                "liq_BC_mole_frac",
                "vap_BC_mole_frac",
                "vap_BC_flow_mol",
                "liq_BC_flow_mol",
                "material_balances",
                "material_flow_linking_constraints",
                    "material_flow_dx_disc_eq"]:
                c.activate()

        # solve for a small lenght if stripper
        if (blk.config.process_type == 'Stripper' or
                blk.config.process_type == 'Regenerator'):
            specified_length = value(blk.length_col)
            blk.length_col.fix(0.6)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info(
            "Step 2: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        print('==============================================================')
        print('STEP 3: Add Mass tranfer terms')
        print('        -->(3a) velocities & Interface pressure')
        print('        -->(3b) Enhancement factor')

        # Initialize :
        #   --> Velocities
        #   --> Interface pressure
        #   --> Enhancement factor

        # velocity
        blk.vel_vap.unfix()
        blk.vel_liq.unfix()
        blk.eq_vel_vap.activate()
        blk.eq_vel_liq.activate()
        for t in blk.t:
            for x in blk.vapor_phase.length_domain:
                blk.vel_vap[t, x].value = value(
                    blk.vapor_phase.properties[t, x].flow_mol /
                    (blk.area_col * blk.vapor_phase.properties[t, x].conc_mol))
            for x in blk.liquid_phase.length_domain:
                blk.vel_liq[t, x].value = value(
                    blk.liquid_phase.properties[t, x].flow_mol /
                    (blk.area_col * blk.liquid_phase.properties[t, x].conc_mol))

        # Interface pressure
        blk.pressure_equil.unfix()
        blk.pressure_at_interface.activate()
        blk.E.fix(1)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info(
            "Step 3a: {}.".format(idaeslog.condition(res)))
        # ----------------------------------------------------------------------
        # Enhancement factor model
        blk.E.unfix()
        for t in blk.t:
            for x in blk.liquid_phase.length_domain:
                blk.E[t, x].value = 100
        blk.yi_MEA.unfix()
        blk.yeq_CO2.unfix()
        blk.E1_eqn.activate()
        blk.E2_eqn.activate()
        blk.yeq_CO2_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info("Step 3 complete: {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------

        print('============================================================')
        print('STEP 4: Isothermal chemical absoption')
        print("Homotopy steps: ")
        print("No mass transfer (0.0) --> (1.0) mass transfer")

        # ISOTHERMAL CHEMICAL ABSORPTION
        blk.N_v.unfix()
        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()
        blk.mass_transfer.activate()
        blk.vapor_phase_mass_transfer_handle.activate()
        blk.liquid_phase_mass_transfer_handle.activate()

        for i in homotopy_steps_m:
            print('homotopy step -->{0:5.2f}'.format(i))
            blk.homotopy_par_m = i
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
                if res.solver.status != SolverStatus.warning:
                    print('')

        init_log.info("Step 4 complete: {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------

        print('==============================================================')
        print('STEP 5: Adiabatic chemical absoption')
        print("Homotopy steps:")
        print("Isothermal (0.0) --> (1.0) Adiabatic ")

        # Unfix temperature
        blk.liquid_phase.properties[:, :].temperature.unfix()
        blk.vapor_phase.properties[:, :].temperature.unfix()
        # unfix heat transfer terms
        blk.Q_v.unfix()
        blk.Q_l.unfix()
        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()
        # unfix enthalpy flow variable terms
        blk.vapor_phase._enthalpy_flow[:, :, :].unfix()
        blk.vapor_phase.enthalpy_flow_dx[:, :, :].unfix()
        blk.liquid_phase._enthalpy_flow[:, :, :].unfix()
        blk.liquid_phase.enthalpy_flow_dx[:, :, :].unfix()

        # activate steady-state energy balance related equations
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in [
                "vapor_phase_heat_transfer",
                "liquid_phase_heat_transfer",
                "vapor_phase_heat_transfer_handle",
                "liquid_phase_heat_transfer_handle",
                "vap_BC_temperature",
                "liq_BC_temperature",
                "enthalpy_flow_linking_constraint",
                "enthalpy_flow_dx_disc_eq",
                    "enthalpy_balances"]:
                c.activate()

        for i in homotopy_steps_h:
            print('homotopy step -->{0:5.2f}'.format(i))
            blk.homotopy_par_h = i
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)

        init_log.info("Step 5 complete: {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------

        # scale up stripper length
        if (blk.config.process_type == 'Stripper' or
                blk.config.process_type == 'Regenerator'):
            packing_height = np.linspace(0.6, specified_length, num=10)
            print('SCALEUP Stripper height')
            for L in packing_height:
                blk.length_col.fix(L)
                print('Packing height = {:6.2f}'.format(L))
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(blk, tee=slc.tee)
                init_log.info("Scaleup: {}.".format(idaeslog.condition(res)))

        if not blk.config.dynamic:
            blk.make_steady_state_column_profile()
            print('=============STEADY-STATE INITIALIZATION COMPLETE=========')

        if blk.config.dynamic:
            print('==========================================================')
            print('STEP 6: unfix Accumulation and Holdup terms')
            print("        --->6a Holdup calculations")
            print("        --->6b Include Accumulation terms")
            # activate holdup constraints
            for c in blk.component_objects(Constraint, descend_into=True):
                if c.local_name in [
                    "vapor_side_area",
                    "liquid_side_area",
                    "material_holdup_calculation",
                        "energy_holdup_calculation"]:
                    c.activate()

            # unfix holdup terms
            blk.vapor_phase.energy_holdup[:, :, :].unfix()
            blk.vapor_phase.material_holdup[:, :, :, :].unfix()
            blk.liquid_phase.energy_holdup[:, :, :].unfix()
            blk.liquid_phase.material_holdup[:, :, :, :].unfix()

            # unfix  CV1D area lumped with phase volumetric holdup fraction
            blk.vapor_phase.area.unfix()
            blk.liquid_phase.area.unfix()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
            init_log.info("Step 6a complete: {}.".format(
                idaeslog.condition(res)))

            # Step 6b:
            # unfix accumulation terms(derivative variables)
            blk.vapor_phase.energy_accumulation[:, :, :].unfix()
            blk.vapor_phase.material_accumulation[:, :, :, :].unfix()
            blk.liquid_phase.energy_accumulation[:, :, :].unfix()
            blk.liquid_phase.material_accumulation[:, :, :, :].unfix()

            # activate constraints for accumulation terms
            for c in blk.component_objects(Constraint, descend_into=True):
                if c.local_name in [
                    "material_accumulation_disc_eq",
                        "energy_accumulation_disc_eq"]:
                    c.activate()

            blk.fix_initial_condition()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
            init_log.info("Step 6 complete: {}.".format(
                idaeslog.condition(res)))

            print('===================INITIALIZATION COMPLETE=================')

            if not blk.config.dynamic:
                blk.make_steady_state_column_profile()
            if blk.config.dynamic:
                blk.make_dynamic_column_profile()

    def fix_initial_condition(blk):
        '''
        Mass balance : Initial condition  is determined by
        fixing n-1 mole fraction and the total molar flowrate
        Energy balance :Initial condition  is determined by
        fixing  the temperature.
        '''
        for x in blk.vapor_phase.length_domain:
            if x != 0:
                blk.vapor_phase.properties[0, x].temperature.fix()
                blk.vapor_phase.properties[0, x].flow_mol.fix()
            for j in blk.vap_comp:
                if (x != 0 and j != 'CO2'):
                    blk.vapor_phase.properties[0, x].mole_frac_comp[j].fix()
        for x in blk.liquid_phase.length_domain:
            if x != 1:
                blk.liquid_phase.properties[0, x].temperature.fix()
                blk.liquid_phase.properties[0, x].flow_mol.fix()
            for j in blk.liq_comp:
                if (x != 1 and j != 'CO2'):
                    blk.liquid_phase.properties[0, x].mole_frac_comp[j].fix()

    def make_steady_state_column_profile(blk):
        normalised_column_height = [x for x in blk.vapor_phase.length_domain]
        simulation_time = [t for t in blk.t]
        # final time
        tf = simulation_time[-1]
        CO2_profile = []
        liquid_temperature_profile = []
        # APPEND RESULTS
        for x in blk.vapor_phase.length_domain:
            CO2_profile.append(
                value(1e-3 * blk.vapor_phase.properties[tf, x].pressure *
                      blk.vapor_phase.properties[tf, x].mole_frac_comp['CO2']))
            liquid_temperature_profile.append(
                value(blk.liquid_phase.properties[tf, x].temperature))

        # plot properties
        fontsize = 18
        labelsize = 18
        fig = plt.figure(figsize=(9, 7))
        ax1 = fig.add_subplot(111)
        ax1.set_title('Steady-state column profile',
                      fontsize=16, fontweight='bold')

        # plot primary axis
        lab1 = ax1.plot(normalised_column_height, CO2_profile,
                        linestyle='--', mec="b", mfc="None",
                        color='b', label='CO$_{2}$ partial pressure [kPa]',
                        marker='o')

        ax1.tick_params(axis='y', labelcolor='b',
                        direction='in', labelsize=labelsize)
        ax1.tick_params(axis='x', direction='in', labelsize=labelsize)

        ax1.set_xlabel('Normalise column  height from bottom',
                       fontsize=fontsize)
        ax1.set_ylabel('P$_{CO_{2}}$  [ kPa]', color='b', fontweight='bold',
                       fontsize=fontsize)
        # plot secondary axis
        ax2 = ax1.twinx()
        lab2 = ax2.plot(normalised_column_height,
                        liquid_temperature_profile,
                        color='g',
                        linestyle='-',
                        label='Liquid temperature profile',
                        marker='s')
        ax2.set_ylabel('T$_{liq}$ [ K ] ', color='g', fontweight='bold',
                       fontsize=fontsize)
        ax2.tick_params(axis='y', labelcolor='g',
                        direction='in', labelsize=labelsize)
        # get the labels
        lab_1 = lab1 + lab2
        labels_1 = [l.get_label() for l in lab_1]
        ax1.legend(lab_1, labels_1, loc='lower center', fontsize=fontsize)
        fig.tight_layout()
        # show graph
        plt.show()

    def make_dynamic_column_profile(blk):
        normalised_column_height = [x for x in blk.vapor_phase.length_domain]
        simulation_time = [t for t in blk.t]
        fluegas_flow = [value(blk.vap_in_flow[t]) for t in blk.t]
        # final time
        tf = simulation_time[-1]
        nf = len(simulation_time)
        # mid-time
        if nf % 2 == 0:
            tm = int(nf / 2)
        else:
            tm = int(nf / 2 + 1)

        CO2_profile_mid = []
        CO2_profile_fin = []
        liquid_temperature_profile_mid = []
        liquid_temperature_profile_fin = []

        # APPEND RESULTS
        for x in blk.vapor_phase.length_domain:
            CO2_profile_mid.append(
                value(1e-3 * blk.vapor_phase.properties[tm, x].pressure *
                      blk.vapor_phase.properties[tm, x].mole_frac_comp['CO2']))
            CO2_profile_fin.append(
                value(1e-3 * blk.vapor_phase.properties[tf, x].pressure *
                      blk.vapor_phase.properties[tf, x].mole_frac_comp['CO2']))

            liquid_temperature_profile_mid.append(
                value(blk.liquid_phase.properties[tm, x].temperature))
            liquid_temperature_profile_fin.append(
                value(blk.liquid_phase.properties[tf, x].temperature))

        # plot properties
        fontsize = 18
        labelsize = 18
        fig = plt.figure(figsize=(12, 7))
        ax1 = fig.add_subplot(211)
        ax1.set_title('Column profile @ {0:6.2f} & {1:6.2f} sec'.format(tm, tf),
                      fontsize=16, fontweight='bold')

        # plot primary axis
        lab1 = ax1.plot(normalised_column_height, CO2_profile_mid,
                        linestyle='--', color='b',
                        label='CO$_{2}$ partial pressure [kPa] @ %d' % tm)
        lab2 = ax1.plot(normalised_column_height, CO2_profile_fin,
                        linestyle='-', color='b',
                        label='CO$_{2}$ partial pressure [kPa] @ %d' % tf)

        ax1.tick_params(axis='y', labelcolor='b',
                        direction='in', labelsize=labelsize)
        ax1.tick_params(axis='x', direction='in', labelsize=labelsize)

        ax1.set_xlabel('Normalise column  height from bottom',
                       fontsize=fontsize)
        ax1.set_ylabel('P$_{CO_{2}}$  [ kPa]', color='b', fontweight='bold',
                       fontsize=fontsize)

        # plot secondary axis
        ax2 = ax1.twinx()
        lab3 = ax2.plot(normalised_column_height,
                        liquid_temperature_profile_mid,
                        color='g', linestyle='--',
                        label='Liquid temperature profile @ {0:6.1f}'.format(tm))
        lab4 = ax2.plot(normalised_column_height,
                        liquid_temperature_profile_fin,
                        color='g', linestyle='-',
                        label='Liquid temperature profile @ {0:6.1f}'.format(tf))
        ax2.set_ylabel('T$_{liq}$ [ K ] ', color='g', fontweight='bold',
                       fontsize=fontsize)
        ax2.tick_params(axis='y', labelcolor='g',
                        direction='in', labelsize=labelsize)
        # get the labels
        lab_1 = lab1 + lab2 + lab3 + lab4
        labels_1 = [l.get_label() for l in lab_1]
        ax1.legend(lab_1, labels_1, fontsize=fontsize)

        # plot flowgas flow
        ax3 = fig.add_subplot(212)
        ax3.plot(simulation_time, fluegas_flow,
                 linestyle='--', mec="g", mfc="None",
                 color='g', label='Fluegas flow [mol/s]',
                 marker='o')
        ax3.tick_params(labelsize=labelsize)
        ax3.set_xlabel('Simulation time (sec)', fontsize=fontsize)
        ax3.set_ylabel(' Fv  [ mol/s]', color='b', fontweight='bold',
                       fontsize=fontsize)
        ax3.legend(['Fluegas flow [mol/s]'], fontsize=fontsize)
        fig.tight_layout()
        plt.show()
