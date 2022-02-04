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
IDAES First Generation (GEN1) General Purpose Solvent Column Model.
"""

# Import Python and third-party libraries
import matplotlib.pyplot as plt
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    check_optimal_termination, Constraint, Expression, Param, Reals,
    NonNegativeReals, value, Var, exp, units as pyunits, SolverStatus)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES Libraries
from idaes.core.util.constants import Constants as CONST
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core import (ControlVolume1DBlock, UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection)
from idaes.core.util import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.control_volume1d import DistributedVars
import idaes.logger as idaeslog


__author__ = "Paul Akula, John Eslick, Anuja Deshpande, Andrew Lee"


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
        default=20,
        domain=int,
        description="Number of finite elements length domain",
        doc="""Number of finite elements to use when discretizing length
domain (default=20)"""))

    CONFIG.declare("length_domain_set", ConfigValue(
        default=[0.0, 1.0],
        domain=list,
        description="List of points in length domain",
        doc="""length_domain_set - (optional) list of point to use to
initialize a new ContinuousSet if length_domain is not
provided (default = [0.0, 1.0])"""))

    CONFIG.declare("transformation_method", ConfigValue(
        default="dae.finite_difference",
        description="Method to use for DAE transformation",
        doc="""Method to use to transform domain. Must be a method recognised
by the Pyomo TransformationFactory,
**default** - "dae.finite_difference".
**Valid values:** {
**"dae.finite_difference"** - Use a finite difference transformation method,
**"dae.collocation"** - use a collocation transformation method}"""))

    CONFIG.declare("collocation_points", ConfigValue(
        default=3,
        domain=int,
        description="Number of collocation points per finite element",
        doc="""Number of collocation points to use per finite element when
discretizing length domain (default=3)"""))

    CONFIG.declare("column_pressure_drop", ConfigValue(
        default=0,
        description="Column pressure drop per unit length in Pa/m",
        doc="Column pressure drop per unit length in Pa/m provided as a value or expression"))

    # Populate the phase side template to default values
    _PhaseCONFIG.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=Bool,
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed, **default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))

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
- a dict (see property package for documentation)
            """))
            
    _PhaseCONFIG.declare("transformation_scheme", ConfigValue(
        default="BACKWARD",
        description="Scheme to use for DAE transformation",
        doc="""Scheme to use when transformating domain. See Pyomo
documentation for supported schemes,
**default** - "BACKWARD".
**Valid values:** {
**"BACKWARD"** - Use a BACKWARD finite difference transformation method,
**"FORWARD""** - Use a FORWARD finite difference transformation method,
**"LAGRANGE-RADAU""** - use a collocation transformation method}"""))

    # Create individual config blocks for vapor(gas) and liquid sides
    CONFIG.declare("vapor_side",
                   _PhaseCONFIG(doc="vapor side config arguments"))

    CONFIG.declare("liquid_side",
                   _PhaseCONFIG(doc="liquid side config arguments"))

    # =========================================================================

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

    # =========================================================================
        """ Set argument values for vapor and liquid sides"""

        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, Liquid flows from 1 to 0
        
        # TODO: Only handling countercurrent flow for now.
        set_direction_vapor = FlowDirection.forward
        set_direction_liquid = FlowDirection.backward

    # =========================================================================
        """ Build Control volume 1D for vapor phase and
            populate vapor control volume"""

        self.vapor_phase = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme":
                self.config.vapor_side.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
            "property_package": self.config.vapor_side.property_package,
            "property_package_args":
                self.config.vapor_side.property_package_args})

        self.vapor_phase.add_geometry(
            flow_direction=set_direction_vapor,
            length_domain_set=self.config.length_domain_set)

        self.vapor_phase.add_state_blocks(
            information_flow=set_direction_vapor,
            has_phase_equilibrium=False)

        self.vapor_phase.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_phase_equilibrium=False,
            has_mass_transfer=True)

        self.vapor_phase.add_energy_balances(
            balance_type=EnergyBalanceType.enthalpyTotal,
            has_heat_transfer=True)

        self.vapor_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=self.config.vapor_side.has_pressure_change)

        self.vapor_phase.apply_transformation()

    # ==========================================================================
        """ Build Control volume 1D for liquid phase and
            populate liquid control volume

        """
        self.liquid_phase = ControlVolume1DBlock(default={
            "transformation_method": self.config.transformation_method,
            "transformation_scheme":
                self.config.liquid_side.transformation_scheme,
            "finite_elements": self.config.finite_elements,
            "collocation_points": self.config.collocation_points,
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "area_definition": DistributedVars.variant,
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
            balance_type=MaterialBalanceType.componentTotal,
            has_phase_equilibrium=False,
            has_mass_transfer=True)

        self.liquid_phase.add_energy_balances(
            balance_type=EnergyBalanceType.enthalpyTotal,
            has_heat_transfer=True)

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

        Args: None

        Returns: None

        """

        # ======================================================================
        # Custom Sets
        vap_comp = self.config.vapor_side.property_package.component_list
        liq_comp = self.config.liquid_side.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp
        solvent_comp_list = \
            self.config.liquid_side.property_package.solvent_set
        solute_comp_list = self.config.liquid_side.property_package.solute_set
        vapor_phase_list_ref = \
            self.config.vapor_side.property_package.phase_list
        liquid_phase_list_ref = \
            self.config.liquid_side.property_package.phase_list

        # Packing  parameters
        self.eps_ref = Param(initialize=0.97,units=None,
                             mutable=True,
                             doc="Packing void space m3/m3")

        self.packing_specific_area = Param(initialize=250,units=pyunits.m**2 / pyunits.m**3,
                           mutable=True,
                           doc="Packing specific surface area (m2/m3)")
        
        self.packing_channel_size = Param(initialize=0.1,units=pyunits.m,
                           mutable=True,
                           doc="Packing channel size (m)")
        
        self.hydraulic_diameter = Expression(expr=4 * self.eps_ref / self.packing_specific_area,
                                 doc="Hydraulic diameter (m)")

        # Add the integer indices along vapor phase length domain
        self.zi = Param(self.vapor_phase.length_domain, mutable=True,
                        doc='''Integer indexing parameter required for transfer
                             across boundaries of a given volume element''')
                             
        # Set the integer indices along vapor phase length domain
        for i, x in enumerate(self.vapor_phase.length_domain, 1):
            self.zi[x] = i

        # Unit Model Design Variables
        # Geometry
        self.diameter_column = Var(domain=Reals,
                                   initialize=0.1,
                                   units=pyunits.m,
                                   doc='Column diameter')
        
        self.area_column = Var(domain=Reals,
                               initialize=0.5,
                               units=pyunits.m**2,
                               doc='Column cross-sectional area')
        
        self.length_column = Var(domain=Reals,
                                 initialize=4.9,
                                 units=pyunits.m,
                                 doc='Column length')

        # Hydrodynamics
        self.velocity_vap = Var(self.flowsheet().time,
                                self.vapor_phase.length_domain,
                                domain=NonNegativeReals,
                                initialize=2,
                                units=pyunits.m / pyunits.s,
                                doc='Vapor superficial velocity')
        
        self.velocity_liq = Var(self.flowsheet().time,
                                self.liquid_phase.length_domain,
                                domain=NonNegativeReals,
                                initialize=0.01,
                                units=pyunits.m / pyunits.s,
                                doc='Liquid superficial velocity')
        
        self.holdup_liq = Var(self.flowsheet().time,
                              self.liquid_phase.length_domain,
                              initialize=0.001,
                              doc='Volumetric liquid holdup [-]')

        def rule_holdup_vap(blk, t, x):
            return blk.eps_ref - blk.holdup_liq[t, x]

        self.holdup_vap = Expression(self.flowsheet().time,
                                     self.vapor_phase.length_domain,
                                     rule=rule_holdup_vap,
                                     doc='Volumetric vapor holdup [-]')
        
        # Define gas velocity at flooding point (m/s)                
        self.gas_velocity_flood = Var(self.flowsheet().time,
                                self.vapor_phase.length_domain,
                                initialize=1,
                                doc='Gas velocity at flooding point')
            
        # Flooding fraction 
        def rule_flood_fraction(blk, t, x):
            return blk.velocity_vap[t, x]/blk.gas_velocity_flood[t, x]
            
        self.flood_fraction = Expression(self.flowsheet().time,
                                    self.vapor_phase.length_domain,
                                    rule=rule_flood_fraction,
                                    doc='Flooding fraction (expected to be below 0.8)')
        
        # Mass and heat transfer terms
        
        # Mass transfer terms
        self.pressure_equil = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            domain=NonNegativeReals,
            initialize=500,
            units=pyunits.Pa,
            doc='Equilibruim pressure of diffusing components at interface')
        
        self.interphase_mass_transfer = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            equilibrium_comp,
            domain=Reals,
            initialize=0.1,
            units=pyunits.mol / (pyunits.s * pyunits.m),
            doc='Rate at which moles of diffusing species transfered into liquid')
        
        self.enhancement_factor = Var(self.flowsheet().time,
                                      self.liquid_phase.length_domain,
                                      units=None,
                                      initialize=160,
                                      doc='Enhancement factor')

        # Heat transfer terms
        self.heat_flux_vap = Var(self.flowsheet().time,
                            self.vapor_phase.length_domain,
                            domain=Reals,
                            initialize=0.0,
                            units=pyunits.J / (pyunits.s * (pyunits.m**3)),
                            doc='Volumetric heat flux in vapor phase')

        # =====================================================================
        # Add performance equations

        # Inter-facial Area model ([m2/m3]):

        self.area_interfacial = Var(self.flowsheet().time,
                                    self.vapor_phase.length_domain,
                                    initialize=0.9,
                                    doc='Specific inter-facial area')

        # ---------------------------------------------------------------------
        # Geometry constraints

        # Column area [m2]
        @self.Constraint(doc="Column cross-sectional area")
        def column_cross_section_area(blk):
            return blk.area_column == (
                CONST.pi * 0.25 * (blk.diameter_column)**2)

        # Area of control volume : vapor side and liquid side
        control_volume_area_definition = ''' column_area * phase_holdup.
        The void fraction of the vapor phase (volumetric vapor holdup) and that
        of the liquid phase(volumetric liquid holdup) are
        lumped into the definition of the cross-sectional area of the
        vapor-side and liquid-side control volume respectively. Hence, the
        cross-sectional area of the control volume changes with time and space.
        '''

        if self.config.dynamic:
            @self.Constraint(self.flowsheet().time,
                             self.vapor_phase.length_domain,
                             doc=control_volume_area_definition)
            def vapor_side_area(bk, t, x):
                return bk.vapor_phase.area[t, x] == (
                    bk.area_column * bk.holdup_vap[t, x])

            @self.Constraint(self.flowsheet().time,
                             self.liquid_phase.length_domain,
                             doc=control_volume_area_definition)
            def liquid_side_area(bk, t, x):
                return bk.liquid_phase.area[t, x] == (
                    bk.area_column * bk.holdup_liq[t, x])
        else:
            self.vapor_phase.area.fix(value(self.area_column))
            self.liquid_phase.area.fix(value(self.area_column))

        # Pressure consistency in phases
        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc='''Mechanical equilibruim: vapor-side pressure
                                    equal liquid -side pressure''')
        def mechanical_equil(bk, t, x):
            return bk.liquid_phase.properties[t, x].pressure == \
                    bk.vapor_phase.properties[t, x].pressure

        # Length of control volume : vapor side and liquid side
        @self.Constraint(doc="Vapor side length")
        def vapor_side_length(blk):
            return blk.vapor_phase.length == blk.length_column

        @self.Constraint(doc="Liquid side length")
        def liquid_side_length(blk):
            return blk.liquid_phase.length == blk.length_column

        # ---------------------------------------------------------------------
        # Hydrodynamic constraints
        # Vapor superficial velocity

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="Vapor superficial velocity")
        def eq_velocity_vap(blk, t, x):
            return blk.velocity_vap[t, x] * blk.area_column * \
                blk.vapor_phase.properties[t, x].dens_mol == \
                blk.vapor_phase.properties[t, x].flow_mol

        # Liquid superficial velocity
        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="Liquid superficial velocity")
        def eq_velocity_liq(blk, t, x):
            return blk.velocity_liq[t, x] * blk.area_column * \
                blk.liquid_phase.properties[t, x].dens_mol == \
                blk.liquid_phase.properties[t, x].flow_mol

        # ---------------------------------------------------------------------
        # Mass transfer coefficients
        # Mass transfer coefficients of diffusing components in vapor phase [mol/m2.s.Pa]
        self.k_v = Var(self.flowsheet().time,
                       self.vapor_phase.length_domain,
                       equilibrium_comp,
                       doc=' Vapor phase mass transfer coefficient')

        # Mass transfer coefficients of diffusing components in liquid phase  [m/s]
        self.k_l = Var(self.flowsheet().time,
                       self.liquid_phase.length_domain,
                       equilibrium_comp,
                       doc='Liquid phase mass transfer coefficient')

        # Intermediate term
        def rule_phi(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                zb = self.vapor_phase.length_domain.at(self.zi[x].value - 1)
                return (blk.enhancement_factor[t, zb] *
                        blk.k_l[t, zb, j] /
                        blk.k_v[t, x, j])

        self.phi = Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            solute_comp_list,
            rule=rule_phi,
            doc='Equilibruim partial pressure intermediate term for solute')

        # Equilibruim partial pressure of diffusing components at interface
        @self.Constraint(self.flowsheet().time,
                          self.vapor_phase.length_domain,
                          equilibrium_comp,
                          doc='''Equilibruim partial pressure of diffusing
                                components at interface''')
        def pressure_at_interface(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.pressure_equil[t, x, j] == 0.0
            else:
                zb = self.vapor_phase.length_domain.at(self.zi[x].value - 1)
                lprops = blk.liquid_phase.properties[t, zb]
                henrycomp = lprops.params.get_component(j).config.henry_component
                if henrycomp is not None and "Liq" in henrycomp:
                    return blk.pressure_equil[t, x, j] == (
                        (blk.vapor_phase.properties[t, x].mole_frac_comp[j] *
                          blk.vapor_phase.properties[
                              t, x].pressure + blk.phi[t, x, j] *
                          lprops.conc_mol_phase_comp_true['Liq',j]) /
                        (1 + blk.phi[t, x, j] /
                          blk.liquid_phase.properties[t, zb].henry['Liq',j]))
                else:
                    return blk.pressure_equil[t, x, j] == (
                        lprops.vol_mol_phase['Liq'] *
                        lprops.conc_mol_phase_comp_true['Liq',j] *
                        lprops.pressure_sat_comp[j])

        # Mass transfer of  diffusing components in vapor phase
        def rule_mass_transfer(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.interphase_mass_transfer[t, x, j] == 0.0
            else:
                return blk.interphase_mass_transfer[t, x, j] == (
                blk.k_v[t, x, j] *
                blk.area_interfacial[t, x] * blk.area_column *
                (blk.vapor_phase.properties[t, x].mole_frac_comp[j] *
                 blk.vapor_phase.properties[t, x].pressure -
                 blk.pressure_equil[t, x, j]))

        self.mass_transfer_vapor = Constraint(self.flowsheet().time,
                                        self.vapor_phase.length_domain,
                                        equilibrium_comp,
                                        rule=rule_mass_transfer,
                                        doc="mass transfer in vapor phase")

        # Liquid phase mass transfer handle
        @self.Constraint(self.flowsheet().time,
                          self.liquid_phase.length_domain,
                          self.liquid_phase.properties.phase_component_set,
                          doc="mass transfer to liquid")
        def liquid_phase_mass_transfer_handle(blk, t, x, p, j):
            if x == self.liquid_phase.length_domain.last():
                return blk.liquid_phase.mass_transfer_term[t, x, p, j] == 0.0
            else:
                zf = self.liquid_phase.length_domain.at(self.zi[x].value + 1)
                if j in equilibrium_comp:
                    return blk.liquid_phase.mass_transfer_term[t, x, p, j] == \
                        blk.interphase_mass_transfer[t, zf, j]
                else:
                    return blk.liquid_phase.mass_transfer_term[t, x, p, j] == \
                        0.0

        # Vapor phase mass transfer handle
        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         self.vapor_phase.properties.phase_component_set,
                         doc="mass transfer from vapor")
        def vapor_phase_mass_transfer_handle(blk, t, x, p, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.vapor_phase.mass_transfer_term[t, x, p, j] == 0.0
            else:
                if j in equilibrium_comp:
                    return blk.vapor_phase.mass_transfer_term[t, x, p, j] == \
                        -blk.interphase_mass_transfer[t, x, j]
                else:
                    return blk.vapor_phase.mass_transfer_term[t, x, p, j] == \
                        0.0

        # Heat transfer coefficients
        # Vapor-liquid heat transfer coefficient [J/m2.s.K]

        self.h_v = Var(self.flowsheet().time,
                       self.vapor_phase.length_domain,
                       initialize=100,
                       doc='''Vapor-liquid heat transfer coefficient''')

        # Vapor-liquid heat transfer coeff modified by Ackmann factor [J/m.s.K]
        def rule_heat_transfer_coeff_Ack(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                Ackmann_factor =\
                    sum(blk.vapor_phase.properties[t, x].cp_mol_phase_comp['Vap',j] *
                     blk.interphase_mass_transfer[t, x, j] for j in equilibrium_comp)
                return Ackmann_factor /\
                    (1 - exp(-Ackmann_factor /
                             (blk.h_v[t, x] * blk.area_interfacial[t, x] *
                              blk.area_column)))
        self.h_v_Ack = Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_heat_transfer_coeff_Ack,
            doc='Vap-Liq heat transfer coefficient corrected by Ackmann factor')

        # Heat flux  - vapor side [J/s.m]
        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="heat transfer - vapor side ")
        def vapor_phase_volumetric_heat_flux(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return blk.heat_flux_vap[t, x] == 0
            else:
                zb = self.vapor_phase.length_domain.at(value(self.zi[x]) - 1)
                return blk.heat_flux_vap[t, x] == blk.h_v_Ack[t, x] * \
                    (blk.liquid_phase.properties[t, zb].temperature -
                      blk.vapor_phase.properties[t, x].temperature)
                    
        # Heat transfer - vapor side [J/s.m]
        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="heat transfer - vapor side ")
        def vapor_phase_heat_transfer(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return blk.vapor_phase.heat[t, x] == 0
            else:
                zb = self.vapor_phase.length_domain.at(value(self.zi[x]) - 1)
                return blk.vapor_phase.heat[t, x] == -blk.heat_flux_vap[t, x] - \
                    (sum(blk.vapor_phase.properties[t, x].enth_mol_phase_comp['Vap',j] *
                      blk.vapor_phase.mass_transfer_term[t, x, 'Vap', j] for j in solute_comp_list)) + \
                    (sum(blk.liquid_phase.properties[t, zb].enth_mol_phase_comp['Liq',j] *
                      blk.liquid_phase.mass_transfer_term[t, zb, 'Liq', j] for j in solvent_comp_list))

        # Heat transfer - liquid side [J/s.m]
        @self.Constraint(self.flowsheet().time,
                          self.liquid_phase.length_domain,
                          doc="heat transfer - liquid side ")
        def liquid_phase_heat_transfer(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.liquid_phase.heat[t, x] == 0
            else:
                zf = self.vapor_phase.length_domain.at(value(self.zi[x]) + 1)
                return blk.liquid_phase.heat[t, x] == -blk.vapor_phase.heat[t, zf]

    # =========================================================================
    # Model initialization routine

    def initialize(blk,
                   vapor_phase_state_args=None,
                   liquid_phase_state_args=None,
                   state_vars_fixed=False,
                   outlvl=idaeslog.NOTSET,
                   solver=None,
                   optarg=None):
        """
        Column initialization.

        Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = None).
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

        dynamic_constraints = [
            "pressure_at_interface",
            "mass_transfer_vapor",
            "liquid_phase_mass_transfer_handle",
            "vapor_phase_mass_transfer_handle",
            "vapor_phase_volumetric_heat_flux",
            "vapor_phase_heat_transfer",
            "liquid_phase_heat_transfer"]

        # ---------------------------------------------------------------------
        # Deactivate unit model level constraints (asides geometry constraints)
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in dynamic_constraints:
                c.deactivate()

        # Fix variables

        # Interface pressure
        blk.pressure_equil.fix()

        # Molar flux
        blk.interphase_mass_transfer.fix(0.0)
        blk.vapor_phase.mass_transfer_term.fix(0.0)
        blk.liquid_phase.mass_transfer_term.fix(0.0)

        # Heat transfer rate
        blk.heat_flux_vap.fix(0.0)
        blk.vapor_phase.heat.fix(0.0)
        blk.liquid_phase.heat.fix(0.0)
        
        # # ---------------------------------------------------------------------
        # Provide state arguments for property package initialization

        init_log.info("Step 1: Property Package initialization")
        
        vap_comp = blk.config.vapor_side.property_package.component_list
        liq_apparent_comp = [c[1] for c in blk.liquid_phase.properties.phase_component_set]
        
        if vapor_phase_state_args is None:
            vapor_phase_state_args = {
                'flow_mol': blk.vapor_inlet.flow_mol[0].value,
                'temperature': blk.vapor_inlet.temperature[0].value,
                'pressure': blk.vapor_inlet.pressure[0].value,
                'mole_frac_comp':
                {j: blk.vapor_inlet.mole_frac_comp[0, j].value 
                 for j in vap_comp}}

        if liquid_phase_state_args is None:
            liquid_phase_state_args = {
                'flow_mol': blk.liquid_inlet.flow_mol[0].value,
                'temperature': blk.liquid_inlet.temperature[0].value,
                'pressure': blk.vapor_inlet.pressure[0].value,
                'mole_frac_comp':
                {j: blk.liquid_inlet.mole_frac_comp[0, j].value 
                 for j in liq_apparent_comp}}

        # Initialize vapor_phase properties block
        vflag = blk.vapor_phase.properties.initialize(
            state_args=vapor_phase_state_args,
            state_vars_fixed=False,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)

        # Initialize liquid_phase properties block
        lflag = blk.liquid_phase.properties.initialize(
            state_args=liquid_phase_state_args,
            state_vars_fixed=False,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True)

        init_log.info("Step 2: Steady-State isothermal mass balance")

        blk.vapor_phase.properties.release_state(flags=vflag)

        blk.liquid_phase.properties.release_state(flags=lflag)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 2: {}.".format(idaeslog.condition(res)))

        assert check_optimal_termination(res)

        # ---------------------------------------------------------------------
        init_log.info('Step 3: Interface equilibrium')

        # Activate interface pressure constraint

        blk.pressure_equil.unfix()
        blk.pressure_at_interface.activate()

        # ----------------------------------------------------------------------

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high(
            "Step 3 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------

        init_log.info('Step 4: Isothermal chemical absoption')
        init_log.info_high("No mass transfer to mass transfer")

        # Unfix mass transfer terms
        blk.interphase_mass_transfer.unfix()

        # Activate mass transfer equation in vapor phase
        blk.mass_transfer_vapor.activate()
        
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        
        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()
        blk.vapor_phase_mass_transfer_handle.activate()
        blk.liquid_phase_mass_transfer_handle.activate()
        
        optarg = {
            "tol": 1e-8,
            "max_iter": 150,
            "bound_push":1e-8}
        opt.options = optarg
                    
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
            if res.solver.status != SolverStatus.warning:
                print('')
        init_log.info_high(
            "Step 4 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info('Step 5: Adiabatic chemical absoption')
        init_log.info_high("Isothermal to Adiabatic ")
        
        # Unfix heat transfer terms
        blk.heat_flux_vap.unfix()
        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()

        # Activate heat transfer and steady-state energy balance related equations
        for c in ["vapor_phase_volumetric_heat_flux",
                  "vapor_phase_heat_transfer",
                  "liquid_phase_heat_transfer"]:
            getattr(blk, c).activate()
            
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:    
            res = opt.solve(blk, tee=slc.tee)

        init_log.info_high(
            "Step 5 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------

        if not blk.config.dynamic:
            init_log.info('Steady-state initialization complete')
        
    def fix_initial_condition(blk):
        """
        Initial condition for material and enthalpy balance.

        Mass balance : Initial condition  is determined by
        fixing n-1 mole fraction and the total molar flowrate

        Energy balance :Initial condition  is determined by
        fixing  the temperature.

        """

        vap_comp = blk.config.vapor_side.property_package.component_list
        liq_comp = blk.config.liquid_side.property_package.component_list
        solute_comp_list = blk.config.liquid_side.property_package.solute_set

        for x in blk.vapor_phase.length_domain:
            if x != 0:
                blk.vapor_phase.properties[0, x].temperature.fix()
                blk.vapor_phase.properties[0, x].flow_mol.fix()
            for j in vap_comp:
                if (x != 0 and j not in solute_comp_list):
                    blk.vapor_phase.properties[0, x].mole_frac_comp[j].fix()
        for x in blk.liquid_phase.length_domain:
            if x != 1:
                blk.liquid_phase.properties[0, x].temperature.fix()
                blk.liquid_phase.properties[0, x].flow_mol.fix()
            for j in liq_comp:
                if (x != 1 and j not in solute_comp_list):
                    blk.liquid_phase.properties[0, x].mole_frac_comp[j].fix()

    def unfix_initial_condition(blk):
        """
        Function to unfix initial condition for material and enthalpy balance.

        """

        vap_comp = blk.config.vapor_side.property_package.component_list
        liq_comp = blk.config.liquid_side.property_package.component_list
        solute_comp_list = blk.config.liquid_side.property_package.solute_set

        for x in blk.vapor_phase.length_domain:
            if x != 0:
                blk.vapor_phase.properties[0, x].temperature.unfix()
                blk.vapor_phase.properties[0, x].flow_mol.unfix()
            for j in vap_comp:
                if (x != 0 and j not in solute_comp_list):
                    blk.vapor_phase.properties[0, x].mole_frac_comp[j].unfix()
        for x in blk.liquid_phase.length_domain:
            if x != 1:
                blk.liquid_phase.properties[0, x].temperature.unfix()
                blk.liquid_phase.properties[0, x].flow_mol.unfix()
            for j in liq_comp:
                if (x != 1 and j not in solute_comp_list):
                    blk.liquid_phase.properties[0, x].mole_frac_comp[j].unfix()
                    
    def make_steady_state_column_profile(blk):
        """
        Steady-state Plot function for Temperature and Solute Pressure profile.

        """

        normalised_column_height = [x for x in blk.vapor_phase.length_domain]
        simulation_time = [t for t in blk.flowsheet().time]

        # final time
        tf = simulation_time[-1]
        
        # solute list
        solute_comp_list = blk.config.liquid_side.property_package.solute_set
        solute_profile = []
        
        liquid_temperature_profile = []
        solute_comp_profile = []

        # APPEND RESULTS
        for j in solute_comp_list:
            for x in blk.vapor_phase.length_domain:
                x_liq = blk.liquid_phase.length_domain.at(blk.zi[x].value)
                solute_comp_profile.append(
                    value(1e-3 * blk.vapor_phase.properties[tf, x].pressure *
                          blk.vapor_phase.properties[tf, x].mole_frac_comp[j]))
                liquid_temperature_profile.append(
                    value(blk.liquid_phase.properties[tf, x_liq].temperature))
            solute_profile.append(solute_comp_profile)

        # plot properties
        fontsize = 18
        labelsize = 18
        fig = plt.figure(figsize=(9, 7))
        ax1 = fig.add_subplot(111)
        ax1.set_title('Steady-state column profile',
                      fontsize=16, fontweight='bold')

        # plot primary axis
        lab1 = ax1.plot(normalised_column_height, solute_profile[0],
                        linestyle='--', mec="b", mfc="None",
                        color='b', label='solute partial pressure [kPa]',
                        marker='o')

        ax1.tick_params(axis='y', labelcolor='b',
                        direction='in', labelsize=labelsize)
        ax1.tick_params(axis='x', direction='in', labelsize=labelsize)

        ax1.set_xlabel('Normalise column  height from bottom',
                       fontsize=fontsize)
        ax1.set_ylabel('P_solute  [ kPa]', color='b', fontweight='bold',
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
        labels_1 = [lb.get_label() for lb in lab_1]
        ax1.legend(lab_1, labels_1, loc='lower center', fontsize=fontsize)
        fig.tight_layout()

        # show graph
        plt.show()

    def make_dynamic_column_profile(blk):
        """
        Dynamic Plot function for Temperature and Solute Pressure profile.

        """

        normalised_column_height = [x for x in blk.vapor_phase.length_domain]
        simulation_time = [t for t in blk.flowsheet().time]
        fluegas_flow = [value(blk.vapor_inlet.flow_mol[t])
                        for t in blk.flowsheet().time]

        # final time
        tf = simulation_time[-1]
        nf = len(simulation_time)

        # mid-time
        if nf % 2 == 0:
            tm = int(nf / 2)
        else:
            tm = int(nf / 2 + 1)

        solute_comp_list = blk.config.liquid_side.property_package.solute_set
        solute_profile_mid = []
        solute_profile_fin = []
        liquid_temperature_profile_mid = []
        liquid_temperature_profile_fin = []
        solute_comp_profile_mid = []
        solute_comp_profile_fin = []

        # APPEND RESULTS
        for j in solute_comp_list:
            for x in blk.vapor_phase.length_domain:
                x_liq = blk.liquid_phase.length_domain.at(blk.zi[x].value)
                solute_comp_profile_mid.append(
                    value(1e-3 * blk.vapor_phase.properties[tm, x].pressure *
                          blk.vapor_phase.properties[tm, x].mole_frac_comp[j]))
                solute_comp_profile_fin.append(
                    value(1e-3 * blk.vapor_phase.properties[tf, x].pressure *
                          blk.vapor_phase.properties[tf, x].mole_frac_comp[j]))
    
                liquid_temperature_profile_mid.append(
                    value(blk.liquid_phase.properties[tm, x_liq].temperature))
                liquid_temperature_profile_fin.append(
                    value(blk.liquid_phase.properties[tf, x_liq].temperature))
            solute_profile_mid.append(solute_comp_profile_mid)
            solute_profile_fin.append(solute_comp_profile_fin)

        # plot properties
        fontsize = 18
        labelsize = 18
        fig = plt.figure(figsize=(12, 7))
        ax1 = fig.add_subplot(211)
        ax1.set_title(
            'Column profile @ {0:6.2f} & {1:6.2f} sec'.format(tm, tf),
            fontsize=16, fontweight='bold')

        # plot primary axis
        lab1 = ax1.plot(normalised_column_height, solute_profile_mid[0],
                        linestyle='--', color='b',
                        label='Solute partial pressure [kPa] @ %d' % tm)
        lab2 = ax1.plot(normalised_column_height, solute_profile_fin[0],
                        linestyle='-', color='b',
                        label='Solute partial pressure [kPa] @ %d' % tf)

        ax1.tick_params(axis='y', labelcolor='b',
                        direction='in', labelsize=labelsize)
        ax1.tick_params(axis='x', direction='in', labelsize=labelsize)

        ax1.set_xlabel('Normalise column  height from bottom',
                       fontsize=fontsize)
        ax1.set_ylabel('P_solute  [ kPa]', color='b', fontweight='bold',
                       fontsize=fontsize)

        # plot secondary axis
        ax2 = ax1.twinx()
        lab3 = ax2.plot(
            normalised_column_height,
            liquid_temperature_profile_mid,
            color='g', linestyle='--',
            label='Liquid temperature profile @ {0:6.1f}'.format(tm))
        lab4 = ax2.plot(
            normalised_column_height,
            liquid_temperature_profile_fin,
            color='g', linestyle='-',
            label='Liquid temperature profile @ {0:6.1f}'.format(tf))
        ax2.set_ylabel('T$_{liq}$ [ K ] ', color='g', fontweight='bold',
                       fontsize=fontsize)
        ax2.tick_params(axis='y', labelcolor='g',
                        direction='in', labelsize=labelsize)
        # get the labels
        lab_1 = lab1 + lab2 + lab3 + lab4
        labels_1 = [lb.get_label() for lb in lab_1]
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
