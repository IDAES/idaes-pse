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
IDAES First Generation (GEN1) MEA Rate-Based Packed Column.

Detailed model equations can be found in the supplimentary information
of the paper :
Akula, Paul; Eslick, John; Bhattacharyya, Debangsu; Miller, David
"Model Development, Validation, and Part-Load Optimization of a
MEA-Based Post-Combustion CO2 Capture Process
Under Part-Load and Variable Capture Operation,
Industrial & Engineering Chemistry Research,2021. (submitted)

"""

# Import Python libraries and third-party
import numpy as np
from enum import Enum

# Import Pyomo libraries
from pyomo.environ import (
    Constraint, Expression, Param, Reals, NonNegativeReals,
    value, Var, exp, SolverStatus, units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool

# Import IDAES Libraries
from idaes.core.util.constants import Constants
from idaes.core import (ControlVolume1DBlock, UnitModelBlockData,
                        declare_process_block_class,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        FlowDirection)
from idaes.core.util import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.generic_models.unit_models.heat_exchanger_1D import \
    HeatExchangerFlowPattern as FlowPattern
from idaes.core.control_volume1d import DistributedVars
import idaes.logger as idaeslog


__author__ = "Paul Akula, John Eslick"


# Set up logger
_log = idaeslog.getLogger(__name__)


class ProcessType(Enum):
    absorber = 1
    stripper = 2


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

    CONFIG.declare("flow_type", ConfigValue(
        default=FlowPattern.countercurrent,
        domain=In(FlowPattern),
        description="Flow configuration of PackedColumn",
        doc="""PackedColumn flow pattern,
**default** - FlowPattern.countercurrent.
**Valid values:** {
**FlowPattern.countercurrent** - countercurrent flow,
**FlowPattern.cocurrent** - cocurrent flow}"""))

    # TODO : Consider removing this
    CONFIG.declare("process_type", ConfigValue(
        default=ProcessType.absorber,
        domain=In(ProcessType),
        description="Flag indicating the type of  process",
        doc="""Flag indicating either absorption or stripping process.
**default** - ProcessType.absorber.
**Valid values:** {
**ProcessType.absorber** - absorption process,
**ProcessType.stripper** - stripping process.}"""))

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

    # ==========================================================================

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

        self.length_column = Var(domain=Reals,
                                 initialize=4.9,
                                 units=pyunits.m,
                                 doc='Column length')

        # ---------------------------------------------------------------------
        # Build control volumnes
        # Vapor phase
        # TODO: Add enthalpy transfer terms to energy balances
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
            length_domain_set=self.config.length_domain_set,
            length_var=self.length_column)

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

        # Liquid phase
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
            has_heat_transfer=True)

        # No liquid phase momentum balance due ot assumption of mechanical
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
        # TODO : Remove this
        dcomp = ['CO2', 'H2O']

        # Packing  parameters
        self.voidage_packing = Param(initialize=0.97,
                                     units=pyunits.dimensionless,
                                     mutable=True,
                                     doc="Packing void space m3/m3")

        self.packing_specific_surface = Param(
            initialize=250,
            units=pyunits.m**2/pyunits.m**3,
            mutable=True,
            doc="Packing specific surface area m2/m3")

        self.diameter_hydraulic = Expression(
            expr=4*self.voidage_packing/self.packing_specific_surface,
            doc="Hydraulic diameter")

        # specific constants for volumetric mass transfer coefficients
        # reference:  Billet and Schultes, 1999
        self.Cv_ref = Var(initialize=0.357,
                          doc='''Vapor packing specific constant in
                                Billet and Schultes' (1999) volumetric
                                mass transfer coefficient correlation''')
        self.Cl_ref = Var(initialize=0.5,
                          doc='''Liquid packing specific constant in
                                Billet and Schultes' (1999) volumetric
                                mass transfer coefficient correlation''')
        self.Cv_ref.fix()
        self.Cl_ref.fix()

        # Add object references - others
        R_ref = Constants.gas_constant

        # Unit Model Parameters/sets
        self._zi = Param(self.vapor_phase.length_domain, mutable=True,
                         doc='''Integer indexing parameter required for
                         transfer across boundaries of a given volume element
                         ''')
        # Set the integer  indices
        for i, x in enumerate(self.vapor_phase.length_domain, 1):
            self._zi[x] = i

        # Continuation parameters for initialization
        # TODO: Fix this - these need to be Vars
        self._homotopy_par_m = Param(
            initialize=0,
            mutable=True,
            units=pyunits.dimensionless,
            doc='Mass transfer homotopy parameter')
        self._homotopy_par_h = Param(
            initialize=0,
            mutable=True,
            units=pyunits.dimensionless,
            doc='Heat transfer homotopy parameter')

        # Interfacial area  parameters
        self.area_interfacial_parA = Var(
            initialize=0.6486,
            units=pyunits.dimensionless,
            doc='''Interfacial area parameter A''')

        self.area_interfacial_parB = Var(
            initialize=0.12,
            units=pyunits.dimensionless,
            doc='''Interfacial area parameter B''')
        self.area_interfacial_parA.fix()
        self.area_interfacial_parB.fix()

        # Holdup  parameters
        self.holdup_parA = Var(initialize=24.2355,
                               units=pyunits.dimensionless,
                               doc='''holdup parameter A''')

        self.holdup_parB = Var(initialize=0.6471,
                               units=pyunits.dimensionless,
                               doc='''holdup parameter B''')
        self.holdup_parA.fix()
        self.holdup_parB.fix()

        # Unit Model Variables
        # Geometry
        self.diameter_column = Var(domain=Reals,
                                   initialize=0.1,
                                   units=pyunits.m,
                                   doc='Column diameter')
        self.area_column = Var(domain=Reals,
                               initialize=0.5,
                               units=pyunits.m**2,
                               doc='Column cross-sectional area')

        # Hydrodynamics
        self.velocity_vap = Var(self.flowsheet().time,
                                self.vapor_phase.length_domain,
                                domain=NonNegativeReals,
                                initialize=2,
                                units=pyunits.m / pyunits.s,
                                doc='Vapor superficial velocity')
        self.velocity_liq = Var(self.flowsheet().time,
                                self.liquid_phase.length_domain,
                                units=pyunits.m / pyunits.s,
                                domain=NonNegativeReals,
                                initialize=0.01,
                                doc='Liquid superficial velocity')
        # mass and heat transfer terms
        # mass transfer
        self.pressure_equil = Var(self.flowsheet().time,
                                  self.vapor_phase.length_domain,
                                  dcomp,
                                  domain=NonNegativeReals,
                                  initialize=500,
                                  units=pyunits.Pa,
                                  doc='''Equilibruim pressure of diffusing
                                      components at the interface ''')
        self.N_v = Var(self.flowsheet().time,
                       self.liquid_phase.length_domain,
                       dcomp,
                       domain=Reals,
                       initialize=0.0,
                       units=pyunits.mol / (pyunits.s * pyunits.m),
                       doc='''Moles of diffusing species transfered
                                     into liquid ''')
        self.enhancement_factor = Var(self.flowsheet().time,
                                      self.liquid_phase.length_domain,
                                      units=pyunits.dimensionless,
                                      domain=NonNegativeReals,
                                      initialize=160,
                                      doc='Enhancement factor')

        self.yi_MEA = Var(self.flowsheet().time,
                          self.liquid_phase.length_domain,
                          domain=NonNegativeReals,
                          initialize=0.5,
                          units=pyunits.dimensionless,
                          doc='''Dimensionless concentration of MEA
                                    at interface ''')
        self.yeq_CO2 = Var(self.flowsheet().time,
                           self.liquid_phase.length_domain,
                           domain=NonNegativeReals,
                           initialize=0.5,
                           units=pyunits.dimensionless,
                           doc='''Dimensionless concentration of CO2
                                      in equilibruim with the bulk''')

        # heat transfer
        self.heat_vap = Var(self.flowsheet().time,
                            self.vapor_phase.length_domain,
                            domain=Reals,
                            initialize=0.0,
                            units=pyunits.J / (pyunits.s * pyunits.m),
                            doc='Heat transfer rate in vapor phase')
        self.heat_liq = Var(self.flowsheet().time,
                            self.vapor_phase.length_domain,
                            domain=Reals,
                            initialize=0.0,
                            units=pyunits.J / (pyunits.s * pyunits.m),
                            doc='Heat transfer rate in liquid phase')

        # =====================================================================
        # Add performance equations

        # Inter-facial Area model ([m2/m3]):
        # reference: Tsai correlation,regressed by Chinen et al. 2018
        def rule_interfacial_area(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                return (blk.packing_specific_surface *
                        blk.area_interfacial_parA *
                        (blk.liquid_phase.properties[t, x].dens_mass /
                         blk.liquid_phase.properties[t, x].surf_tens *
                         (blk.velocity_liq[t, x])**(4.0 / 3.0)
                         )**blk.area_interfacial_parB)

        self.area_interfacial = Expression(self.flowsheet().time,
                                           self.vapor_phase.length_domain,
                                           rule=rule_interfacial_area,
                                           doc='Specific inter-facial area')

        # liquid holdup model
        # reference: Tsai correlation,regressed by Chinen et al. 2018
        def rule_holdup_liq(blk, t, x):
            return (blk.holdup_parA *
                    (blk.velocity_liq[t, x] *
                     (blk.liquid_phase.properties[t, x].visc_d /
                      blk.liquid_phase.properties[t, x].dens_mass)**(0.333)
                     )**blk.holdup_parB)

        self.holdup_liq = Expression(self.flowsheet().time,
                                     self.liquid_phase.length_domain,
                                     rule=rule_holdup_liq,
                                     doc='Volumetric liquid holdup [-]')

        # vapor holdup model
        # reference: Tsai correlation,regressed by Chinen et al. 2018
        def rule_holdup_vap(blk, t, x):
            return blk.voidage_packing - blk.holdup_liq[t, x]

        self.holdup_vap = Expression(self.flowsheet().time,
                                     self.vapor_phase.length_domain,
                                     rule=rule_holdup_vap,
                                     doc='Volumetric vapor holdup [-]')

        # ---------------------------------------------------------------------
        # Geometry contraints

        # Column area [m2]
        @self.Constraint(doc="Column cross-sectional area")
        def column_cross_section_area(blk):
            return blk.area_column == (
                Constants.pi * 0.25 * (blk.diameter_column)**2)

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
            def vapor_side_area(blk, t, x):
                return blk.vapor_phase.area[t, x] == (blk.area_column *
                                                      blk.holdup_vap[t, x])

            @self.Constraint(self.flowsheet().time,
                             self.liquid_phase.length_domain,
                             doc=control_volume_area_definition)
            def liquid_side_area(blk, t, x):
                return blk.liquid_phase.area[t, x] == (blk.area_column *
                                                       blk.holdup_liq[t, x])
        else:
            self.vapor_phase.area.fix(value(self.area_column))
            self.liquid_phase.area.fix(value(self.area_column))

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc='''Mechanical equilibruim: vapor-side pressure
                                equal liquid -side pressure''')
        def mechanical_equil(bk, t, x):
            return bk.liquid_phase.properties[t, x].pressure == \
                bk.vapor_phase.properties[t, x].pressure

        # ---------------------------------------------------------------------
        # Hydrodynamic contraints
        # Vapor superficial velocity

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="Vapor superficial velocity")
        def eq_velocity_vap(blk, t, x):
            return blk.velocity_vap[t, x] * blk.area_column * \
                blk.vapor_phase.properties[t, x].conc_mol == \
                blk.vapor_phase.properties[t, x].flow_mol

        # Liquid superficial velocity
        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="Liquid superficial velocity")
        def eq_velocity_liq(blk, t, x):
            return blk.velocity_liq[t, x] * blk.area_column * \
                blk.liquid_phase.properties[t, x].conc_mol == \
                blk.liquid_phase.properties[t, x].flow_mol

        # ---------------------------------------------------------------------
        # Mass transfer coefficients, Billet and Schultes (1999) correlation,
        # where parameters are regressed by Chinen et al. (2018).

        # vapor mass transfer coefficients for diffusing components
        # [mol/m2.s.Pa]
        def rule_mass_transfer_coeff_vap(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                return (
                    1/(R_ref * blk.vapor_phase.properties[t, x].temperature) *
                    blk.Cv_ref/(blk.holdup_vap[t, x])**0.5 *
                    (blk.packing_specific_surface /
                     blk.diameter_hydraulic)**0.5 *
                    (blk.vapor_phase.properties[t, x].diffus[j])**(2/3) *
                    (blk.vapor_phase.properties[t, x].visc_d /
                     blk.vapor_phase.properties[t, x].dens_mass)**(1/3) *
                    ((blk.velocity_vap[t, x] *
                      blk.vapor_phase.properties[t, x].dens_mass) /
                     (blk.packing_specific_surface *
                      blk.vapor_phase.properties[t, x].visc_d))**(3/4))

        self.k_v = Expression(self.flowsheet().time,
                              self.vapor_phase.length_domain,
                              dcomp,
                              rule=rule_mass_transfer_coeff_vap,
                              doc='Vapor mass transfer coefficient')

        # mass transfer coefficients of CO2 in liquid phase  [m/s]
        def rule_mass_transfer_coeff_CO2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return (blk.Cl_ref*12**(1/6) *
                        (blk.velocity_liq[t, x] *
                         blk.liquid_phase.properties[t, x].diffus['CO2'] /
                         (blk.diameter_hydraulic * blk.holdup_liq[t, x]))**0.5)

        self.k_l_CO2 = Expression(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            rule=rule_mass_transfer_coeff_CO2,
            doc='''CO2 mass transfer coefficient in solvent''')

        # mass tranfer terms
        def rule_phi(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                zb = self.vapor_phase.length_domain[self._zi[x].value - 1]
                return (blk.enhancement_factor[t, zb] *
                        blk.k_l_CO2[t, zb] /
                        blk.k_v[t, x, 'CO2'])

        self.phi = Expression(self.flowsheet().time,
                              self.vapor_phase.length_domain,
                              rule=rule_phi,
                              doc='''CO2 Equilibruim partial pressure
                                   intermediate  term''')

        # Equilibruim partial pressure of diffusing components at interface
        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         dcomp,
                         doc='''Equilibruim partial pressure of diffusing
                                components at interface''')
        def pressure_at_interface(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.pressure_equil[t, x, j] == 0.0
            else:
                zb = self.vapor_phase.length_domain[self._zi[x].value - 1]
                if j == 'H2O':
                    return blk.pressure_equil[t, x, j] == (
                        blk.liquid_phase.properties[t, zb].vol_mol *
                        blk.liquid_phase.properties[
                            t, zb].conc_mol_comp_true[j] *
                        blk.liquid_phase.properties[t, zb].pressure_sat[j])
                elif j == 'CO2':
                    return blk.pressure_equil[t, x, j] == (
                        (blk.vapor_phase.properties[t, x].mole_frac_comp[j] *
                         blk.vapor_phase.properties[
                             t, x].pressure + blk.phi[t, x] *
                         blk.liquid_phase.properties[
                             t, zb].conc_mol_comp_true[j]) /
                        (1 + blk.phi[t, x] /
                         blk.liquid_phase.properties[t, zb].henry_N2O_analogy))

        # mass transfer of  diffusing components
        def rule_mass_transfer(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return blk.N_v[t, x, j] == 0.0
            else:
                return blk.N_v[t, x, j] == (blk.k_v[t, x, j] *
                                            blk.area_interfacial[t, x] * blk.area_column *
                                            (blk.vapor_phase.properties[t, x].mole_frac_comp[j] *
                                             blk.vapor_phase.properties[t, x].pressure -
                                             blk.pressure_equil[t, x, j])) * blk._homotopy_par_m

        self.mass_transfer = Constraint(self.flowsheet().time,
                                        self.vapor_phase.length_domain,
                                        dcomp, rule=rule_mass_transfer,
                                        doc="mass transfer to liquid")

        # mass tranfer term handle
        # liquid side

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            self.config.liquid_side.property_package._phase_component_set,
            doc="mass transfer to liquid")
        def liquid_phase_mass_transfer_handle(blk, t, x, p, j):
            if x == self.liquid_phase.length_domain.last():
                return blk.liquid_phase.mass_transfer_term[t, x, p, j] == 0.0
            else:
                zf = self.vapor_phase.length_domain[self._zi[x].value + 1]
                if j == 'MEA':
                    return blk.liquid_phase.mass_transfer_term[t, x, p, j] == \
                        0.0
                else:
                    return blk.liquid_phase.mass_transfer_term[t, x, p, j] == \
                        blk.N_v[t, zf, j]
        # vapor side
        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            self.config.vapor_side.property_package._phase_component_set,
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

        # Heat transfer coefficients, Chilton Colburn  analogy
        # Vapor-liquid heat transfer coefficient [J/m2.s.K]

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

        self.h_v = Expression(self.flowsheet().time,
                              self.vapor_phase.length_domain,
                              rule=rule_heat_transfer_coeff,
                              doc='''vap-liq heat transfer coefficient''')

        # Vapor-liquid heat transfer coefficient modified by Ackmann factor [J/m.s.K]
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
                             (blk.h_v[t, x] * blk.area_interfacial[t, x] * blk.area_column)))
        self.h_v_Ack = Expression(self.flowsheet().time,
                                  self.vapor_phase.length_domain,
                                  rule=rule_heat_transfer_coeff_Ack,
                                  doc='''vap-liq heat transfer coefficient corrected
                                         by Ackmann factor''')

        # heat transfer vapor  side [J/s.m]
        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="heat transfer - vapor side ")
        def vapor_phase_heat_transfer(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return blk.heat_vap[t, x] == 0
            else:
                zb = self.vapor_phase.length_domain[value(self._zi[x]) - 1]
                return blk.heat_vap[t, x] == blk.h_v_Ack[t, x] * \
                    (blk.liquid_phase.properties[t, zb].temperature -
                     blk.vapor_phase.properties[t, x].temperature) * \
                    blk._homotopy_par_h

        # heat transfer liquid side [J/s.m]
        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="heat transfer - liquid side ")
        def liquid_phase_heat_transfer(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.heat_liq[t, x] == 0
            else:
                zf = self.vapor_phase.length_domain[value(self._zi[x]) + 1]
                return blk.heat_liq[t, x] == blk.heat_vap[t, zf] + \
                    (blk.liquid_phase.properties[t, x].habs * blk.N_v[t, zf, 'CO2'] -
                     blk.liquid_phase.properties[t, x].hvap * blk.N_v[t, zf, 'H2O']) *\
                    blk._homotopy_par_h

        # heat transfer handle
        # vapor  heat transfer handle

        @self.Constraint(self.flowsheet().time,
                         self.vapor_phase.length_domain,
                         doc="vapor - heat transfer handle")
        def vapor_phase_heat_transfer_handle(blk, t, x):
            return blk.vapor_phase.heat[t, x] == blk.heat_vap[t, x]

        # liquid  heat transfer handle
        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc="liquid - heat transfer handle")
        def liquid_phase_heat_transfer_handle(blk, t, x):
            return blk.liquid_phase.heat[t, x] == -blk.heat_liq[t, x]

        # Enhancement factor model
        # reference: Jozsef Gaspar,Philip Loldrup Fosbol, (2015)
        # self.yi_MEA[z] is equivalent to sqrt(yi_MEA) in the document

        def rule_conc_mol_comp_interface_CO2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                zf = self.liquid_phase.length_domain[self._zi[x].value + 1]
                return blk.pressure_equil[t, zf, 'CO2'] /\
                    blk.liquid_phase.properties[t, x].henry_N2O_analogy

        self.conc_mol_comp_CO2_eq = Expression(self.flowsheet().time,
                                               self.liquid_phase.length_domain,
                                               rule=rule_conc_mol_comp_interface_CO2,
                                               doc='''Concentration of CO2
                                                at the interface ]''')

        def rule_Hatta(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return (blk.liquid_phase.properties[t, x].k2_rxn *
                        blk.liquid_phase.properties[t, x].conc_mol_comp_true['MEA'] *
                        blk.liquid_phase.properties[t, x].diffus['CO2'])**0.5 /\
                    blk.k_l_CO2[t, x]

        self.Hatta = Expression(self.flowsheet().time,
                                self.liquid_phase.length_domain,
                                rule=rule_Hatta,
                                doc='Hatta number')

        def rule_yb_CO2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Expression.Skip
            else:
                return blk.liquid_phase.properties[t, x].conc_mol_comp_true['CO2'] /\
                    blk.conc_mol_comp_CO2_eq[t, x]

        self.yb_CO2 = Expression(self.flowsheet().time,
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

        self.instant_E = Expression(self.flowsheet().time,
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

        self.yi_MEACOO = Expression(self.flowsheet().time,
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

        self.yi_MEAH = Expression(self.flowsheet().time,
                                  self.liquid_phase.length_domain,
                                  rule=rule_yi_MEAH,
                                  doc='Dimensionless concentration of MEA+')

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc='''dimensionless concentration of CO2
                                at equilibruim with the bulk ''')
        def yeq_CO2_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.yeq_CO2[t, x] == 0.0
            else:
                return blk.yeq_CO2[t, x] * blk.yi_MEA[t, x]**4 == \
                    blk.yb_CO2[t, x] * blk.yi_MEAH[t, x] * blk.yi_MEACOO[t, x]

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc='''Enhancement factor model Eqn 1 ''')
        def E1_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.enhancement_factor[t, x] == 1
            else:
                return (blk.enhancement_factor[t, x] - 1) * (1 - blk.yb_CO2[t, x]) == \
                    (blk.instant_E[t, x] - 1) * (1 - blk.yi_MEA[t, x]**2)

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc='''Enhancement factor model Eqn 2 ''')
        def E2_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.yi_MEA[t, x] == 0
            else:
                return blk.enhancement_factor[t, x] * (1 - blk.yb_CO2[t, x]) == \
                    blk.Hatta[t, x] * blk.yi_MEA[t, x] * \
                    (1 - blk.yeq_CO2[t, x])

        @self.Constraint(self.flowsheet().time,
                         self.liquid_phase.length_domain,
                         doc='Enhancement factor lower bound ')
        def E3_eqn(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return 1 - blk.enhancement_factor[t, x] <= 0.0

        if self.config.dynamic:
            self.fix_initial_condition()

    # ==========================================================================
    # Model initialization routine

    def initialize(blk,
                   vapor_phase_state_args=None,
                   liquid_phase_state_args=None,
                   state_vars_fixed=False,
                   homotopy_steps_m=None,
                   homotopy_steps_h=None,
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
            homotopy_steps_m : List of continuations steps between 0 and 1
                               for turning mass transfer constrainst gradually
            homotopy_steps_h : List of continuations steps between 0 and 1
                               for turning heat transfer constraints gradually
            optarg : solver options dictionary object (default=None, use
                     default solver options)
            solver : str indicating which solver to use during initialization
                    (default = None, use IDAES default solver)

        """

        # Set up logger for initialization and solve
        init_log = idaeslog.getInitLogger(blk.name, outlvl, tag="unit")
        solve_log = idaeslog.getSolveLogger(blk.name, outlvl, tag="unit")

        # Set solver options
        # TODO: Work out why using default solver here doubles test run time
        opt = get_solver(solver, optarg)

        if homotopy_steps_m is None:
            homotopy_steps_m = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1]

        if homotopy_steps_h is None:
            homotopy_steps_h = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

        dynamic_constraints = [
            "vapor_side_area",
            "liquid_side_area",
            "eq_velocity_vap",
            "eq_velocity_liq",
            "E1_eqn",
            "E2_eqn",
            "pressure_at_interface",
            "mass_transfer",
            "liquid_phase_mass_transfer_handle",
            "vapor_phase_mass_transfer_handle",
            "vapor_phase_heat_transfer",
            "liquid_phase_heat_transfer",
            "vapor_phase_heat_transfer_handle",
            "liquid_phase_heat_transfer_handle",
            "yeq_CO2_eqn",
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
            "mechanical_equil",
            "pressure_balance"]

        # ---------------------------------------------------------------------
        # Deactivate unit model level constraints (asides geometry constraints)
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in dynamic_constraints:
                c.deactivate()

        # Fix some variables
        # Hydrodynamics - velocity
        blk.velocity_liq.fix()
        blk.velocity_vap.fix()

        # interface pressure
        blk.pressure_equil.fix()

        # flux
        blk.N_v.fix(0.0)
        blk.vapor_phase.mass_transfer_term.fix(0.0)
        blk.liquid_phase.mass_transfer_term.fix(0.0)

        # Enhancement factor model
        blk.enhancement_factor.fix()
        blk.yi_MEA.fix()
        blk.yeq_CO2.fix()

        # heat transfer
        blk.heat_vap.fix(0.0)
        blk.heat_liq.fix(0.0)
        blk.vapor_phase.heat.fix(0.0)
        blk.liquid_phase.heat.fix(0.0)
        # area
        blk.vapor_phase.area.fix(value(blk.area_column))
        blk.liquid_phase.area.fix(value(blk.area_column))

        # other variables
        # Pressure_dx
        blk.vapor_phase.pressure_dx[:, :].fix(0.0)

        # vapor side flow terms
        blk.vapor_phase._enthalpy_flow.fix(1.0)
        blk.vapor_phase.enthalpy_flow_dx[:, :, :].fix(0.0)
        blk.vapor_phase._flow_terms.fix(1.0)
        blk.vapor_phase.material_flow_dx[:, :, :, :].fix(0.0)

        # liquid side flow terms
        blk.liquid_phase._enthalpy_flow.fix(1.0)
        blk.liquid_phase.enthalpy_flow_dx[:, :, :].fix(0.0)
        blk.liquid_phase._flow_terms.fix(1.0)
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
            blk.unfix_initial_condition()

        # ---------------------------------------------------------------------
        # get values for state variables for initialization
        if vapor_phase_state_args is None:
            if blk.config.process_type == ProcessType.absorber:
                vapor_phase_state_args = {
                    'flow_mol': blk.vapor_inlet.flow_mol[0].value,
                    'temperature': blk.vapor_inlet.temperature[0].value,
                    'pressure': blk.vapor_inlet.pressure[0].value,
                    'mole_frac_comp':
                    {'H2O': blk.vapor_inlet.mole_frac_comp[0, 'H2O'].value,
                     'CO2': blk.vapor_inlet.mole_frac_comp[0, 'CO2'].value,
                     'N2': blk.vapor_inlet.mole_frac_comp[0, 'N2'].value,
                     'O2': blk.vapor_inlet.mole_frac_comp[0, 'O2'].value}}
            elif blk.config.process_type == ProcessType.stripper:
                vapor_phase_state_args = {
                    'flow_mol': blk.vapor_inlet.flow_mol[0].value,
                    'temperature': blk.vapor_inlet.temperature[0].value,
                    'pressure': blk.vapor_inlet.pressure[0].value,
                    'mole_frac_comp':
                    {'H2O': blk.vapor_inlet.mole_frac_comp[0, 'H2O'].value,
                     'CO2': blk.vapor_inlet.mole_frac_comp[0, 'CO2'].value}}

        if liquid_phase_state_args is None:
            liquid_phase_state_args = {
                'flow_mol': blk.liquid_inlet.flow_mol[0].value,
                'temperature': blk.liquid_inlet.temperature[0].value,
                'pressure': blk.vapor_inlet.pressure[0].value,
                'mole_frac_comp':
                {'H2O': blk.liquid_inlet.mole_frac_comp[0, 'H2O'].value,
                 'CO2': blk.liquid_inlet.mole_frac_comp[0, 'CO2'].value,
                 'MEA': blk.liquid_inlet.mole_frac_comp[0, 'MEA'].value}}

        init_log.info("STEP 1: Property Package initialization")

        # Initialize vapor_phase block
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

        init_log.info("STEP 2: Steady-State ISOTHERMAL MASS BALANCE")
        init_log.info_high('No mass transfer ')
        init_log.info_high('No heat transfer')

        # unfix flow variable terms
        # vapor side

        blk.vapor_phase.properties.release_state(flags=vflag)
        blk.vapor_phase.properties[:, :].temperature.fix()
        blk.vapor_phase.properties[:, :].pressure.fix()

        blk.vapor_phase._flow_terms[:, :, :, :].unfix()
        blk.vapor_phase.material_flow_dx[:, :, :, :].unfix()
        # liquid-side
        blk.liquid_phase.properties.release_state(flags=lflag)
        blk.liquid_phase.properties[:, :].temperature.fix()
        blk.liquid_phase.properties[:, :].pressure.fix()

        blk.liquid_phase._flow_terms[:, :, :, :].unfix()
        blk.liquid_phase.material_flow_dx[:, :, :, :].unfix()

        # activate mass balance related equations
        # liquid control volume

        for c in [
                "material_balances",
                "material_flow_linking_constraints",
                "material_flow_dx_disc_eq"]:
            getattr(blk.liquid_phase, c).activate()
        # vapor control volume
        for c in [
                "material_balances",
                "material_flow_linking_constraints",
                "material_flow_dx_disc_eq"]:
            getattr(blk.vapor_phase, c).activate()

        # solve for a small length if stripper
        if (blk.config.process_type == ProcessType.stripper):
            _specified_length = value(blk.length_column)
            blk.length_column.fix(0.6)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 2: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info('STEP 3: Add Mass tranfer terms')
        init_log.info_high('(3a) Velocities & Interface pressure')
        init_log.info_high('(3b) Enhancement factor')

        # Initialize : Velocities, Interface pressure, Enhancement factor

        # velocity
        blk.velocity_vap.unfix()
        blk.velocity_liq.unfix()
        blk.eq_velocity_vap.activate()
        blk.eq_velocity_liq.activate()
        for t in blk.flowsheet().time:
            for x in blk.vapor_phase.length_domain:
                blk.velocity_vap[t, x].value = value(
                    blk.vapor_phase.properties[t, x].flow_mol /
                    (blk.area_column * blk.vapor_phase.properties[t, x].conc_mol))
            for x in blk.liquid_phase.length_domain:
                blk.velocity_liq[t, x].value = value(
                    blk.liquid_phase.properties[t, x].flow_mol /
                    (blk.area_column * blk.liquid_phase.properties[t, x].conc_mol))

        # Interface pressure
        blk.pressure_equil.unfix()
        blk.pressure_at_interface.activate()
        blk.enhancement_factor.fix(1)

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 3a: {}.".format(idaeslog.condition(res)))
        # ----------------------------------------------------------------------
        # Enhancement factor model
        blk.enhancement_factor.unfix()
        for t in blk.flowsheet().time:
            for x in blk.liquid_phase.length_domain:
                blk.enhancement_factor[t, x].value = 100
        blk.yi_MEA.unfix()
        blk.yeq_CO2.unfix()
        blk.E1_eqn.activate()
        blk.E2_eqn.activate()
        blk.yeq_CO2_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 3 complete: {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------

        init_log.info('STEP 4: Isothermal chemical absoption')
        init_log.info_high("Homotopy steps: ")
        init_log.info_high("No mass transfer (0.0) --> (1.0) mass transfer")

        # ISOTHERMAL CHEMICAL ABSORPTION
        blk.N_v.unfix()
        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()
        blk.mass_transfer.activate()
        blk.vapor_phase_mass_transfer_handle.activate()
        blk.liquid_phase_mass_transfer_handle.activate()

        for i in homotopy_steps_m:
            init_log.info('homotopy step -->{0:5.2f}'.format(i))
            blk._homotopy_par_m = i
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
                if res.solver.status != SolverStatus.warning:
                    print('')

        init_log.info_high("Step 4 complete: {}.".format(idaeslog.condition(res)))
        # ---------------------------------------------------------------------
        for c in ["mechanical_equil"]:
            getattr(blk, c).activate()
        for c in ["pressure_balance", "pressure_dx_disc_eq"]:
            getattr(blk.vapor_phase, c).activate()

        blk.vapor_phase.pressure_dx[:, :].unfix()

        # Unfix pressure
        for t in blk.flowsheet().time:
            for x in blk.vapor_phase.length_domain:
                # Unfix all vapor pressure variables except at the inlet
                if (blk.vapor_phase.properties[t, x].config.defined_state
                        is False):
                    blk.vapor_phase.properties[t, x].pressure.unfix()
            for x in blk.liquid_phase.length_domain:
                    blk.liquid_phase.properties[t, x].pressure.unfix()

        # ---------------------------------------------------------------------
        init_log.info('STEP 5: Adiabatic chemical absoption')
        init_log.info_high("Homotopy steps:")
        init_log.info_high("Isothermal (0.0) --> (1.0) Adiabatic ")

        # Unfix temperature
        for t in blk.flowsheet().time:
            for x in blk.vapor_phase.length_domain:
                # Unfix all vapor temperature variables except at the inlet
                if (blk.vapor_phase.properties[t, x].config.defined_state
                        is False):
                    blk.vapor_phase.properties[t, x].temperature.unfix()
            for x in blk.liquid_phase.length_domain:
                # Unfix all liquid temperature variables except at the inlet
                if (blk.liquid_phase.properties[t, x].config.defined_state
                        is False):
                    blk.liquid_phase.properties[t, x].temperature.unfix()

        # unfix heat transfer terms
        blk.heat_vap.unfix()
        blk.heat_liq.unfix()
        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()

        # unfix enthalpy flow variable terms
        blk.vapor_phase._enthalpy_flow[:, :, :].unfix()
        blk.vapor_phase.enthalpy_flow_dx[:, :, :].unfix()
        blk.liquid_phase._enthalpy_flow[:, :, :].unfix()
        blk.liquid_phase.enthalpy_flow_dx[:, :, :].unfix()

        # activate steady-state energy balance related equations
        # unit model
        for c in [
                "vapor_phase_heat_transfer",
                "liquid_phase_heat_transfer",
                "vapor_phase_heat_transfer_handle",
                "liquid_phase_heat_transfer_handle"]:
            getattr(blk, c).activate()

        # liquid control volume
        for c in [
                "enthalpy_flow_linking_constraint",
                "enthalpy_flow_dx_disc_eq",
                "enthalpy_balances"]:
            getattr(blk.liquid_phase, c).activate()

        # vapor control volume
        for c in [
                "enthalpy_flow_linking_constraint",
                "enthalpy_flow_dx_disc_eq",
                "enthalpy_balances"]:
            getattr(blk.vapor_phase, c).activate()

        for i in homotopy_steps_h:
            init_log.info('homotopy step -->{0:5.2f}'.format(i))
            blk._homotopy_par_h = i
            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)

        init_log.info_high("Step 5 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        # scale up at this if stripper
        if (blk.config.process_type == ProcessType.stripper):
            packing_height = np.linspace(0.6, _specified_length, num=10)
            init_log.info_high('SCALEUP Stripper height')
            for L in packing_height:
                blk.length_column.fix(L)
                print('Packing height = {:6.2f}'.format(L))
                with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                    res = opt.solve(blk, tee=slc.tee)
                init_log.info_high("Scaleup: {}.".format(idaeslog.condition(res)))

        if not blk.config.dynamic:
            init_log.info('STEADY-STATE INITIALIZATION COMPLETE')

        if blk.config.dynamic:
            init_log.info('STEP 6: unfix Accumulation and Holdup terms')
            init_log.info_high("6a Holdup calculations")
            init_log.info_high("6b Include Accumulation terms")

            # activate holdup constraints
            # unit model
            for c in [
                    "vapor_side_area",
                    "liquid_side_area"]:
                getattr(blk, c).activate()
            # liquid control volume
            for c in [
                    "material_holdup_calculation",
                    "energy_holdup_calculation"]:
                getattr(blk.liquid_phase, c).activate()
            # vapor control volume
            for c in [
                    "material_holdup_calculation",
                    "energy_holdup_calculation"]:
                getattr(blk.vapor_phase, c).activate()

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
            init_log.info_high("Step 6a complete: {}.".format(
                idaeslog.condition(res)))

            # Step 6b:
            # unfix accumulation terms(derivative variables)
            blk.vapor_phase.energy_accumulation[:, :, :].unfix()
            blk.vapor_phase.material_accumulation[:, :, :, :].unfix()
            blk.liquid_phase.energy_accumulation[:, :, :].unfix()
            blk.liquid_phase.material_accumulation[:, :, :, :].unfix()

            # activate constraints for accumulation terms
            # liquid control volume
            for c in [
                    "material_accumulation_disc_eq",
                    "energy_accumulation_disc_eq"]:
                getattr(blk.liquid_phase, c).activate()

            # vapor control volume
            for c in [
                    "material_accumulation_disc_eq",
                    "energy_accumulation_disc_eq"]:
                getattr(blk.vapor_phase, c).activate()

            blk.fix_initial_condition()

            with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
                res = opt.solve(blk, tee=slc.tee)
            init_log.info_high("Step 6 complete: {}.".format(
                idaeslog.condition(res)))
            init_log.info('INITIALIZATION COMPLETE')

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

        for x in blk.vapor_phase.length_domain:
            if x != 0:
                blk.vapor_phase.properties[0, x].temperature.fix()
                blk.vapor_phase.properties[0, x].flow_mol.fix()
            for j in vap_comp:
                if (x != 0 and j != 'CO2'):
                    blk.vapor_phase.properties[0, x].mole_frac_comp[j].fix()
        for x in blk.liquid_phase.length_domain:
            if x != 1:
                blk.liquid_phase.properties[0, x].temperature.fix()
                blk.liquid_phase.properties[0, x].flow_mol.fix()
            for j in liq_comp:
                if (x != 1 and j != 'CO2'):
                    blk.liquid_phase.properties[0, x].mole_frac_comp[j].fix()

    def unfix_initial_condition(blk):
        """
        Function to unfix initial condition for material and enthalpy balance.

        """

        vap_comp = blk.config.vapor_side.property_package.component_list
        liq_comp = blk.config.liquid_side.property_package.component_list

        for x in blk.vapor_phase.length_domain:
            if x != 0:
                blk.vapor_phase.properties[0, x].temperature.unfix()
                blk.vapor_phase.properties[0, x].flow_mol.unfix()
            for j in vap_comp:
                if (x != 0 and j != 'CO2'):
                    blk.vapor_phase.properties[0, x].mole_frac_comp[j].unfix()
        for x in blk.liquid_phase.length_domain:
            if x != 1:
                blk.liquid_phase.properties[0, x].temperature.unfix()
                blk.liquid_phase.properties[0, x].flow_mol.unfix()
            for j in liq_comp:
                if (x != 1 and j != 'CO2'):
                    blk.liquid_phase.properties[0, x].mole_frac_comp[j].unfix()
