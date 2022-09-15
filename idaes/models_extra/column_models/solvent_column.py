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
IDAES General Purpose Packed Solvent Column Model.
"""

# Import Pyomo libraries
from pyomo.environ import (
    Constraint,
    Expression,
    Param,
    Reals,
    NonNegativeReals,
    Var,
    SolverStatus,
    TerminationCondition,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock, ConfigValue, Bool, In
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES Libraries
from idaes.core.util.constants import Constants
from idaes.core import (
    ControlVolume1DBlock,
    UnitModelBlockData,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    FlowDirection,
    MaterialFlowBasis,
    useDefault,
    DistributedVars,
)
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError, InitializationError
import idaes.logger as idaeslog
from idaes.core.util.config import DefaultBool


__author__ = "Paul Akula, John Eslick, Anuja Deshpande, Andrew Lee"


# Set up logger
_log = idaeslog.getLogger(__name__)


@declare_process_block_class("PackedColumn")
class PackedColumnData(UnitModelBlockData):
    """
    Packed Column Model Class.
    """

    CONFIG = ConfigBlock()
    # For now, only support steady-state
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Solvent Columns do not yet support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            description="Holdup construction flag",
            doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - construct holdup terms,
**False** - do not construct holdup terms}""",
        ),
    )

    # Configuration template for phase specific  arguments
    _PhaseCONFIG = ConfigBlock()

    CONFIG.declare(
        "finite_elements",
        ConfigValue(
            default=20,
            domain=int,
            description="Number of finite elements length domain",
            doc="""Number of finite elements to use when discretizing length
domain (default=20)""",
        ),
    )

    CONFIG.declare(
        "length_domain_set",
        ConfigValue(
            default=[0.0, 1.0],
            domain=list,
            description="List of points in length domain",
            doc="""length_domain_set - (optional) list of point to use to
initialize a new ContinuousSet if length_domain is not
provided (default = [0.0, 1.0])""",
        ),
    )

    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed, **default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )

    # Populate the phase side template to default values
    _PhaseCONFIG.declare(
        "property_package",
        ConfigValue(
            default=None,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a ParameterBlock object""",
        ),
    )

    _PhaseCONFIG.declare(
        "property_package_args",
        ConfigValue(
            default={},
            description="Arguments for constructing vapor property package",
            doc="""A dict of arguments to be passed to the PropertyBlockData
and used when constructing these
(default = 'use_parent_value')
- 'use_parent_value' - get package from parent (default = None)
- a dict (see property package for documentation)
            """,
        ),
    )

    # Create individual config blocks for vapor(gas) and liquid sides
    CONFIG.declare("vapor_phase", _PhaseCONFIG(doc="vapor side config arguments"))

    CONFIG.declare("liquid_phase", _PhaseCONFIG(doc="liquid side config arguments"))

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

        # Check phase lists match assumptions
        if self.config.vapor_phase.property_package.phase_list != ["Vap"]:
            raise ConfigurationError(
                f"{self.name} SolventCondenser model requires that the vapor "
                f"phase property package have a single phase named 'Vap'"
            )
        if self.config.liquid_phase.property_package.phase_list != ["Liq"]:
            raise ConfigurationError(
                f"{self.name} SolventCondenser model requires that the liquid "
                f"phase property package have a single phase named 'Liq'"
            )

        # Check for at least one common component in component lists
        if not any(
            j in self.config.vapor_phase.property_package.component_list
            for j in self.config.liquid_phase.property_package.component_list
        ):
            raise ConfigurationError(
                f"{self.name} SolventCondenser model requires that the liquid "
                f"and vapor phase property packages have at least one "
                f"common component."
            )

        # Geometry
        self.diameter_column = Var(
            domain=Reals, initialize=0.1, units=pyunits.m, doc="Column diameter"
        )

        self.area_column = Var(
            domain=Reals,
            initialize=0.5,
            units=pyunits.m**2,
            doc="Column cross-sectional area",
        )

        self.length_column = Var(
            domain=Reals, initialize=4.9, units=pyunits.m, doc="Column length"
        )

        @self.Constraint(doc="Column cross-sectional area")
        def column_cross_section_area_eqn(blk):
            return blk.area_column == (Constants.pi * 0.25 * (blk.diameter_column) ** 2)

        # =====================================================================
        """ Set argument values for vapor and liquid sides"""

        # Set flow directions for the control volume blocks
        # Gas flows from 0 to 1, Liquid flows from 1 to 0

        # TODO : Add support for different discretization schemes in future
        # TODO : Add support for dynamics in future (if required)
        # TODO: Only handling countercurrent flow for now.
        set_direction_vapor = FlowDirection.forward
        set_direction_liquid = FlowDirection.backward

        # =====================================================================
        """ Build Control volume 1D for vapor phase and
            populate vapor control volume"""

        self.vapor_phase = ControlVolume1DBlock(
            transformation_method="dae.finite_difference",
            transformation_scheme="BACKWARD",
            finite_elements=self.config.finite_elements,
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            area_definition=DistributedVars.variant,
            property_package=self.config.vapor_phase.property_package,
            property_package_args=self.config.vapor_phase.property_package_args,
        )

        self.vapor_phase.add_geometry(
            flow_direction=set_direction_vapor,
            length_domain_set=self.config.length_domain_set,
            length_var=self.length_column,
        )

        self.vapor_phase.add_state_blocks(
            information_flow=set_direction_vapor, has_phase_equilibrium=False
        )

        self.vapor_phase.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_phase_equilibrium=False,
            has_mass_transfer=True,
        )

        self.vapor_phase.add_energy_balances(
            balance_type=EnergyBalanceType.enthalpyTotal,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        self.vapor_phase.add_momentum_balances(
            balance_type=MomentumBalanceType.pressureTotal,
            has_pressure_change=self.config.has_pressure_change,
        )

        self.vapor_phase.apply_transformation()

        # =====================================================================
        """ Build Control volume 1D for liquid phase and
            populate liquid control volume

        """
        self.liquid_phase = ControlVolume1DBlock(
            transformation_method="dae.finite_difference",
            transformation_scheme="FORWARD",
            finite_elements=self.config.finite_elements,
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            area_definition=DistributedVars.variant,
            property_package=self.config.liquid_phase.property_package,
            property_package_args=self.config.liquid_phase.property_package_args,
        )

        self.liquid_phase.add_geometry(
            flow_direction=set_direction_liquid,
            length_domain_set=self.config.length_domain_set,
            length_var=self.length_column,
        )

        self.liquid_phase.add_state_blocks(
            information_flow=set_direction_liquid, has_phase_equilibrium=False
        )

        self.liquid_phase.add_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            has_phase_equilibrium=False,
            has_mass_transfer=True,
        )

        self.liquid_phase.add_energy_balances(
            balance_type=EnergyBalanceType.enthalpyTotal,
            has_heat_transfer=True,
            has_enthalpy_transfer=True,
        )

        self.liquid_phase.apply_transformation()

        # Add Ports for vapor side
        self.add_inlet_port(name="vapor_inlet", block=self.vapor_phase)
        self.add_outlet_port(name="vapor_outlet", block=self.vapor_phase)

        # Add Ports for liquid side
        self.add_inlet_port(name="liquid_inlet", block=self.liquid_phase)
        self.add_outlet_port(name="liquid_outlet", block=self.liquid_phase)

        # TODO : Fix units

        # ======================================================================
        # Unit level sets
        vap_comp = self.config.vapor_phase.property_package.component_list
        liq_comp = self.config.liquid_phase.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp

        # Get units of measurement - liquid phase units will be used as basis
        t_init = self.flowsheet().time.first()
        if (
            self.liquid_phase.properties[t_init, 0].get_material_flow_basis()
            != self.vapor_phase.properties[t_init, 0].get_material_flow_basis()
        ):
            raise ConfigurationError(
                f"{self.name} vapor and liquid property packages must use the "
                f"same material flow basis."
            )

        vunits = (
            self.config.vapor_phase.property_package.get_metadata().get_derived_units
        )
        lunits = (
            self.config.liquid_phase.property_package.get_metadata().get_derived_units
        )
        flow_basis = self.liquid_phase.properties[t_init, 0].get_material_flow_basis()
        if flow_basis == MaterialFlowBasis.molar:
            fb = "flow_mole"
        elif flow_basis == MaterialFlowBasis.molar:
            fb = "flow_mass"
        else:
            raise ConfigurationError(
                f"{self.name} SolventCondenser only supports mass or molar "
                f"basis for MaterialFlowBasis."
            )

        # Hydrodynamics and cacking parameters
        self.eps_ref = Param(
            initialize=0.97,
            units=pyunits.dimensionless,
            mutable=True,
            doc="Packing void space m3/m3",
        )

        self.packing_specific_area = Param(
            initialize=250,
            units=lunits("length") ** 2 / lunits("length") ** 3,
            mutable=True,
            doc="Packing specific surface area",
        )

        self.packing_channel_size = Param(
            initialize=0.1,
            units=lunits("length"),
            mutable=True,
            doc="Packing channel size",
        )

        self.hydraulic_diameter = Expression(
            expr=4 * self.eps_ref / self.packing_specific_area, doc="Hydraulic diameter"
        )

        self.area_interfacial = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=0.9,
            units=(pyunits.m) ** 2 / (pyunits.m) ** 3,
            doc="Specific interfacial area",
        )

        # Liquid and vapor holdups
        self.holdup_liq = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            initialize=0.001,
            units=pyunits.dimensionless,
            doc="Volumetric liquid holdup [-]",
        )

        # TODO : Consider making this a Var & Constraint
        def rule_holdup_vap(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Expression.Skip
            else:
                zb = self.vapor_phase.length_domain.prev(x)
                return blk.eps_ref - blk.holdup_liq[t, zb]

        self.holdup_vap = Expression(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            rule=rule_holdup_vap,
            doc="Volumetric vapor holdup [-]",
        )

        # Area of control volume : vapor side and liquid side
        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Vapor phase cross-sectional area constraint",
        )
        def vapor_phase_area(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return blk.vapor_phase.area[t, x] == (blk.eps_ref * blk.area_column)
            else:
                return blk.vapor_phase.area[t, x] == (
                    blk.area_column * blk.holdup_vap[t, x]
                )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Liquid phase cross-sectional area constraint",
        )
        def liquid_phase_area(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return blk.liquid_phase.area[t, x] == (blk.eps_ref * blk.area_column)
            else:
                return blk.liquid_phase.area[t, x] == (
                    blk.area_column * blk.holdup_liq[t, x]
                )

        # =====================================================================
        # Add performance equations
        # Pressure equality in phases
        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Mechanical equilibruim constraint",
        )
        def mechanical_equilibrium(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                return blk.liquid_phase.properties[t, x].pressure == pyunits.convert(
                    blk.vapor_phase.properties[t, x].pressure,
                    to_units=lunits("pressure"),
                )

        # ---------------------------------------------------------------------
        # Mass transfer relationships
        self.mass_transfer_coeff_vap = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            units=(
                lunits("amount")
                / lunits("pressure")
                / lunits("length") ** 2
                / lunits("time")
            ),
            doc="Vapor phase mass transfer coefficient",
        )

        # Equilibruim partial pressure of components at interface
        self.pressure_equil = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            domain=NonNegativeReals,
            initialize=500,
            units=lunits("pressure"),
            doc="Equilibruim pressure of components at interface",
        )

        # Mass transfer constraints
        self.interphase_mass_transfer = Var(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            equilibrium_comp,
            domain=Reals,
            initialize=0.1,
            units=lunits("amount") / (lunits("time") * lunits("length")),
            doc="Interphase mass transfer rate",
        )

        self.liquid_phase_mass_transfer_model()

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="Interphase mass transfer constraint",
        )
        def interphase_mass_transfer_eqn(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            elif j in equilibrium_comp:
                return blk.interphase_mass_transfer[t, x, j] == (
                    blk.mass_transfer_coeff_vap[t, x, j]
                    * blk.area_interfacial[t, x]
                    * blk.area_column
                    * (
                        blk.vapor_phase.properties[t, x].mole_frac_comp[j]
                        * pyunits.convert(
                            blk.vapor_phase.properties[t, x].pressure,
                            to_units=lunits("pressure"),
                        )
                        - blk.pressure_equil[t, x, j]
                    )
                )
            else:
                return blk.interphase_mass_transfer[t, x, j] == 0.0

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            self.liquid_phase.properties.component_list,
            doc="Liquid phase mass transfer constraint",
        )
        def liquid_mass_transfer_eqn(blk, t, x, j):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            elif j in equilibrium_comp:
                zf = self.vapor_phase.length_domain.next(x)
                return (
                    blk.liquid_phase.mass_transfer_term[t, x, "Liq", j]
                    == blk.interphase_mass_transfer[t, zf, j]
                )
            else:
                return blk.liquid_phase.mass_transfer_term[t, x, "Liq", j] == 0.0

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            self.vapor_phase.properties.component_list,
            doc="Vapor phase mass transfer constraint",
        )
        def vapor_mass_transfer_eqn(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            elif j in equilibrium_comp:
                return (
                    pyunits.convert(
                        blk.vapor_phase.mass_transfer_term[t, x, "Vap", j],
                        to_units=lunits("amount") / lunits("time") / lunits("length"),
                    )
                    == -blk.interphase_mass_transfer[t, x, j]
                )
            else:
                return blk.vapor_phase.mass_transfer_term[t, x, "Vap", j] == 0.0

        # ---------------------------------------------------------------------
        # Heat transfer relationships
        self.heat_transfer_coeff = Var(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            initialize=100,
            units=lunits("power") / lunits("temperature") / lunits("length"),
            doc="Vapor-liquid heat transfer coefficient",
        )

        # Heat transfer
        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Heat transfer calculation",
        )
        def heat_transfer_eqn1(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                return pyunits.convert(
                    blk.vapor_phase.heat[t, x],
                    to_units=lunits("power") / lunits("length"),
                ) == (
                    blk.heat_transfer_coeff[t, x]
                    * (
                        blk.liquid_phase.properties[t, zb].temperature
                        - pyunits.convert(
                            blk.vapor_phase.properties[t, x].temperature,
                            to_units=lunits("temperature"),
                        )
                    )
                )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Heat transfer balance",
        )
        def heat_transfer_eqn2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                zf = self.vapor_phase.length_domain.next(x)
                return blk.liquid_phase.heat[t, x] == -pyunits.convert(
                    blk.vapor_phase.heat[t, zf],
                    to_units=lunits("power") / lunits("length"),
                )

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            doc="Enthalpy transfer calculation",
        )
        def enthalpy_transfer_eqn1(blk, t, x):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                return blk.vapor_phase.enthalpy_transfer[t, x] == (
                    (
                        sum(
                            blk.vapor_phase.properties[t, x].enth_mol_phase_comp[
                                "Vap", j
                            ]
                            * blk.vapor_phase.mass_transfer_term[t, x, "Vap", j]
                            for j in equilibrium_comp
                        )
                    )
                )

        @self.Constraint(
            self.flowsheet().time,
            self.liquid_phase.length_domain,
            doc="Enthalpy transfer equality",
        )
        def enthalpy_transfer_eqn2(blk, t, x):
            if x == self.liquid_phase.length_domain.last():
                return Constraint.Skip
            else:
                zf = self.vapor_phase.length_domain.next(x)
                return blk.liquid_phase.enthalpy_transfer[t, x] == -pyunits.convert(
                    blk.vapor_phase.enthalpy_transfer[t, zf],
                    to_units=lunits("power") / lunits("length"),
                )

    def liquid_phase_mass_transfer_model(self):
        """
        Liquid phase mass transfer sub-model for calculating equilbrium
        partial pressure used in driving force.

        For the generic model, it is assumed that liquid phase mass transfer
        is sufficiently fast that the equilibrium partial pressure is equal to
        the fugacity of each component in the bulk.

        Derived Classes should overload this method as required.
        """
        vap_comp = self.config.vapor_phase.property_package.component_list
        liq_comp = self.config.liquid_phase.property_package.component_list
        equilibrium_comp = vap_comp & liq_comp

        @self.Constraint(
            self.flowsheet().time,
            self.vapor_phase.length_domain,
            equilibrium_comp,
            doc="Equilibrium partial pressure at interface",
        )
        def pressure_at_interface(blk, t, x, j):
            if x == self.vapor_phase.length_domain.first():
                return Constraint.Skip
            else:
                zb = self.liquid_phase.length_domain.prev(x)
                lprops = blk.liquid_phase.properties[t, zb]
                return blk.pressure_equil[t, x, j] == (lprops.fug_phase_comp["Liq", j])

    # =========================================================================
    # Scaling routine
    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()

        # ---------------------------------------------------------------------
        # Scale variables
        for (t, x, j), v in self.pressure_equil.items():
            if iscale.get_scaling_factor(v) is None:
                sf_pe = iscale.get_scaling_factor(
                    self.pressure_equil, default=None, warning=True
                )
                if sf_pe is None:
                    sf_pe = iscale.get_scaling_factor(
                        self.liquid_phase.properties[t, x].fug_phase_comp["Liq", j],
                        default=None,
                        warning=True,
                    )
                if sf_pe is None:
                    sf_pe = iscale.get_scaling_factor(
                        self.liquid_phase.properties[t, x].pressure,
                        default=1,
                        warning=True,
                    )

                iscale.set_scaling_factor(v, sf_pe)

        for (t, x), v in self.vapor_phase.heat.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.heat_transfer_coeff, default=1e-3, warning=True
                )
                iscale.set_scaling_factor(v, sf)

        for (t, x), v in self.liquid_phase.heat.items():
            if iscale.get_scaling_factor(v) is None:
                sf = iscale.get_scaling_factor(
                    self.heat_transfer_coeff, default=1e-3, warning=True
                )
                iscale.set_scaling_factor(v, sf)

        # ---------------------------------------------------------------------
        # Scale constraints
        for (t, x), v in self.mechanical_equilibrium.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.properties[t, x].pressure, default=1, warning=True
                ),
            )

        for (t, x, j), v in self.pressure_at_interface.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.pressure_equil[t, x, j], default=1, warning=False
                ),
            )

        for (t, x, j), v in self.interphase_mass_transfer_eqn.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.interphase_mass_transfer[t, x, j], default=1, warning=True
                ),
            )

        for (t, x, j), v in self.liquid_mass_transfer_eqn.items():
            try:
                sf = iscale.get_scaling_factor(
                    self.interphase_mass_transfer[t, x, j], default=1, warning=False
                )
            except KeyError:
                # This implies a non-volatile component
                sf = 1
            iscale.constraint_scaling_transform(v, sf)

        for (t, x, j), v in self.vapor_mass_transfer_eqn.items():
            try:
                sf = iscale.get_scaling_factor(
                    self.interphase_mass_transfer[t, x, j], default=1, warning=False
                )
            except KeyError:
                # This implies a non-volatile component
                sf = 1
            iscale.constraint_scaling_transform(v, sf)

        for (t, x), v in self.heat_transfer_eqn1.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.heat[t, x], default=1, warning=True
                ),
            )

        for (t, x), v in self.heat_transfer_eqn2.items():
            iscale.constraint_scaling_transform(
                v,
                iscale.get_scaling_factor(
                    self.vapor_phase.heat[t, x], default=1, warning=True
                ),
            )

    # =========================================================================
    # Model initialization routine
    def initialize(
        blk,
        vapor_phase_state_args=None,
        liquid_phase_state_args=None,
        state_vars_fixed=False,
        outlvl=idaeslog.NOTSET,
        solver=None,
        optarg=None,
    ):
        """
        Standard Packed Column initialization.

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

        unit_constraints = [
            "pressure_at_interface",
            "interphase_mass_transfer_eqn",
            "liquid_mass_transfer_eqn",
            "vapor_mass_transfer_eqn",
            "heat_transfer_eqn1",
            "heat_transfer_eqn2",
            "enthalpy_transfer_eqn1",
            "enthalpy_transfer_eqn2",
        ]

        # ---------------------------------------------------------------------
        # Deactivate unit model level constraints
        for c in blk.component_objects(Constraint, descend_into=True):
            if c.local_name in unit_constraints:
                c.deactivate()

        # Fix variables

        # Interface pressure
        blk.pressure_equil.fix()

        # Molar flux
        blk.interphase_mass_transfer.fix(0.0)
        blk.vapor_phase.mass_transfer_term.fix(0.0)
        blk.liquid_phase.mass_transfer_term.fix(0.0)

        # Heat transfer rate
        blk.vapor_phase.heat.fix(0.0)
        blk.liquid_phase.heat.fix(0.0)
        blk.vapor_phase.enthalpy_transfer.fix(0.0)
        blk.liquid_phase.enthalpy_transfer.fix(0.0)

        # ---------------------------------------------------------------------
        # Provide state arguments for property package initialization

        init_log.info("Step 1: Property Package initialization")
        # Initialize vapor_phase properties block
        vflag = blk.vapor_phase.initialize(
            state_args=vapor_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        # Initialize liquid_phase properties block
        lflag = blk.liquid_phase.initialize(
            state_args=liquid_phase_state_args,
            outlvl=outlvl,
            optarg=optarg,
            solver=solver,
            hold_state=True,
        )

        init_log.info("Step 2: Steady-State isothermal mass balance")

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 2: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 3: Interface equilibrium")

        # Activate interface pressure constraint
        blk.pressure_equil.unfix()
        blk.pressure_at_interface.activate()

        for k in blk.pressure_at_interface:
            calculate_variable_from_constraint(
                blk.pressure_equil[k], blk.pressure_at_interface[k]
            )

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 3 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 4a: Isothermal absoption")
        init_log.info_high("Calculating mass flux")

        # Unfix mass transfer terms
        blk.interphase_mass_transfer.unfix()

        # Activate mass transfer equation in vapor phase
        blk.interphase_mass_transfer_eqn.activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        init_log.info("Step 4b: Isothermal chemical absoption")
        init_log.info_high("Adding mass transfer to material balances")

        blk.vapor_phase.mass_transfer_term.unfix()
        blk.liquid_phase.mass_transfer_term.unfix()
        blk.vapor_mass_transfer_eqn.activate()
        blk.liquid_mass_transfer_eqn.activate()

        # Fix this
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)
        init_log.info_high("Step 4 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        init_log.info("Step 5: Adiabatic chemical absoption")
        init_log.info_high("Isothermal to Adiabatic ")

        # Unfix heat transfer terms
        blk.vapor_phase.heat.unfix()
        blk.liquid_phase.heat.unfix()
        blk.vapor_phase.enthalpy_transfer.unfix()
        blk.liquid_phase.enthalpy_transfer.unfix()

        # Activate heat transfer equations
        for c in [
            "heat_transfer_eqn1",
            "heat_transfer_eqn2",
            "enthalpy_transfer_eqn1",
            "enthalpy_transfer_eqn2",
        ]:
            getattr(blk, c).activate()

        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = opt.solve(blk, tee=slc.tee)

        init_log.info_high("Step 5 complete: {}.".format(idaeslog.condition(res)))

        # ---------------------------------------------------------------------
        blk.vapor_phase.release_state(flags=vflag)
        blk.liquid_phase.release_state(flags=lflag)

        if (
            res.solver.termination_condition != TerminationCondition.optimal
            or res.solver.status != SolverStatus.ok
        ):
            raise InitializationError(
                f"{blk.name} failed to initialize successfully. Please check "
                f"the output logs for more information."
            )
