#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Standard IDAES Gibbs reactor model.
"""
# Import Pyomo libraries
from pyomo.environ import Constraint, Param, Reals, Reference, Set, units, value, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In, ListOf, Bool

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.constants import Constants

__author__ = "Jinliang Ma, Andrew Lee"


class GibbsReactorScaler(CustomScalerBase):
    """
    Scaler for Gibbs Reactor units.

    Due to the nature of Gibbs Reactors, scaling is highly dependent on the outlet
    concentrations which cannot be predicted a priori, thus we rely on users to
    provide the best initial guesses they can for the outlet concentrations.

    """

    UNIT_SCALING_FACTORS = {
        # "QuantityName: (reference units, scaling factor)
        "Delta Pressure": (units.Pa, 1e-3),
        "Heat": (units.J / units.s, 1e-6),
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Variable scaling routine for Gibbs reactors.

        Due to the nature of Gibbs Reactors, scaling is highly dependent on the outlet
        concentrations which cannot be predicted a priori, thus we rely on users to
        provide the best initial guesses they can for the outlet concentrations.

        Args:
            model: instance of GibbsReactor to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            None
        """
        if submodel_scalers is None:
            submodel_scalers = {}

        # Step 1: Property scaling

        # Step 1a: propagate any existing scaling from inlet to outlet
        # This is likely a very poor approximation for concentrations, but we expect users
        # to provide better scaling factors manually (as there is no way for us to know).
        # This will do a first pass fill in however if the user does not provide any information.
        self.propagate_state_scaling(
            target_state=model.control_volume.properties_out,
            source_state=model.control_volume.properties_in,
            overwrite=overwrite,
        )

        # Step 1b: Call Scalers for state blocks
        # Inlet properties
        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_in",
            submodel_scalers=submodel_scalers,
            method="variable_scaling_routine",
            overwrite=overwrite,
        )
        # Outlet properties
        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_out",
            submodel_scalers=submodel_scalers,
            method="variable_scaling_routine",
            overwrite=overwrite,
        )

        # Step 2: Scaling Gibbs reactor variables
        # Control volume variables - support only heat and deltaP
        if hasattr(model.control_volume, "heat"):
            for v in model.control_volume.heat.values():
                self.scale_variable_by_units(v, overwrite=overwrite)
        if hasattr(model.control_volume, "deltaP"):
            for v in model.control_volume.deltaP.values():
                self.scale_variable_by_units(v, overwrite=overwrite)

        # Lagrangian multipliers
        # Best guess scaling for these is R*T, need to convert units
        p_units = (
            model.control_volume.config.property_package.get_metadata().get_derived_units
        )
        for (t, _), v in model.lagrange_mult.items():
            tsf = self.get_scaling_factor(
                model.control_volume.properties_out[t].temperature
            )
            if tsf is not None:
                nominal_t = 1 / tsf
            else:
                nominal_t = 500
            lsf = value(
                1
                / units.convert(
                    Constants.gas_constant * nominal_t * p_units("temperature"),
                    to_units=p_units("energy_mole"),
                )
            )
            self.set_variable_scaling_factor(v, lsf, overwrite=overwrite)

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: dict = None
    ):
        """
        Routine to apply scaling factors to constraints in model.

        Constraints will be scaled based on nominal Jacobian norms, and thus will
        be heavily dependent on variable scaling.

        Args:
            model: model to be scaled
            overwrite: whether to overwrite existing scaling factors
            submodel_scalers: dict of Scalers to use for sub-models, keyed by submodel local name

        Returns:
            None
        """
        # Step 1: Call Scalers for state blocks
        # Inlet properties
        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_in",
            submodel_scalers=submodel_scalers,
            method="constraint_scaling_routine",
            overwrite=overwrite,
        )
        # Outlet properties
        self.call_submodel_scaler_method(
            model=model,
            submodel="control_volume.properties_out",
            submodel_scalers=submodel_scalers,
            method="constraint_scaling_routine",
            overwrite=overwrite,
        )

        # Step 2: Scale all the control volume constraints
        for cd in model.control_volume.component_data_objects(
            ctype=Constraint, descend_into=False
        ):
            self.scale_constraint_by_nominal_value(
                cd, scheme="inverse_sum", overwrite=overwrite
            )

        # Step 3: Scale local constraints
        # Scale Gibbs minimization constraints
        for cd in model.gibbs_minimization.values():
            self.scale_constraint_by_nominal_value(
                cd, scheme="inverse_sum", overwrite=overwrite
            )

        # Scale inert species balance if they are present
        if hasattr(model, "inert_species_balance"):
            for cd in model.inert_species_balance.values():
                self.scale_constraint_by_nominal_value(
                    cd, scheme="inverse_sum", overwrite=overwrite
                )


@declare_process_block_class("GibbsReactor")
class GibbsReactorData(UnitModelBlockData):
    """
    Standard Gibbs Reactor Unit Model Class

    This model assume all possible reactions reach equilibrium such that the
    system partial molar Gibbs free energy is minimized.
    Since some species mole flow rate might be very small,
    the natural log of the species molar flow rate is used.
    Instead of specifying the system Gibbs free energy as an objective
    function, the equations for zero partial derivatives of the grand function
    with Lagrangian multiple terms with respect to product species mole flow
    rates and the multiples are specified as constraints.
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Gibbs reactors do not support dynamic models, thus this must be
False.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag",
            doc="""Gibbs reactors do not have defined volume, thus this must be
False.""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
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
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
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
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )
    CONFIG.declare(
        "inert_species",
        ConfigValue(
            default=[],
            domain=ListOf(str),
            description="List of inert species",
            doc="List of species which do not take part in reactions.",
        ),
    )

    default_scaler = GibbsReactorScaler

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(GibbsReactorData, self).build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=False)

        # Validate list of inert species
        for i in self.config.inert_species:
            if i not in self.control_volume.properties_in.component_list:
                raise ConfigurationError(
                    "{} invalid component in inert_species argument. {} is "
                    "not in the property package component list.".format(self.name, i)
                )

        self.control_volume.add_total_element_balances()

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()

        # Add performance equations
        # Add Lagrangian multiplier variables
        # Only need multipliers for species which participate in reactions
        l_set = []
        for e in self.control_volume.config.property_package.element_list:
            skip = True
            for j in self.control_volume.properties_in.component_list:
                if j in self.config.inert_species:
                    continue
                elif self.control_volume.properties_out.params.element_comp[j][e] != 0:
                    skip = False

            if not skip:
                l_set.append(e)
        self.lagrange_set = Set(initialize=l_set)

        e_units = self.control_volume.config.property_package.get_metadata().get_derived_units(
            "energy_mole"
        )
        self.lagrange_mult = Var(
            self.flowsheet().time,
            self.lagrange_set,
            domain=Reals,
            initialize=100,
            doc="Lagrangian multipliers",
            units=e_units,
        )

        # TODO : Remove this once scaling is properly implemented
        self.gibbs_scaling = Param(default=1, mutable=True)

        # Use Lagrangian multiple method to derive equations for Out_Fi
        # Use RT*lagrange as the Lagrangian multiple such that lagrange is in
        # a similar order of magnitude as log(Yi)

        @self.Constraint(
            self.flowsheet().time,
            self.control_volume.properties_in.phase_component_set,
            doc="Gibbs energy minimisation constraint",
        )
        def gibbs_minimization(b, t, p, j):
            # Use natural log of species mole flow to avoid Pyomo solver
            # warnings of reaching infeasible point
            if j in self.config.inert_species:
                return Constraint.Skip
            return 0 == b.gibbs_scaling * (
                b.control_volume.properties_out[t].gibbs_mol_phase_comp[p, j]
                + sum(
                    b.lagrange_mult[t, e]
                    * b.control_volume.properties_out[t].params.element_comp[j][e]
                    for e in self.lagrange_set
                )
            )

        if len(self.config.inert_species) > 0:

            @self.Constraint(
                self.flowsheet().time,
                self.control_volume.properties_in.phase_list,
                self.config.inert_species,
                doc="Inert species balances",
            )
            def inert_species_balance(b, t, p, j):
                # Add species balances for inert components
                cv = b.control_volume
                e_comp = cv.properties_out[t].config.parameters.element_comp

                # Check for linear dependence with element balances
                # If an inert species is the only source of element e,
                # the inert species balance would be linearly dependent on the
                # element balance for e.
                dependent = True

                if len(self.control_volume.properties_in.phase_list) > 1:
                    # Multiple phases avoid linear dependency
                    dependent = False
                else:
                    for e in self.control_volume.config.property_package.element_list:
                        if e_comp[j][e] == 0:
                            # Element e not in component j, no effect
                            continue
                        else:
                            for i in self.control_volume.properties_in.component_list:
                                if i == j:
                                    continue
                                else:
                                    # If comp j shares element e with comp i
                                    # cannot be linearly dependent
                                    if e_comp[i][e] != 0:
                                        dependent = False

                if (
                    not dependent
                    and (p, j) in self.control_volume.properties_in.phase_component_set
                ):
                    return 0 == (
                        cv.properties_in[t].get_material_flow_terms(p, j)
                        - cv.properties_out[t].get_material_flow_terms(p, j)
                    )
                else:
                    return Constraint.Skip

        # Set references to balance terms at unit level
        if (
            self.config.has_heat_transfer is True
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            self.heat_duty = Reference(self.control_volume.heat[:])
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.control_volume.deltaP[:])

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}
