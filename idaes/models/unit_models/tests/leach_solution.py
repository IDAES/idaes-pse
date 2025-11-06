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
A simple property package of sulfuric acid solution used for testing unit models'
compatibility with inherent reactions.

We do not want to make this property package an example because it uses hours for 
units of time instead of seconds. As a result, the derived units for pressure are
kg/(m * hr**2) instead of kg/(m * s**2), i.e. Pascals, which causes trouble when
pressure drop terms are included.
"""
from pyomo.environ import (
    ComponentMap,
    Constraint,
    Param,
    Set,
    Var,
    units as pyunits,
)
from idaes.core import (
    declare_process_block_class,
    StateBlockData,
    StateBlock,
    PhysicalParameterBlock,
    Phase,
    Component,
    MaterialFlowBasis,
)
from idaes.core.scaling import CustomScalerBase
from idaes.core.util.initialization import fix_state_vars

__author__ = "Andrew Lee, Douglas Allan"


class LeachSolutionScaler(CustomScalerBase):
    """
    Scaler object for leach solution properties
    """

    DEFAULT_SCALING_FACTORS = {
        "flow_vol": 0.1,
        "conc_mass_comp[H2O]": 1e-6,
        "conc_mass_comp[H]": 1e-1,
        "conc_mass_comp[HSO4]": 1e-2,
        "conc_mass_comp[SO4]": 1e-1,
    }

    def variable_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        self.scale_variable_by_default(model.flow_vol, overwrite=overwrite)
        for var in model.conc_mass_comp.values():
            self.scale_variable_by_default(var, overwrite=overwrite)

        for idx, var in model.conc_mol_comp.items():
            self.scale_variable_by_definition_constraint(
                var, model.molar_concentration_constraint[idx], overwrite=overwrite
            )

    def constraint_scaling_routine(
        self, model, overwrite: bool = False, submodel_scalers: ComponentMap = None
    ):
        for idx, con in model.molar_concentration_constraint.items():
            # Constraint is on a mass basis
            self.scale_constraint_by_component(
                con, model.conc_mass_comp[idx], overwrite=overwrite
            )
        if not model.config.defined_state:
            self.scale_constraint_by_component(
                model.h2o_concentration,
                model.conc_mass_comp["H2O"],
                overwrite=overwrite,
            )
            self.scale_constraint_by_nominal_value(
                model.hso4_dissociation, overwrite=overwrite
            )


@declare_process_block_class("LeachSolutionParameters")
class LeachSolutionParameterData(PhysicalParameterBlock):
    """
    Parameter block for simple property package that has an inherent reaction
    """

    def build(self):
        super().build()

        self.liquid = Phase()

        # Solvent
        self.H2O = Component()

        # Acid related species
        self.H = Component()
        self.HSO4 = Component()
        self.SO4 = Component()

        self.mw = Param(
            self.component_list,
            units=pyunits.kg / pyunits.mol,
            initialize={
                "H2O": 18e-3,
                "H": 1e-3,
                "HSO4": 97e-3,
                "SO4": 96e-3,
            },
        )

        # Inherent reaction for partial dissociation of HSO4
        self._has_inherent_reactions = True
        self.inherent_reaction_idx = Set(initialize=["Ka2"])
        self.inherent_reaction_stoichiometry = {
            ("Ka2", "liquid", "H"): 1,
            ("Ka2", "liquid", "HSO4"): -1,
            ("Ka2", "liquid", "SO4"): 1,
            ("Ka2", "liquid", "H2O"): 0,
        }
        self.Ka2 = Param(
            initialize=10**-1.99,
            mutable=True,
            units=pyunits.mol / pyunits.L,
        )

        # Assume dilute acid, density of pure water
        self.dens_mol = Param(
            initialize=1,
            units=pyunits.kg / pyunits.litre,
            mutable=True,
        )

        self._state_block_class = LeachSolutionStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "flow_vol": {"method": None},
                "conc_mass_comp": {"method": None},
                "conc_mol_comp": {"method": None},
                "dens_mol": {"method": "_dens_mol"},
            }
        )
        obj.add_default_units(
            {
                "time": pyunits.hour,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


class _LeachSolutionStateBlock(StateBlock):
    default_scaler = LeachSolutionScaler

    def fix_initialization_states(self):
        """
        Fixes state variables for state blocks.

        Returns:
            None
        """
        # Fix state variables
        fix_state_vars(self)

        # Deactivate inherent reactions
        for k in self:
            if not self[k].config.defined_state:
                self[k].h2o_concentration.deactivate()
                self[k].hso4_dissociation.deactivate()


@declare_process_block_class(
    "LeachSolutionStateBlock", block_class=_LeachSolutionStateBlock
)
class LeachSolutionStateBlockData(StateBlockData):
    """
    State block for simple property package that has an inherent reaction
    """

    def build(self):
        super().build()

        self.flow_vol = Var(
            units=pyunits.L / pyunits.hour,
            bounds=(1e-8, None),
        )
        self.conc_mass_comp = Var(
            self.params.component_list,
            units=pyunits.mg / pyunits.L,
            bounds=(1e-10, None),
        )
        self.conc_mol_comp = Var(
            self.params.component_list,
            units=pyunits.mol / pyunits.L,
            initialize=1e-5,
            bounds=(1e-8, None),
        )

        # Concentration conversion constraint
        @self.Constraint(self.params.component_list)
        def molar_concentration_constraint(b, j):
            return (
                pyunits.convert(
                    b.conc_mol_comp[j] * b.params.mw[j],
                    to_units=pyunits.mg / pyunits.litre,
                )
                == b.conc_mass_comp[j]
            )

        if not self.config.defined_state:
            # Concentration of H2O based on assumed density
            self.h2o_concentration = Constraint(
                expr=self.conc_mass_comp["H2O"] == 1e6 * pyunits.mg / pyunits.L
            )
            # Equilibrium for partial dissociation of HSO4
            self.hso4_dissociation = Constraint(
                expr=self.conc_mol_comp["HSO4"] * self.params.Ka2
                == self.conc_mol_comp["H"] * self.conc_mol_comp["SO4"]
            )

    def get_material_flow_terms(self, p, j):
        # Note conversion to mol/hour
        if j == "H2O":
            # Assume constant density of 1 kg/L
            return self.flow_vol * self.params.dens_mol / self.params.mw[j]
        else:
            # Need to convert from moles to mass
            return pyunits.convert(
                self.flow_vol * self.conc_mass_comp[j] / self.params.mw[j],
                to_units=pyunits.mol / pyunits.hour,
            )

    def get_material_flow_basis(self):
        return MaterialFlowBasis.molar

    def define_state_vars(self):
        return {
            "flow_vol": self.flow_vol,
            "conc_mass_comp": self.conc_mass_comp,
        }
