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
Methods for creating additional state variables for electrolyte systems
"""
from pyomo.environ import Constraint, NonNegativeReals, Reference, units as pyunits, Var

from idaes.models.properties.modular_properties.base.generic_property import StateIndex
from idaes.models.properties.modular_properties.base.utility import (
    get_bounds_from_config,
)
from idaes.core.util.exceptions import BurntToast
from idaes.core.util.misc import add_object_reference
from idaes.core.base.components import IonData, ApparentData
import idaes.logger as idaeslog
import idaes.core.util.scaling as iscale

# Set up logger
_log = idaeslog.getLogger(__name__)


def define_electrolyte_state(b):
    if b.params.config.state_components == StateIndex.true:
        _true_species_state(b)
    elif b.params.config.state_components == StateIndex.apparent:
        _apparent_species_state(b)
    else:
        raise BurntToast(
            "{} - unrecognized value for state_components "
            "argument - this should never happen. Please "
            "contact the IDAES developers".format(b.name)
        )


def calculate_electrolyte_scaling(b):
    if b.params.config.state_components == StateIndex.true:
        _true_species_scaling(b)
    elif b.params.config.state_components == StateIndex.apparent:
        _apparent_species_scaling(b)
    else:
        raise BurntToast(
            "{} - unrecognized value for state_components "
            "argument - this should never happen. Please "
            "contact the IDAES developers".format(b.name)
        )


def _apparent_species_state(b):
    # Create references to base state vars
    add_object_reference(b, "flow_mol_apparent", b.flow_mol)
    b.flow_mol_phase_apparent = Reference(b.flow_mol_phase)
    b.flow_mol_phase_comp_apparent = Reference(b.flow_mol_phase_comp)
    b.mole_frac_phase_comp_apparent = Reference(b.mole_frac_phase_comp)

    # Get units and bounds for true species state
    units = b.params.get_metadata().derived_units
    f_bounds, f_init = get_bounds_from_config(b, "flow_mol", units.FLOW_MOLE)

    # Create true species state vars
    b.flow_mol_phase_comp_true = Var(
        b.params.true_phase_component_set,
        initialize=f_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc="Phase-component molar flowrates of true species",
        units=units.FLOW_MOLE,
    )

    b.mole_frac_phase_comp_true = Var(
        b.params.true_phase_component_set,
        initialize=1 / len(b.params.true_species_set),
        bounds=(1e-20, 1.001),
        doc="Phase-component molar fractions of true species",
        units=pyunits.dimensionless,
    )

    # Check for inherent reactions and add apparent extent terms if required
    if b.has_inherent_reactions:
        b.apparent_inherent_reaction_extent = Var(
            b.params.inherent_reaction_idx,
            initialize=0,
            units=units.FLOW_MOLE,
            doc="Apparent extent of inherent reactions",
        )

    def appr_to_true_species(b, p, j):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_aqueous_phase:
            if isinstance(cobj, IonData):
                e = 0
                for a in b.params._apparent_set:
                    aobj = b.params.get_component(a)
                    if j in aobj.config.dissociation_species:
                        e += (
                            aobj.config.dissociation_species[j]
                            * b.flow_mol_phase_comp_apparent[p, a]
                        )
            else:
                e = b.flow_mol_phase_comp_apparent[p, j]

            # Next, check for inherent reactions
            if b.has_inherent_reactions:
                for r in b.params.inherent_reaction_idx:
                    # Get stoichiometric coeffiicient for inherent reactions
                    gamma = b.params.inherent_reaction_stoichiometry[r, p, j]

                    if gamma != 0:
                        e += gamma * b.apparent_inherent_reaction_extent[r]

            return b.flow_mol_phase_comp_true[p, j] == e
        else:
            return (
                b.flow_mol_phase_comp_apparent[p, j] == b.flow_mol_phase_comp_true[p, j]
            )

    b.appr_to_true_species = Constraint(
        b.params.true_phase_component_set,
        rule=appr_to_true_species,
        doc="Apparent to true species conversion",
    )

    def true_species_mole_fractions(b, p, j):
        return (
            b.mole_frac_phase_comp_true[p, j]
            * sum(
                b.flow_mol_phase_comp_true[p, k]
                for k in b.params.true_species_set
                if (p, k) in b.params.true_phase_component_set
            )
            == b.flow_mol_phase_comp_true[p, j]
        )

    b.true_mole_frac_constraint = Constraint(
        b.params.true_phase_component_set,
        rule=true_species_mole_fractions,
        doc="Calculation of true species mole fractions",
    )


def _apparent_species_scaling(b):
    for p, j in b.params.true_phase_component_set:
        sf = iscale.get_scaling_factor(
            b.flow_mol_phase_comp_true[p, j], default=1, warning=True
        )

        iscale.constraint_scaling_transform(
            b.appr_to_true_species[p, j], sf, overwrite=False
        )
        iscale.constraint_scaling_transform(
            b.true_mole_frac_constraint[p, j], sf, overwrite=False
        )


def _true_species_state(b):
    # Create references to base state vars
    add_object_reference(b, "flow_mol_true", b.flow_mol)
    b.flow_mol_phase_true = Reference(b.flow_mol_phase)
    b.flow_mol_phase_comp_true = Reference(b.flow_mol_phase_comp)
    b.mole_frac_phase_comp_true = Reference(b.mole_frac_phase_comp)

    # Get units and bounds for apparent species state
    units = b.params.get_metadata().derived_units
    f_bounds, f_init = get_bounds_from_config(b, "flow_mol", units.FLOW_MOLE)

    # Create apparent species state vars
    b.flow_mol_phase_comp_apparent = Var(
        b.params.apparent_phase_component_set,
        initialize=f_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc="Phase-component molar flowrates of apparent species",
        units=units.FLOW_MOLE,
    )

    b.mole_frac_phase_comp_apparent = Var(
        b.params.apparent_phase_component_set,
        initialize=1 / len(b.params.apparent_species_set),
        bounds=(1e-20, 1.001),
        doc="Phase-component molar fractions of apparent species",
        units=pyunits.dimensionless,
    )

    def true_to_appr_species(b, p, j):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_aqueous_phase():
            total_charge = sum(
                b.flow_mol_phase_comp_true[p, c]
                * b.params.get_component(c).config.charge
                for c in b.params.cation_set
            )

            if isinstance(cobj, ApparentData):
                # Need to recompose composition
                ions = list(cobj.config.dissociation_species.keys())

                e = (
                    b.flow_mol_phase_comp_true[p, ions[0]]
                    * abs(b.params.get_component(ions[0]).config.charge)
                    * b.flow_mol_phase_comp_true[p, ions[1]]
                    / (total_charge * cobj.config.dissociation_species[ions[1]])
                )
            elif j == "H2O":
                # Special case for water
                # Need to generalise to cover other cases with weak acids
                e = b.flow_mol_phase_comp_true[p, j]

                if "H+" in b.params.cation_set and "OH-" in b.params.anion_set:
                    e += (
                        b.flow_mol_phase_comp_true[p, "H+"]
                        * b.flow_mol_phase_comp_true[p, "OH-"]
                        / total_charge
                    )
                elif "H3O+" in b.params.cation_set and "OH-" in b.params.anion_set:
                    # Assume H3O+
                    e += (
                        2
                        * b.flow_mol_phase_comp_true[p, "H3O+"]
                        * b.flow_mol_phase_comp_true[p, "OH-"]
                        / total_charge
                    )
            else:
                e = b.flow_mol_phase_comp_true[p, j]

            return b.flow_mol_phase_comp_apparent[p, j] == e
        else:
            return (
                b.flow_mol_phase_comp_apparent[p, j] == b.flow_mol_phase_comp_true[p, j]
            )

    b.true_to_appr_species = Constraint(
        b.params.apparent_phase_component_set,
        rule=true_to_appr_species,
        doc="True to apparent species conversion",
    )

    def apparent_species_mole_fractions(b, p, j):
        return (
            b.mole_frac_phase_comp_apparent[p, j]
            * sum(
                b.flow_mol_phase_comp_apparent[p, k]
                for k in b.params.apparent_species_set
                if (p, k) in b.params.apparent_phase_component_set
            )
            == b.flow_mol_phase_comp_apparent[p, j]
        )

    b.appr_mole_frac_constraint = Constraint(
        b.params.apparent_phase_component_set,
        rule=apparent_species_mole_fractions,
        doc="Calculation of apparent species mole fractions",
    )


def _true_species_scaling(b):
    pass
