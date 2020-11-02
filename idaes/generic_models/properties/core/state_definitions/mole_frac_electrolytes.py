##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Methods for setting up molefractions for electrolyte systems
"""
from pyomo.environ import Constraint, Reference, Var

from idaes.core.phases import AqueousPhase
from idaes.generic_models.properties.core.generic.generic_property import \
    StateIndex
from idaes.core.util.exceptions import BurntToast, ConfigurationError


def create_mole_frac_vars(b):
    # First, create Vars for true and apparent mole fractions
    b.mole_frac_phase_comp_true = Var(
        b.params.true_phase_component_set,
        initialize=1/len(b.params.true_species_set),
        bounds=(0, None),
        doc='True phase mole fractions',
        units=None)

    b.mole_frac_phase_comp_apparent = Var(
        b.params.apparent_phase_component_set,
        initialize=1/len(b.params.apparent_species_set),
        bounds=(0, None),
        doc='Apparent phase mole fractions',
        units=None)

    # Make a Reference for the master mole fraction
    if b.params.config["state_components"] == StateIndex.true:
        b.mole_frac_phase_comp = Reference(b.mole_frac_phase_comp_true)
    elif b.params.config["state_components"] == StateIndex.apparent:
        b.mole_frac_phase_comp = Reference(b.mole_frac_phase_comp_apparent)
    else:
        raise BurntToast(
            "{} unrecognized value for configuration argument "
            "'state_components'; this should never happen. Please contact "
            "the IDAES developers with this bug.".format(b.name))

    # Finally, create constraints that link the two mole fraction vars
    # Determine which set of constraints to construct
    if b.config.defined_state:
        # For inlets and similar cases, constraints are based on state vars
        cons_set = b.params.config["state_components"]
    else:
        # For other state blocks, constraints based on material balance form
        # TODO: Implement config argument for this
        # For now, just use state vars
        cons_set = b.params.config["state_components"]

    # Constraints for non-aqueous phases - mole fractions are equal
    def rule_non_aqueous(b, p, j):
        if isinstance(b.params.get_phase(p), AqueousPhase):
            # Aqueous phases will be covered in other constraints
            return Constraint.Skip
        else:
            # For non-aqueous phases, equate mole fractions
            return b.mole_frac_phase_comp_true[p, j] == \
                b.mole_frac_phase_comp_apparent[p, j]
    b._non_aqueous_mole_frac_equality = Constraint(
        b.phase_component_set,
        rule=rule_non_aqueous,
        doc="Equating true and apparent mole fractions in non-aqueous phases")

    if cons_set == StateIndex.true:
        # True to apparent conversion
        anion = b.params.anion_set
        cation = b.params.cation_set
        molec = b.params.solvent_set|b.params.solute_set

        # TODO: Need different extents for each aqueous phase
        b._extent_apparent = Var(b.params._apparent_set,
                                 initialize=0,
                                 doc="Extent of apparent reactions")

        if "H3O+" in cation and "OH-" in anion:
            # Water self-ionization uses hydronium form
            gamma_H2O = 2
        elif "H+" in cation and "OH-" in anion:
            # Water self-ionization uses proton form
            gamma_H2O = 1
        else:
            # No water self-ionization
            gamma_H2O = 0

        if gamma_H2O > 0:
            b._extent_apparent_H2O = Var(initialize=0,
                                         doc="Extent of water self-ionization")

        def true_to_apparent_conversion(b, p, j):
            p_obj = b.params.get_phase(p)
            c_obj = b.params.get_component(j)
            if not isinstance(p_obj, AqueousPhase):
                # Non-aqueous phases covered earlier
                return Constraint.Skip
            elif not c_obj._is_phase_valid(p_obj):
                # Component not valid in aqueous phase
                return Constraint.Skip
            elif j == "H2O" and gamma_H2O > 0:
                return (b.mole_frac_phase_comp_apparent[p, j] ==
                        b.mole_frac_phase_comp_true[p, j] +
                        gamma_H2O*b._extent_apparent_H2O)
            elif j in molec:
                return (b.mole_frac_phase_comp_apparent[p, j] ==
                        b.mole_frac_phase_comp_true[p, j])
            else:
                # Can only write one constraint for between H+/H3O+ and OH-
                # Select based on acidic or basic config argument
                if (j in ["H+", "H3O+"] and
                        p_obj.config.equation_of_state_options["pH_range"] ==
                        "basic"):
                    return Constraint.Skip
                if (j == "OH-" and
                        p_obj.config.equation_of_state_options["pH_range"] ==
                        "acidic"):
                    return Constraint.Skip

                if j not in anion and j not in cation:
                    lhs = b.mole_frac_phase_comp_apparent[p, j]
                    rhs = 0
                else:
                    lhs = 0
                    rhs = b.mole_frac_phase_comp_true[p, j]

                for i in b.params._apparent_set:
                    i_comp = b.params.get_component(i)
                    if j in i_comp.config.dissociation_species.keys():
                        gamma = i_comp.config.dissociation_species[j]
                        rhs = rhs + gamma*b._extent_apparent[i]

                return lhs == rhs

        b._aqueous_mole_frac_equality = Constraint(
            b.phase_list,
            b.params.component_list,  # in this case, we want ALL components
            rule=true_to_apparent_conversion,
            doc="Relating true and apparent mole fractions in aqueous phases")

    elif cons_set == StateIndex.apparent:
        # Apparent to true conversion
        all_aqueous = (b.params.anion_set |
                       b.params.cation_set |
                       b.params.solvent_set |
                       b.params.solute_set)
    else:
        raise BurntToast(
            "{} Something went wrong trying to create mole fraction conversion"
            " constraints; this should never happen. Please contact "
            "the IDAES developers with this bug.".format(b.name))
