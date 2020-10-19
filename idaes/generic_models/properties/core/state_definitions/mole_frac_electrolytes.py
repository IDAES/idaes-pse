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

        for p in b.phase_list:
            p_obj = b.params.get_phase(p)
            if isinstance(p_obj, AqueousPhase):
                try:
                    basis = p_obj.config.equation_of_state_options[
                        "true_to_apparent_conversion"]
                except KeyError:
                    basis = "anion"

                # Add balances for all molecular species
                def rule_molecular(b, j):
                    return b.mole_frac_phase_comp_true[p, j] == \
                        b.mole_frac_phase_comp_apparent[p, j]
                b._molecular_mole_frac_equality = Constraint(
                    molec,
                    rule=rule_molecular,
                    doc="Equating mole fractions of molecular species")

                if basis == "anion":
                    print(1)
                elif basis == "cation":
                    print(2)
                else:
                    raise ConfigurationError(
                        "{} Unrecognized value for equation of state option "
                        "true_to_apparent_conversion: {}. Expected 'anion' or "
                        "'cation'".format(b.name, basis))
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
