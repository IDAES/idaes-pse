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
Methods for creating additional state variables for electrolyte systems
"""
from pyomo.environ import Constraint, NonNegativeReals, Reference, Var

from idaes.generic_models.properties.core.generic.generic_property import \
    StateIndex
from idaes.generic_models.properties.core.generic.utility import \
    get_bounds_from_config
from idaes.core.util.exceptions import BurntToast
from idaes.core.util.misc import add_object_reference
from idaes.core.components import IonData
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


def define_electrolyte_state(b):
    if b.params.config.state_components == StateIndex.true:
        _true_species_state(b)
    elif b.params.config.state_components == StateIndex.apparent:
        _apparent_species_state(b)
    else:
        raise BurntToast("{} - unrecognized value for state_components "
                         "argument - this should never happen. Please "
                         "contact the IDAES developers".format(b.name))


def _apparent_species_state(b):
    # Create references to base state vars
    add_object_reference(b, "flow_mol_apparent", b.flow_mol)
    b.flow_mol_phase_apparent = Reference(b.flow_mol_phase)
    b.flow_mol_phase_comp_apparent = Reference(b.flow_mol_phase_comp)
    b.mole_frac_phase_comp_apparent = Reference(b.mole_frac_phase_comp)

    # Get units and bounds for true species state
    units = b.params.get_metadata().derived_units
    f_bounds, f_init = get_bounds_from_config(
        b, "flow_mol", units["flow_mole"])

    # Create true species state vars
    b.flow_mol_phase_comp_true = Var(
        b.params.true_phase_component_set,
        initialize=f_init,
        domain=NonNegativeReals,
        bounds=f_bounds,
        doc="Phase-component molar flowrates of true species",
        units=units["flow_mole"])

    # Check for inherent reactions and add apparent extent terms if required
    if b.has_inherent_reactions:
        b.inherent_reaction_extent = Var(
            b.params.inherent_reaction_idx,
            initialize=0,
            units=units["flow_mole"],
            doc="Apparent extent of inherent reactions")

    def appr_to_true_species(b, p, j):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_aqueous_phase:
            if isinstance(cobj, IonData):
                e = 0
                for a in b.params._apparent_set:
                    aobj = b.params.get_component(a)
                    if j in aobj.config.dissociation_species:
                        e += (aobj.config.dissociation_species[j] *
                              b.flow_mol_phase_comp_apparent[p, a])
            else:
                e = b.flow_mol_phase_comp_apparent[p, j]

            # Next, check for inherent reactions
            if b.has_inherent_reactions:
                for r in b.params.inherent_reaction_idx:
                    # Get stoichiometric coeffiicient for inherent reactions
                    gamma = b.params.inherent_reaction_stoichiometry[r, p, j]

                    if gamma != 0:
                        e += gamma*b.inherent_reaction_extent[r]

            return b.flow_mol_phase_comp_true[p, j] == e
        else:
            return b.flow_mol_phase_comp_apparent[p, j] == \
                b.flow_mol_phase_comp_true[p, j]

    b.appr_to_true_species = Constraint(
        b.params.true_phase_component_set,
        rule=appr_to_true_species,
        doc="Apparent to true species conversion")


def _true_species_state(b):
    # Create references to base state vars
    add_object_reference(b, "flow_mol_true", b.flow_mol)
    b.flow_mol_phase_true = Reference(b.flow_mol_phase)
    b.flow_mol_phase_comp_true = Reference(b.flow_mol_phase_comp)
    b.mole_frac_phase_comp_true = Reference(b.mole_frac_phase_comp)
