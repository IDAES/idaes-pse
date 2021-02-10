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
from pyomo.environ import Reference

from idaes.generic_models.properties.core.generic.generic_property import StateIndex
from idaes.core.util.exceptions import BurntToast
from idaes.core.util.misc import add_object_reference
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


def define_electrolyte_state(b):
    if b.params.config.state_components == StateIndex.true:
        # Create references to base state vars
        add_object_reference(b, "flow_mol_true", b.flow_mol)
        b.flow_mol_phase_true = Reference(b.flow_mol_phase)
        b.flow_mol_phase_comp_true = Reference(b.flow_mol_phase_comp)
        b.mole_frac_phase_comp_true = Reference(b.mole_frac_phase_comp)
    elif b.params.config.state_components == StateIndex.apparent:
        # Create references to base state vars
        add_object_reference(b, "flow_mol_apparent", b.flow_mol)
        b.flow_mol_phase_apparent = Reference(b.flow_mol_phase)
        b.flow_mol_phase_comp_apparent = Reference(b.flow_mol_phase_comp)
        b.mole_frac_phase_comp_apparent = Reference(b.mole_frac_phase_comp)
    else:
        raise BurntToast("{} - unrecognized value for state_components "
                         "argument - this should never happen. Please "
                         "contact the IDAES developers".format(b.name))
