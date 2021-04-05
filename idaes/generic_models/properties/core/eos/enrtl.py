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
Methods for eNRTL activity coefficient method.

Only applicable to liquid/electrolyte phases
"""

from pyomo.environ import Reals, units as pyunits, Var

from .eos_base import EoSBase
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class ENRTL(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    def common(b, pobj):
        pass
        # def _return_component_list(b):
        #     if b.config.species_basis == "true":
        #         return b.config.parameters.true_species_set
        #     elif b.config.species_basis == "apparent":
        #         return b.config.parameters.apparent_species_set
        #     else:
        #         raise BurntToast(
        #             "{} unrecognised value for configuration argument "
        #             "'species_basis'. This should never happen. Please "
        #             "contact the IDAES developers with this bug"
        #             .format(b.name))
        # b._return_component_list = types.MethodType(_return_component_list, b)

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    def build_parameters(b):
        param_block = b.parent_block()

        # Build symmetric indexing set
        sym_set = []
        for i in param_block.apparent_species_set:
            for j in param_block.apparent_species_set:
                if (j, i) not in sym_set:
                    sym_set.append((i, j))

        # Get user provided values for alpha (if present)
        try:
            alpha_data = param_block.config.parameter_data[
                b.local_name+"_alpha"]
        except KeyError:
            alpha_data = {}

        # Check for unused parameters in alpha_data
        for (i, j) in alpha_data.keys():
            if (i, j) not in sym_set and (j, i) not in sym_set:
                raise ConfigurationError(
                    "{} eNRTL alpha parameter provided for invalid "
                    "component pair {}. Please check typing and only provide "
                    "parameters for apparent species pairs."
                    .format(b.name, (i, j)))

        def alpha_init(b, i, j):
            if (i, j) in alpha_data.keys():
                v = alpha_data[(i, j)]
                # Check for non-symmetric value assignment
                if (j, i) in alpha_data.keys():
                    if alpha_data[(j, i)] != v:
                        raise ConfigurationError(
                            "{} eNRTL alpha parameter assigned non-symmetric "
                            "value for pair {}. Please assign only one value "
                            "for component pair.".format(b.name, (i, j)))
                    else:
                        _log.info("eNRTL alpha value provided for both {} and "
                                  "{}. It is only necessary to provide a "
                                  "value for one of these due to symmetry."
                                  .format((i, j), (j, i)))
            elif(j, i) in alpha_data.keys():
                v = alpha_data[(j, i)]
            elif ((i in param_block.solvent_set or
                   i in param_block.solute_set) and
                  (j in param_block.solvent_set or
                   j in param_block.solute_set)):
                # Molecular-molecular interaction, default value is 0.3
                v = 0.3
            else:
                # All other intereactions have default value 0.2
                v = 0.2
            return v

        b.add_component(
            'alpha',
            Var(sym_set,
                within=Reals,
                initialize=alpha_init,
                doc='Symmetric non-randomness parameters',
                units=pyunits.dimensionless))

    @staticmethod
    def dens_mol_phase(b, p):
        return 55e3

    @staticmethod
    def enth_mol_phase(b, p):
        return 1e2*b.temperature

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        return 1e2*b.temperature
