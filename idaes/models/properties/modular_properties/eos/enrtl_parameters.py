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
Sub-methods for eNRTL activity coefficient method.

Includes temperature dependance rules for alpha and tau
"""
from pyomo.environ import Reals, units as pyunits, Var

from idaes.core.util.exceptions import BurntToast, ConfigurationError
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


class ConstantAlpha(object):
    @staticmethod
    def build_parameters(b):
        param_block = b.parent_block()

        # Get user provided values for alpha (if present)
        try:
            alpha_data = param_block.config.parameter_data[b.local_name + "_alpha"]
        except KeyError:
            alpha_data = {}

        # Check for unused parameters in alpha_data
        for (i, j) in alpha_data.keys():
            if (i, j) not in b.component_pair_set_symmetric and (
                j,
                i,
            ) not in b.component_pair_set_symmetric:
                raise ConfigurationError(
                    "{} eNRTL alpha parameter provided for invalid "
                    "component pair {}. Please check typing and only provide "
                    "parameters for valid species pairs.".format(b.name, (i, j))
                )

        def alpha_init(b, i, j):
            if (i, j) in alpha_data.keys():
                v = alpha_data[(i, j)]
                # Check for non-symmetric value assignment
                if (j, i) in alpha_data.keys():
                    if alpha_data[(j, i)] != v:
                        raise ConfigurationError(
                            "{} eNRTL alpha parameter assigned non-symmetric "
                            "value for pair {}. Please assign only one value "
                            "for component pair.".format(b.name, (i, j))
                        )
                    else:
                        _log.info(
                            "eNRTL alpha value provided for both {} and "
                            "{}. It is only necessary to provide a "
                            "value for one of these due to symmetry.".format(
                                (i, j), (j, i)
                            )
                        )
            elif (j, i) in alpha_data.keys():
                v = alpha_data[(j, i)]
            elif (i in param_block.solvent_set or i in param_block.solute_set) and (
                j in param_block.solvent_set or j in param_block.solute_set
            ):
                # Molecular-molecular interaction, default value is 0.3
                v = 0.3
            else:
                # All other intereactions have default value 0.2
                v = 0.2
            return v

        b.add_component(
            "alpha",
            Var(
                b.component_pair_set_symmetric,
                within=Reals,
                initialize=alpha_init,
                doc="Symmetric non-randomness parameters",
                units=pyunits.dimensionless,
            ),
        )

    @staticmethod
    def return_expression(b, pobj, i, j, T):
        if (i, j) in pobj.alpha:
            return pobj.alpha[i, j]
        elif (j, i) in pobj.alpha:
            return pobj.alpha[j, i]
        elif i == j:
            return 0.2
        else:
            raise BurntToast(
                "{} alpha rule encountered unexpected index {}. Please contact"
                "the IDAES Developers with this bug.".format(b.name, (i, j))
            )


class ConstantTau(object):
    @staticmethod
    def build_parameters(b):
        param_block = b.parent_block()

        # Get user provided values for tau (if present)
        try:
            tau_data = param_block.config.parameter_data[b.local_name + "_tau"]
        except KeyError:
            tau_data = {}

        # Check for unused parameters in tau_data
        for (i, j) in tau_data.keys():
            if (i, j) not in b.component_pair_set:
                raise ConfigurationError(
                    "{} eNRTL tau parameter provided for invalid "
                    "component pair {}. Please check typing and only provide "
                    "parameters for valid species pairs.".format(b.name, (i, j))
                )

        def tau_init(b, i, j):
            if (i, j) in tau_data.keys():
                v = tau_data[(i, j)]
            else:
                # Default interaction value is 0
                v = 0
            return v

        b.add_component(
            "tau",
            Var(
                b.component_pair_set,
                within=Reals,
                initialize=tau_init,
                doc="Binary interaction energy parameters",
                units=pyunits.dimensionless,
            ),
        )

    @staticmethod
    def return_expression(b, pobj, i, j, T):
        if (i, j) in pobj.tau:
            return pobj.tau[i, j]
        elif (j, i) in pobj.tau:
            return pobj.tau[j, i]
        elif i == j:
            return 0
        else:
            raise BurntToast(
                "{} tau rule encountered unexpected index {}. Please contact"
                "the IDAES Developers with this bug.".format(b.name, (i, j))
            )
