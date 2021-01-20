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
import types

from .eos_base import EoSBase
from idaes.core.util.exceptions import BurntToast
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class ENRTL(EoSBase):
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

    def build_parameters(b):
        pass
