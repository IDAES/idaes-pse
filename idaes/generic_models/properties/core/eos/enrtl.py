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
Methods for eNRTL activity coeffient method.

Only applicable to liquid/electrolyte phases
"""
from .eos_base import EoSBase
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class ENRTL(EoSBase):
    def common(b, pobj):
        pass

    def build_parameters(b):
        pass
