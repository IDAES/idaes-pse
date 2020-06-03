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
from idaes.core.property_base import PhysicalParameterBlock


class IndexMePlease1(PhysicalParameterBlock):

    @classmethod
    def define_metadata(cls, m):
        m.add_default_units({'temperature': 'K'})
        m.add_properties({'pressure': {'units': 'Pa', 'method': 'foo'},
                          'temperature': {'method': 'bar'}})
