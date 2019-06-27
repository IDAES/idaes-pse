##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for update_workshop_materials
"""

import idaes.util.update_workshop_materials as up

def test_update_workshop_materials():
    # download the module from pyomo.org
    im = up.download_and_import_install_module()

    # test that the module loaded and can be called
    assert im.get_test_string() == 'Successfully called get_test_string()'

if __name__ == '__main__':
    test_update_workshop_materials()
    
