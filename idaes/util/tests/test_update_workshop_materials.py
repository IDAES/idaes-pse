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
    # Note: This tests that the install_idaes_workshop_materials succeeds
    # We purposely do NOT test code within the file (don't want
    # to run code outside our repository on testing)

    # download the module from pyomo.org
    download_dest = up.download_install_module()

    download_succeeded = False
    with open(download_dest, 'r') as fd:
        for line in fd:
            if line.startswith('def execute():'):
                download_succeeded = True
                break

    assert download_succeeded == True

if __name__ == '__main__':
    test_update_workshop_materials()
    
