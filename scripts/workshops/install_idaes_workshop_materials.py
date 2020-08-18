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
This is a script that will be placed on pyomo.org so it can be downloaded and executed to update and install any workshop material required during the workshop.

TODO: This script should be moved to a new IDAES repository.
"""
import shutil

def execute():
    shutil.copytree('./Module_0_Welcome', '../../../jovyan/Module_0_Welcome')
    shutil.copytree('./Module_1_Flash_Unit', '../../../jovyan/Module_1_Flash_Unit')
    shutil.copytree('./Module_2_Flowsheet', '../../../jovyan/Module_2_Flowsheet')
    shutil.copytree('./Module_3_Custom_Unit_Model', '../../../jovyan/Module_3_Custom_Unit_Model')

