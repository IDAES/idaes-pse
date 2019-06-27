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
This script downloads a python file from pyomo.org that will allow us to update 
the workshop material easily during a workshop.

The file install_idaes_workshop_materials.py is downloaded from pyomo.org, imported,
 and the method execute() is then called to do whatever actions are necessary.

To explain what happens:
- update_workshop_materials.py is a module in the idaes/util folder that does the work
- install_idaes_workshop_materials.py is posted on a site (for now pyomo.org, but later could be a repository on github).

1) JupyterHub: User executes update_workshop_materials.ipynb which calls to update_workshop_materials.py in idaes/util
2) update_workshop_materials.py downloads another python file (install_idaes_workshop_materials.py)
   from pyomo.org and calls "execute()" from that module.
3) install_idaes_workshop_materials.py does whatever is necessary to get the workshop materials 
   into the user folder (and perform any updates necessary).

TODO: This should probably be changed to get a zip file from an IDAES repository rather than install_idaes_workshop_materials.py from pyomo.org
"""

import os
import logging
import pyomo.common.fileutils as futils
import pyomo.common.download as dload
from pyutilib.misc import import_file

_install_idaes_workshop_materials_url = 'http://www.pyomo.org/s/install_idaes_workshop_materials.py'

def download_and_import_install_module():
    """
    Downloads install_idaes_workshop_materials.py from pyomo.org and imports the module.
    """
    download_dir = futils.this_file_dir()
    download_dir = os.path.join(download_dir, '../../examples/workshops')
    download_dest = os.path.join(download_dir, 'install_idaes_workshop_materials.py')
    
    print("\n\n"
          "#######################################################################################\n"
          "# Downloading: {}\n"
          "# to: {}\n"
          "#######################################################################################\n"
          "".format(_install_idaes_workshop_materials_url, download_dest)
    )

    if not os.path.isdir(download_dir):
        raise NameError('Unable to locate download directory: {}'.format(download_dir))

    try:
        downloader = dload.FileDownloader()
        downloader.set_destination_filename(download_dest)
        downloader.get_binary_file(_install_idaes_workshop_materials_url)
    except:
        print("\n\n***\nFailed to download: {}\n***".format(_install_idaes_workshop_materials_url))
        raise

    print('... download complete')
    print('... installing workshop materials')
    install_module = import_file(download_dest)
    return install_module

def download_and_execute_install_module():
    """
    Downloads, imports, and calls execute on the install_idaes_workshop_materials.py file from pyomo.org
    """
    install_module = download_and_import_install_module()
    install_module.execute()
    
if __name__ == '__main__':
    download_and_execute_install_module()

