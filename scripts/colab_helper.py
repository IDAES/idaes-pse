#!/usr/bin/env python
###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
###############################################################################
"""
Install IDAES, Ipopt, and other solvers on Google Colab

Created by Alex Dowling (adowling@nd.edu) and Jeff Kantor at the University of Notre Dame.

To use this script, add the following to a code block in a Jupyter notebook:

```
# Ipopt installer
import sys
# If running on Google Colab, install Ipopt via IDAES
if "google.colab" in sys.modules:
    !wget "https://raw.githubusercontent.com/idaes-pse/main/scripts/colab_helper.py"
    import colab_helper
    colab_helper.install_idaes()
    colab_helper.install_ipopt()

# Otherwise, attempt to load idaes which should inlcude Ipopt
# See https://idaes-pse.readthedocs.io/en/stable/tutorials/getting_started/index.html
# for instructions on running idaes get-extensions
else:
    try:
        import idaes
        # Provided idaes extensions are installed, importing idaes provides access to
        print("Successfully loaded idaes.")
    except:
        print("IDAES is not installed!!!")
```

"""

import shutil
import sys
import os.path
import os
import urllib

import subprocess

def _check_available(executable_name): return (shutil.which(executable_name) or os.path.isfile(executable_name)) 

def package_available(package_name):
    
    if package_name == "glpk":
        return _check_available("gpsol")        
    else:
        return _check_available(package_name)

def on_colab(): return "google.colab" in sys.modules

def command_with_output(command):
    r = subprocess.getoutput(command)
    print(r)


def install_idaes():
    ''' Installs latest version of IDAES-PSE via pip
    
    '''

    # Check if idaes is available. If not, install it
    if not package_available("idaes"):
        # Tip: "!pip" means run the 'pip' command outside the notebook.
        print("Installing idaes via pip...")
        os.system("pip install -q idaes_pse")
        assert package_available("idaes"), "idaes was not successfully installed."
        print("idaes was successfully installed")
        os.system('idaes --version')
    else:
        print("IDAES found! No need to install.")

def install_ipopt():
    ''' Install Ipopt and possibly other solvers. 
    
    If running on Colab, this will install Ipopt, k_aug, and other COIN-OR
    solvers via idaes get-extensions.
    '''

    # Check if Ipopt (solver) is available. If not, install it.
    if not package_available("ipopt"):
        # Check if we are running on Google Colab
        if on_colab():
            # Install idaes solvers
            print("Running idaes get-extensions to install Ipopt, k_aug, and more...")
            os.system("idaes get-extensions")

            # Add symbolic link for idaes solvers
            # These are needed for Colab to find the solvers
            os.system("ln -s /root/.idaes/bin/ipopt ipopt")
            os.system("ln -s /root/.idaes/bin/k_aug k_aug")
            os.system("ln -s /root/.idaes/bin/couenne couenne")
            os.system("ln -s /root/.idaes/bin/bonmin bonmin")
            os.system("ln -s /root/.idaes/bin/cbc cbc")
            os.system("ln -s /root/.idaes/bin/clp clp")
            os.system("ln -s /root/.idaes/bin/ipopt_l1 ipopt_l1")
            os.system("ln -s /root/.idaes/bin/dot_sens dot_sens")
            
            command_with_output('./ipopt -v')
            command_with_output('./k_aug -v')
            command_with_output('./couenne -v')
            command_with_output('./bonmin -v')
            #command_with_output('./cbc -v')
            #command_with_output('./clp -v')
            command_with_output('./ipopt_l1 -v')
            

        # Check again if Ipopt is available
        if not package_available("ipopt"):

            # As a backup, try copying the executable from AMPL
            if on_colab():
                print("Installing Ipopt via zip file...")
                os.system('wget -N -q "https://ampl.com/dl/open/ipopt/ipopt-linux64.zip"')
                os.system('!unzip -o -q ipopt-linux64')
            # Otherwise, try conda
            else:
                try:
                    print("Installing Ipopt via conda...")
                    os.system('conda install -c conda-forge ipopt')
                except:
                    pass
        

    else:
        print("Ipopt found! No need to install.")
        

    # Verify Ipopt is now available
    assert package_available("ipopt"), "Ipopt is not available"
    
    print("ipopt was successfully installed")
    
    solvers = ["k_aug", "cbc", "clp", "bonmin", "couenne", "ipopt_l1"]
    for s in solvers:
        if package_available(s):
            print(s,"was successfuly installed")
    print(" ")