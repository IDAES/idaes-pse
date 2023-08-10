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

Created by Alex Dowling (adowling@nd.edu) and Jeff Kantor at the University of Notre Dame
with input from John Siirola at Sandia National Laboratories.

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

def _print_solver_version(solvername):
    subprocess.run([sys.executable, "-m", solvername, "-v"], check=True)


def install_idaes():
    ''' Installs latest version of IDAES-PSE via pip
    
    '''

    try:
        import idaes
        print("idaes was found! No need to install.")
    except ImportError:
        print("Installing idaes via pip...")
        subprocess.run([sys.executable, "-m", "pip", "install", "-q", "idaes_pse"], check=True)
        print("idaes was successfully installed")
        subprocess.run(["idaes", "--version"], check=True)

def install_ipopt(try_conda_as_backup=False):
    ''' Install Ipopt and possibly other solvers. 
    
    If running on Colab, this will install Ipopt, k_aug, and other COIN-OR
    solvers via idaes get-extensions.
    '''

    # Check if Ipopt (solver) is available. If not, install it.
    if not package_available("ipopt"):
        print("Running idaes get-extensions to install Ipopt, k_aug, and more...")
        subprocess.run(["idaes", "get-extensions"], check=True)

    # Check again if Ipopt is available. If not, try conda
    if try_conda_as_backup and not package_available("ipopt"):
        print("Installing Ipopt via conda...")
        subprocess.run([sys.executable, "-m", "conda", "install", "-c", 
                        "conda-forge","ipopt"], check=True)
    else:
        print("Ipopt found! No need to install.")
        

    """
    # Verify Ipopt is now available
    assert package_available("ipopt"), "Ipopt is not available"
    
    print("ipopt was successfully installed")
    
    solvers = ["k_aug", "cbc", "clp", "bonmin", "couenne", "ipopt_l1"]
    for s in solvers:
        if package_available(s):
            print(s,"was successfuly installed")
    print(" ")
    """

def create_solver_symbolic_links():
    """ Create symbolic links to use solvers on Colab
    """

    for s in ["ipopt", "k_aug", "couenne", "bonmin", "cbc", 
              "clp", "ipopt_l1", "dot_sens"]:
        subprocess.run([sys.executable, "-m", "ln", "/root/.idaes/bin/"+s, s], check=True)

def print_solver_versions():

    # This does not work for cbc and clp. Not sure why
    for s in ["ipopt", "k_aug", "couenne", "bonmin", 
              "ipopt_l1", "dot_sens"]:
        subprocess.run([sys.executable, "-m", s, "-v"], check=True)