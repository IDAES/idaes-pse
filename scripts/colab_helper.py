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
# If running on Google Colab, install Pyomo and Ipopt via IDAES
if "google.colab" in sys.modules:
    !wget "https://raw.githubusercontent.com/IDAES/idaes-pse/main/scripts/colab_helper.py"
    import colab_helper
    colab_helper.install_idaes()
    colab_helper.install_ipopt()
```

For testing purposes, you may need to use:
"https://raw.githubusercontent.com/adowling2/idaes-pse/colab-install-script/scripts/colab_helper.py"

"""

__version__ = "2023.08.10"

import shutil
import sys
import os.path
import os
import re

import subprocess


def _check_available(executable_name):
    """Utility to check in an executable is available"""
    return shutil.which(executable_name) or os.path.isfile(executable_name)


def package_available(package_name):
    """Utility to check if a package/executable is available

    This supports customization, e.g., glpk, for special package names
    """

    if package_name == "glpk":
        return _check_available("gpsol")
    else:
        return _check_available(package_name)


def on_colab():
    """Utility returns True if executed on Colab, False otherwise"""
    return "google.colab" in sys.modules


def install_idaes(verbose=False):
    """Installs latest version of IDAES-PSE via pip

    Argument:
        verbose: bool, if True, display console output from pip install

    """

    try:
        import idaes

        print("idaes was found! No need to install.")
    except ImportError:
        print("Installing idaes via pip...")
        v = subprocess.run(
            [sys.executable, "-m", "pip", "install", "-q", "idaes_pse"],
            check=True,
            capture_output=True,
            text=True,
        )
        if verbose:
            print(v.stdout)
            print(v.stderr)
        print("idaes was successfully installed")
        v = subprocess.run(
            ["idaes", "--version"], check=True, capture_output=True, text=True
        )
        print(v.stdout)
        print(v.stderr)


def install_ipopt(verbose=False, try_conda_as_backup=False):
    """Install Ipopt and possibly other solvers.

    If running on Colab, this will install Ipopt, k_aug, and other COIN-OR
    solvers via idaes get-extensions.

    Arguments:
        verbose: bool, if True, display console output from idaes get-extensions and conda
        try_conda_as_backup: bool, if True, install ipopt via conda if idaes get-extensions fails
    """

    # Check if Ipopt (solver) is available. If not, install it.
    if not package_available("ipopt"):
        print("Running idaes get-extensions to install Ipopt, k_aug, and more...")
        v = subprocess.run(
            ["idaes", "get-extensions"], check=True, capture_output=True, text=True
        )
        if verbose:
            print(v.stdout)
            print(v.stderr)
        _update_path()
        print("Checking solver versions:")
        _print_solver_versions()

    # Check again if Ipopt is available. If not, try conda
    if try_conda_as_backup and not package_available("ipopt"):
        print("Installing Ipopt via conda...")
        v = subprocess.run(
            [sys.executable, "-m", "conda", "install", "-c", "conda-forge", "ipopt"],
            check=True,
            capture_output=True,
            text=True,
        )
        if verbose:
            print(v.stdout)
            print(v.stderr)
        print("Checking ipopt version:")
        _print_single_solver_version("ipopt")


def _update_path():
    """Add idaes executables to PATH"""
    if not re.search(re.escape("/root/.idaes/bin/"), os.environ["PATH"]):
        os.environ["PATH"] = "/root/.idaes/bin/:" + os.environ["PATH"]


def _print_single_solver_version(solvername):
    """Print the version for a single solver
    Arg:
        solvername: solver executable name (string)
    """
    v = subprocess.run([solvername, "-v"], check=True, capture_output=True, text=True)
    print(v.stdout)
    print(v.stderr)


def _print_solver_versions():
    """Print versions of solvers in idaes get-extensions

    This is the primary check that solvers installed correctly and are callable
    """

    # This does not work for cbc and clp; calling --version with these solvers,
    # enters their scripting language mode.
    for s in ["ipopt", "k_aug", "couenne", "bonmin", "ipopt_l1", "dot_sens"]:
        _print_single_solver_version(s)
