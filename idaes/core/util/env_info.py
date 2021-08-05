#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################

__Author__ = "John Eslick"

import sys
import platform
from collections import OrderedDict
import pkg_resources

import pyomo

import idaes.ver as ver

class EnvironmentInfo():
    """Get information about IDAES and the environment IDAES is running in,
    including OS and Critial dependency versions."""

    def __init__(self):
        # Get idaes version from ver module.  This works even if you just
        # check a new version our from github, have IDAES installed in-place
        # and don't reinstall, which is likely mode for a lot of developers
        self.git_hash = ver.gh
        self.package_version = ver.package_version
        self.version_string = ver.__version__
        # Get Python info
        self.python_version = sys.version
        # Get os info
        self.os_platform = platform.system()
        self.os_release = platform.release()
        self.os_version = platform.version()
        # Get pyomo
        self.pyomo_version = pyomo.version.__version__
        # Get dependency info
        self.dependency_versions = OrderedDict()
        reqs = pkg_resources.get_distribution("idaes-pse").requires()
        for dep in [x.name for x in reqs]:
            if dep == "pyomo":
                continue # pyomo is special
            try:
                self.dependency_versions[dep] = \
                    pkg_resources.get_distribution(dep).version
            except pkg_resources.DistributionNotFound:
                self.dependency_versions[dep] = None
        # Extra packages, users must install these for esoteric features
        self.extra_versions = OrderedDict()
        extras = ["seaborn"]
        for dep in extras:
            try:
                self.extra_versions[dep] = \
                    pkg_resources.get_distribution(dep).version
            except pkg_resources.DistributionNotFound:
                self.extra_versions[dep] = None


    def display_dict(self):
        """Return a dictionary that is in a format that allows easy printing"""
        d = OrderedDict([
            ("IDAES", OrderedDict([
                # for version, I'm breaking off the +label, which is the git hash
                ("IDAES Version", self.version_string.split("+")[0]),
                # and displaying it here
                ("IDAES Git Hash", self.git_hash),
            ])),
            ("Pyomo", OrderedDict([
                ("Pyomo Version", self.pyomo_version)
            ])),
            ("OS", OrderedDict([
                ("Platform", self.os_platform),
                ("Release", self.os_release),
                ("Version", self.os_version),
            ])),
            ("Python", OrderedDict([
                ("Python Version", self.python_version)
            ])),
            ("Dependencies", OrderedDict()),
            ("Extras", OrderedDict()),
        ])
        for k, v in self.dependency_versions.items():
            d["Dependencies"][k] = v if v is not None else "Not Installed"
        for k, v in self.extra_versions.items():
            d["Extras"][k] = v if v is not None else "Not Installed"
        return d
