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

__author__ = "John Eslick"

import sys
import platform
import pkg_resources
import json

import idaes
import pyomo
import pyomo.environ as pyo

import idaes.ver as ver


class EnvironmentInfo:
    """Get information about IDAES and the environment IDAES is running in,
    including OS and Critial dependency versions."""

    known_solvers = [
        "ipopt",
        "ipopt_sens",
        "ipopt_l1",
        "bonmin",
        "couenne",
        "cbc",
        "k_aug",
        "dot_sens",
    ]

    extras = [
        "seaborn",
    ]

    def __init__(self, additional_solvers=()):
        # Get idaes version from ver module.  This works even if you just
        # check a new version our from github, have IDAES installed in-place
        # and don't reinstall, which is likely mode for a lot of developers
        self.git_hash = ver.gh
        self.package_version = ver.package_version
        self.version_string = ver.__version__
        self.bin_directory = idaes.bin_directory
        self.data_directory = idaes.data_directory
        self.global_config = idaes._global_config_file
        self.local_config = idaes._local_config_file
        # Get Python info
        self.python_version = sys.version
        self.python_executable = sys.executable
        # Get os info
        self.os_platform = platform.system()
        self.os_release = platform.release()
        self.os_version = platform.version()
        # Get pyomo
        self.pyomo_version = pyomo.version.__version__
        # Get dependency info
        self.dependency_versions = {}
        reqs = pkg_resources.get_distribution("idaes-pse").requires()
        for dep in [x.name for x in reqs]:
            if dep == "pyomo":
                continue  # pyomo is special
            try:
                self.dependency_versions[dep] = pkg_resources.get_distribution(
                    dep
                ).version
            except pkg_resources.DistributionNotFound:
                self.dependency_versions[dep] = None
        # Extra packages, users must install these for esoteric features
        self.extra_versions = {}
        for dep in self.extras:
            try:
                self.extra_versions[dep] = pkg_resources.get_distribution(dep).version
            except pkg_resources.DistributionNotFound:
                self.extra_versions[dep] = None
        self.solver_versions = {}
        for s in self.known_solvers + list(additional_solvers):
            slv = pyo.SolverFactory(s, validate=False)
            try:
                a = slv.available()
            except AttributeError:
                a = False
            if not a:
                self.solver_versions[s] = None
            else:
                v = slv.version()
                if v is None:
                    v = "Unknown Version Installed"
                self.solver_versions[s] = ".".join([f"{x}" for x in v])

    def to_json(self, fname=None):
        if fname is not None:
            with open(fname, "w") as f:
                json.dump(self.to_dict(), f, indent=4)
        else:
            return json.dumps(self.to_dict(), indent=4)

    def to_dict(self):
        """Return a dictionary that is in a format that allows easy printing"""
        d = {
            "IDAES": {
                # break off the +label, which is the git hash
                "Version": self.version_string.split("+")[0],
                "Git Hash": self.git_hash,
                "Binary Directory": self.bin_directory,
                "Data Directory": self.data_directory,
                "Global Config": self.global_config,
            },
            "Pyomo": {
                "Version": self.pyomo_version,
            },
            "Python": {
                "Version": self.python_version,
                "Executable": self.python_executable,
            },
            "OS": {
                "Platform": self.os_platform,
                "Release": self.os_release,
                "Version": self.os_version,
            },
            "Dependencies": {
                k: v if v is not None else "Not Installed"
                for k, v in sorted(self.dependency_versions.items())
            },
            "Extras": {
                k: v if v is not None else "Not Installed"
                for k, v in sorted(self.extra_versions.items())
            },
            "Solvers": {
                k: v if v is not None else "Not Installed"
                for k, v in sorted(self.solver_versions.items())
            },
        }
        return d
