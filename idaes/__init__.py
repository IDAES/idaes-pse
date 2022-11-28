# coding: utf-8
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
"""__init__.py for idaes module

Set up logging for the idaes module, and import plugins.
"""
import os
import copy
import logging


def _handle_optional_compat_activation(
    env_var: str = "IDAES_ACTIVATE_V1_COMPAT",
):
    _log = logging.getLogger("idaes_v1_compat")
    found_in_env = os.environ.get(env_var, None)
    if found_in_env:
        _log.warning(
            "Found environment variable %s=%s. Activating IDAES V1 compatibility.",
            env_var,
            found_in_env,
        )
        try:
            from _idaes_v1_compat import activate
        except ImportError:
            _log.error("Required package _idaes_v1_compat not found")
        else:
            activate()


_handle_optional_compat_activation()


from . import config

from .ver import __version__  # noqa

_log = logging.getLogger(__name__)

# Standard locations for config file, binary libraries and executables, ...
data_directory, bin_directory, testing_directory = config.get_data_directory()
# To avoid a circular import the config module doesn't import idaes, but
# some functions in the config module that are executed later use this
# these directories are static from here on.
config.data_directory = data_directory
config.bin_directory = bin_directory
config.testing_directory = testing_directory

# Set the path for the global and local config files
if data_directory is not None:
    _global_config_file = os.path.join(data_directory, "idaes.conf")
else:
    _global_config_file = None
_local_config_file = "idaes.conf"

# Create the general IDAES configuration block, with default config
cfg = config._new_idaes_config_block()
config.reconfig(cfg)
# read global config and overwrite provided config options
config.read_config(_global_config_file, cfg=cfg)
# read local config and overwrite provided config options
config.read_config(_local_config_file, cfg=cfg)

# Setup the environment so solver executables can be run
config.setup_environment(bin_directory, cfg.use_idaes_solvers)

# Debug log for basic testing of the logging config
_log.debug("'idaes' logger debug test")


def _create_data_dir():
    """Create the IDAES directory to store data files in."""
    config.create_dir(data_directory)


def _create_bin_dir(bd=None):
    """Create the IDAES directory to store executable files in.

    Args:
        bd: alternate binary directory, used for testing
    """
    _create_data_dir()
    if bd is None:
        bd = bin_directory
    config.create_dir(bd)


def _create_testing_dir():
    """Create an idaes testing directory"""
    _create_data_dir()
    config.create_dir(testing_directory)


if data_directory is not None:
    try:
        _create_data_dir()
    except FileNotFoundError:
        pass  # the standard place for this doesn't exist, shouldn't be a show stopper

    try:
        _create_bin_dir()
    except FileNotFoundError:
        pass  # the standard place for this doesn't exist, shouldn't be a show stopper

    try:
        _create_testing_dir()
    except FileNotFoundError:
        pass  # the standard place for this doesn't exist, shouldn't be a show stopper


def reconfig():
    return config.reconfig(cfg)


def read_config(val):
    return config.read_config(val=val, cfg=cfg)


def write_config(path, default=False):
    _cfg = None if default else cfg
    return config.write_config(path=path, cfg=_cfg)


class temporary_config_ctx(object):
    def __enter__(self):
        self.orig_config = copy.deepcopy(cfg)

    def __exit__(self, exc_type, exc_value, traceback):
        global cfg
        cfg = self.orig_config
        reconfig()
