# coding: utf-8
"""__init__.py for idaes module

Set up logging for the idaes module, and import plugins.
"""
import os

from . import config as cfg
import logging

from .ver import __version__  # noqa

_log = logging.getLogger(__name__)

_convergence_classes = {}

# Standard locations for config file, binary libraries and executables, ...
data_directory, bin_directory = config.get_data_directory()

# Set the path for the global and local config files
_global_config_file = os.path.join(data_directory, "idaes.conf")
_local_config_file = "idaes.conf"

# Create the general IDAES configuration block, with default config
_config = config.new_idaes_config_block()
# read global config and overwrite provided config options
config.read_config(_global_config_file, _config)
# read local config and overwrite provided config options
config.read_config(_local_config_file, _config)

# Setup the environment so solver executables can be run
config.setup_environment(bin_directory, _config.use_idaes_solvers)

# Debug log for basic testing of the logging config
_log.debug("'idaes' logger debug test")


def _create_data_dir():
    """Create the IDAES directory to store data files in."""
    config.create_dir(data_directory)

def _create_bin_dir():
    """Create the IDAES directory to store executable files in."""
    _create_data_dir()
    config.create_dir(bin_directory)

cfg = _config

def register_convergence_class(name):
    def _register_convergence_class(cls):
        _convergence_classes[name] = ".".join([cls.__module__, cls.__name__])
        return cls
    return _register_convergence_class
