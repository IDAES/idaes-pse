# coding: utf-8
""" __init__.py for idaes module

Set up logging for the idaes module, and import plugins.
"""
import os
import logging.config
import pyomo.common.plugin
import idaes.config
import toml

from .ver import __version__  # noqa

_log = logging.getLogger(__name__)

# Create the general IDAES configuration block
_config = config.new_idaes_config_block()

# Standard locations for config file, binary libraries and executables, ...
try:
    if os.name == 'nt':  # Windows
        data_directory = os.path.join(os.environ['LOCALAPPDATA'], "idaes")
    else:  # any other OS
        data_directory = os.path.join(os.environ['HOME'], ".idaes")
except AttributeError:
    data_directory = None

# Standard location for executable binaries.
if data_directory is not None:
    bin_directory = os.path.join(data_directory, "bin")
else:
    bin_directory = None

# Standard location for IDAES library files.
if data_directory is not None:
    lib_directory = os.path.join(data_directory, "lib")
else:
    lib_directory = None

def _create_data_dir():
    """Create the IDAES directory to store data files in."""
    config.create_dir(data_directory)

def _create_bin_dir():
    """Create the IDAES directory to store executable files in."""
    _create_data_dir()
    config.create_dir(bin_directory)

def _create_lib_dir():
    """Create the IDAES directory to store library files in."""
    _create_data_dir()
    config.create_dir(lib_directory)

# Set default configuration.  Used TOML string to serve as an example for
# and definitive guide for IDAES configuration files.
config.read_config(toml.loads(config.default_config), _config)

# Try to read the global IDAES config file.
# Set where to look for config files
_global_config_file = os.path.join(data_directory, "idaes.conf")
_local_config_file = "idaes.conf"

# Try to read global config then local
config.read_config(_global_config_file, _config)
config.read_config(_local_config_file, _config)
_log.debug("'idaes' logger debug test")

if _config.use_idaes_solvers:
    # Add IDAES stuff to the path unless you configure otherwise
    os.environ['PATH'] = os.pathsep.join([bin_directory, os.environ['PATH']])
    if os.name == 'nt':  # Windows (this is to find MinGW libs)
        os.environ['PATH'] = os.pathsep.join([os.environ['PATH'], lib_directory])
    else: # Linux and OSX, so far no need for this, but maybe in future
        __orig_ld = os.environ.get('LD_LIBRARY_PATH', '')
        os.environ['LD_LIBRARY_PATH'] = os.pathsep.join(
            [__orig_ld, lib_directory])

# Load plugins, could read a config file later by calling _read_config, but
# plugins only automatiaclly import when 'idaes' is imported. Could call
# _import_packages again later though

# This make "idaes" the current plugin environment while importing plugins here
pyomo.common.plugin.push("idaes")  # Add idaes plugin environment at top of stack
# Import plugins standard IDAES plugins, non-optional plugins
config.import_packages(_config["plugins"]["required"], optional=False)
# Import contrib plugins, failure to import these is non-fatal.
config.import_packages(_config["plugins"]["optional"], optional=True)
# go back to the previous plugin environment (what it was before pushing "idaes")
pyomo.common.plugin.pop()
