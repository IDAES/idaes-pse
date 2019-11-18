# coding: utf-8
""" __init__.py for idaes module

Set up logging for the idaes module, and import plugins.
"""
import os
import logging.config
import importlib
import toml
import pyomo.common.plugin
import pyomo.common.config

from .ver import __version__  # noqa

_log = logging.getLogger(__name__)

# Create the general IDAES configuration block
_config = pyomo.common.config.ConfigBlock("idaes", implicit=False)
_config.declare(
    "logging",
    pyomo.common.config.ConfigBlock(
        implicit=True,
        description="Logging configuration dictionary",
        doc="This stores the logging configuration. See the Python "
        "logging.config.dictConfig() documentation for details.",
    ),
)
_config.declare(
    "plugins",
    pyomo.common.config.ConfigBlock(
        implicit=False,
        description="Plugin search configuration",
        doc="Plugin search configuration",
    ),
)
_config.plugins.declare(
    "required",
    pyomo.common.config.ConfigValue(
        default=[],
        description="Modules with required plugins",
        doc="This is a string list of modules from which to load plugins. "
        "This will look in {module}.plugins for things to load. Exceptions"
        "raised while attempting to load these plugins are considered fatal. "
        "This is used for core plugins.",
    ),
)
_config.plugins.declare(
    "optional",
    pyomo.common.config.ConfigValue(
        default=[],
        description="Modules with optional plugins to load",
        doc="This is a string list of modules from which to load plugins. "
        "This will look in {module}.plugins for things to load. Exceptions "
        "raised while attempting to load these plugins will be logged but "
        "are nonfatal. This is used for contrib plugins.",
    ),
)

_config.declare(
    "use_idaes_solvers",
    pyomo.common.config.ConfigValue(
        default=True,
        description="Add the IDAES bin directory to the path.",
        doc="Add the IDAES bin directory to the path such that solvers provided "
        "by IDAES will be used in preference to previously installed solvers.",
    ),
)

# Standard locations for config file, binary libraries and executables, ...
try:
    if os.name == 'nt':  # Windows
        data_directory = os.path.join(os.environ['LOCALAPPDATA'], "idaes")
    else:  # any other OS
        data_directory = os.path.join(os.environ['HOME'], ".idaes")
except AttributeError:
    data_directory = None

# Standard location for executable binaries.
if data_directory != None:
    bin_directory = os.path.join(data_directory, "bin")
else:
    bin_directory = None

# Standard location for IDAES library files.
if data_directory != None:
    lib_directory = os.path.join(data_directory, "lib")
else:
    lib_directory = None

def _create_data_dir():
    """Create the IDAES directory to store data files in."""
    if os.path.exists(data_directory):
        return
    else:
        os.mkdir(data_directory)

def _create_bin_dir():
    """Create the IDAES directory to store executable files in."""
    _create_data_dir()
    if os.path.exists(bin_directory):
        return
    else:
        os.mkdir(bin_directory)

def _create_lib_dir():
    """Create the IDAES directory to store library files in."""
    _create_data_dir()
    if os.path.exists(lib_directory):
        return
    else:
        os.mkdir(lib_directory)

# Could create directories here, but I'll make that happen when the user does
# something that requires them.  For now the commandline utility commends will
# cause the directories to be made.

def _read_config(config):
    """Read either a TOML formatted config file or a configuration dictionary.
    Args:
        config: A config file path or dict
    Returns:
        None
    """
    config_file = None
    if config is None:
        return
    elif isinstance(config, dict):
        pass  # don't worry this catches ConfigBlock too it seems
    else:
        config_file = config
        try:
            with open(config_file, "r") as f:
                config = toml.load(f)
        except IOError:  # don't require config file
            _log.debug("Config file {} not found (this is okay)".format(config))
            return
    _config.set_value(config)
    logging.config.dictConfig(_config["logging"])
    if config_file is not None:
        _log.debug("Read config {}".format(config_file))


def _import_packages(packages, optional=True):
    """Import plugin package, condensed from pyomo.environ.__init__.py
    Args:
        packages: list of pacakges in which to look for plugins
        optional: true, log ImportError but continue; false, raise if ImportError
    Returns:
        None
    """
    for name in packages:
        pname = name + '.plugins'  # look in plugins sub-package
        try:
            pkg = importlib.import_module(pname)
        except ImportError as e:
            _log.exception("failed to import plugin: {}".format(pname))
            if not optional:
                raise e
        if hasattr(pkg, 'load'):  # run load function for a module if it exists
            pkg.load()


# Set default configuration.  Used TOML string to serve as an example for
# and definitive guide for IDAES configuration files.
_read_config(
    toml.loads(
        """
use_idaes_solvers = true
[plugins]
  required = []
  optional = []
[logging]
  version = 1
  disable_existing_loggers = false
  [logging.formatters.f1]
    format = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"
  [logging.handlers.console]
    class = "logging.StreamHandler"
    formatter = "f1"
    stream = "ext://sys.stdout"
  [logging.loggers.idaes]
    level = "INFO"
    handlers = ["console"]
"""
    )
)

# Try to read the global IDAES config file.
# Set where to look for config files
_global_config_file = os.path.join(data_directory, "idaes.conf")
_local_config_file = "idaes.conf"

# Try to read global config then local
_read_config(_global_config_file)
_read_config(_local_config_file)
_log.debug("'idaes' logger debug test")

if _config["use_idaes_solvers"]:
    # Add IDAES stuff to the path unless you configure otherwise
    os.environ['PATH'] = os.pathsep.join([bin_directory, os.environ['PATH']])
    if os.name == 'nt':  # Windows (this is to find MinGW libs)
        os.environ['PATH'] = os.pathsep.join([os.environ['PATH'], lib_directory])
    else: # Linux and OSX, so far no need for this, but maybe in future
        os.environ['LD_LIBRARY_PATH'] = os.pathsep.join(
            [os.environ['LD_LIBRARY_PATH'], lib_directory])

# Load plugins, could read a config file later by calling _read_config, but
# plugins only automatiaclly import when 'idaes' is imported. Could call
# _import_packages again later though

# This make "idaes" the current plugin environment while importing plugins here
pyomo.common.plugin.push("idaes")  # Add idaes plugin environment at top of stack
# Import plugins standard IDAES plugins, non-optional plugins
_import_packages(_config["plugins"]["required"], optional=False)
# Import contrib plugins, failure to import these is non-fatal.
_import_packages(_config["plugins"]["optional"], optional=True)
# go back to the previous plugin environment (what it was before pushing "idaes")
pyomo.common.plugin.pop()
