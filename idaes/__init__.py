""" __init__.py for idaes module

Set up logging for the idaes module, and import plugins.
"""
import os
import logging.config
from importlib import import_module as _do_import
import toml
import pyomo.common.plugin

# Set default configuration.  Used TOML string to serve as an example for
# and definative guid for IDAES configuration files.
_config = toml.loads("""
[plugins]
required = []
optional = []
[logging]
version = 1
disable_existing_loggers = false
formatters.f1.format = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
formatters.f1.datefmt = "%Y-%m-%d %H:%M:%S"
handlers.console.class = "logging.StreamHandler"
handlers.console.formatter = "f1"
handlers.console.stream = "ext://sys.stderr"
loggers.idaes.level = "INFO"
loggers.idaes.handlers = ["console"]
""")

# Use default logging config until config files have been read.
logging.config.dictConfig(_config["logging"])
_log = logging.getLogger(__name__)

# Try to read the global IDAES config file.
try:
    if os.name == 'nt': # Windows
        config_file = os.path.join(os.environ['APPDATA'], '.idaes/idaes.conf')
    else: # any other OS
        config_file = os.path.join(os.environ['HOME'], '.idaes/idaes.conf')
    with open(config_file, "r") as f:
        _config.update(toml.load(f))
    _log.debug("Read global idaes.conf")
except IOError:
    # don't require config file
    _log.debug("Global idaes.conf not found (this is okay)")
except AttributeError:
    _log.debug("Could not find path for global idaes.conf (this is okay)")

# Try working directory for an IDAES config.
try:
    with open("idaes.conf", "r") as f:
        _config.update(toml.load(f))
    _log.debug("Read local idaes.conf")
except IOError:
    # don't require config file
    _log.debug("Local idaes.conf not found (this is okay)")

logging.config.dictConfig(_config["logging"])
_log.debug("'idaes' logger debug test")

##
## Load pugins. This code was taken and condensed from Pyomo
## specifically pyomo.environ.__init__.py
##

def _import_packages(packages, optional=True):
    """Import plugin package
    Args:
        packages: list of pacakges in which to look for plugins
        optional: true, log ImportError but contine; false, raise if ImportError
    """
    for name in packages:
        pname = name + '.plugins' # look in plugins sub-package
        try:
            pkg = _do_import(pname)
        except ImportError as e:
            _log.exception("failed to import plugin: {}".format(pname))
            if not optional:
                raise e
        if hasattr(pkg, 'load'): # run load function for a module if it exists
            pkg.load()

# This make "idaes" the current plugin envirnment while importing plugins here
pyomo.common.plugin.push("idaes") # Add idaes plugin environment at top of stack
# Import plugins standard IDAES plugins, non-optional plugins (currently none)
_import_packages(_config["plugins"]["required"], optional=False)
# Import contrib plugins, failure to import these is non-fatal. (currently none)
_import_packages(_config["plugins"]["optional"], optional=True)
# go back to the previous plugin environment (what it was before pushing "idaes")
pyomo.common.plugin.pop()
