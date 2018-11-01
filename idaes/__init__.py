""" __init__.py for idaes module

Set up logging for the idaes module, and import plugins.
"""
import logging.config
import json
from importlib import import_module as _do_import

import pyomo.common.plugin

#TODO<jce> read global idaes config if available
#          the location of the config would probably be:
#          * %APPDATA%/.idaes/idaes.config (Windows)
#          * $HOME/.idaes/idaes.config (not Windows)

_logging_config_dict = {
    "version":1,
    "disable_existing_loggers":False,
    "formatters":{
        "f1":{
            "format":"%(asctime)s - %(levelname)s - %(name)s - %(message)s",
            "datefmt":"%Y-%m-%d %H:%M:%S"}},
    "handlers":{
        "console":{
            "class":"logging.StreamHandler",
            "formatter":"f1",
            "stream":"ext://sys.stderr"}},
    "loggers":{
        "idaes":{
            "level":"INFO",
            "handlers":["console"]}}}
logging.config.dictConfig(_logging_config_dict)

try: # Also provide or to switch to ini config file or YAML?
    with open("logging.json", "r") as f:
        _logging_config_dict = json.load(f)
    logging.config.dictConfig(_logging_config_dict)
except IOError:
    pass # don't require a config file

_log = logging.getLogger(__name__)
_log.debug("Set up 'idaes' logger")


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
        pname = name+'.plugins' # look in plugins sub-package
        try:
            pkg = _do_import(pname)
        except ImportError as e:
            _log.exception("failed to import plugin: {}".format(pname))
            if not optional:
                raise e
        if hasattr(pkg, 'load'): # run load function for a package if it exists
            pkg.load()

# This make "idaes" the current plugin envirnment while importing plugins here
pyomo.common.plugin.push("idaes") # Add idaes plugin environment at top of stack
# Import plugins standard IDAES plugins, non-optional plugins (currently none)
_import_packages([], optional=False)
# Import contrib plugins, failure to import these is non-fatal. (currently none)
_import_packages([], optional=True)
# go back to the previous plugin environment (what it was before pushing "idaes")
pyomo.common.plugin.pop()
