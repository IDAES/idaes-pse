""" __init__.py for idaes module

Set up logging for the idaes module.
"""
import logging.config
import json

config_dict = {
    "version":1,
    "disable_existing_loggers":False,
    "formatters":{
        "f1":{
            "format":"%(asctime)s - %(levelname)s - %(name)s - %(message)s",
            "datefmt":"%Y-%m-%d %H:%M"}},
    "handlers":{
        "console":{
            "class":"logging.StreamHandler",
            "formatter":"f1",
            "stream":"ext://sys.stdout"}},
    "loggers":{
        "idaes":{
            "level":"INFO",
            "handlers":["console"]}}}
logging.config.dictConfig(config_dict)

try: # Also provide or to switch to ini config file or YAML?
    with open("logging.json", "r") as f:
        config_dict = json.load(f)
    logging.config.dictConfig(config_dict)
except IOError:
    pass # don't require a config file
except:
    logging.getLogger("idaes").exception("Error reading logger config file")

logging.getLogger("idaes").debug("Set up 'idaes' logger")


##
## Load pugins. This code was taken and condensed from Pyomo
## specifically pyomo.environ.__init__.py
##
import sys as _sys
from importlib import import_module as _do_import
import pyomo.common.plugin

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
            logging.getLogger(__name__).exception(
                "failed to import plugins: {}".format(pname))
            if not optional:
                raise e
        if hasattr(pkg, 'load'): # run load function for a package if it exists
            pkg.load()

# This make "idaes" the current plugin envirnment while importing plugins here
pyomo.common.plugin.push("idaes") # Add idaes plugin environment at top of stack
# The next line imports plugins that must import successfully. Currently only
# imports from idaes.core, but other locations can be added as needed.
_import_packages(['idaes.core'], optional=False)
# The next line does nothing.  It's a place holder for eventual addition of
# plugins to contrib. With these packages, failure to import is okay.
_import_packages([], optional=True)
# go back to the previous plugin environment (what it was before pushing "idaes")
pyomo.common.plugin.pop()
