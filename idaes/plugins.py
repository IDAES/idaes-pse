import importlib
import logging

_log = logging.getLogger(__name__)

def import_packages(packages, optional=True):
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
