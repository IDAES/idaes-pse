##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import pyomo.common.config
import logging.config
import json
import yaml
import os
import importlib

_log = logging.getLogger(__name__)
# Default release version if no options provided for get-extensions
default_binary_release = "2.4.2"
# Where to download releases from get-extensions
release_base_url = "https://github.com/IDAES/idaes-ext/releases/download"
# Where to get release checksums
release_checksum_url = \
    "https://raw.githubusercontent.com/IDAES/idaes-ext/main/releases/sha256sum_{}.txt"
# Map some platform names to others for get-extensions
binary_platform_map = {
    "rhel6": "centos6",
    "rhel7": "centos7",
    "rhel8": "centos8",
    "linux": "centos7",
}
# Set of known platforms with available binaries and descriptions of them
known_binary_platform = {
    "auto":"Auto-select windows, darwin or linux",
    "windows":"Microsoft Windows (built on verion 1909)",
    "darwin": "OSX (currently not available)",
    "linux": "Linux (maps to {})".format(binary_platform_map["linux"]),
    "centos6": "CentOS 6",
    "centos7": "CentOS 7",
    "centos8": "CentOS 8",
    "rhel6": "Red Hat Enterprise Linux 6",
    "rhel7": "Red Hat Enterprise Linux 7",
    "rhel8": "Red Hat Enterprise Linux 8",
    "ubuntu1804": "Ubuntu 18.04",
    "ubuntu1910": "Ubuntu 19.10",
    "ubuntu2004": "Ubuntu 20.04",
}

# Store the original environment variable values so we can revert changes
orig_environ = {
    "PATH": os.environ.get("PATH", ""),
    "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
    "DYLD_LIBRARY_PATH": os.environ.get("DYLD_LIBRARY_PATH", ""),
}

# Default configuration json string.  Store the default config like this because
# it provides a nice example of how to write a config file and test to make sure
# they get read properly.
default_config = """
{
    "use_idaes_solvers":true,
    "default_solver":"ipopt",
    "logger_capture_solver":true,
    "logger_tags":[
        "framework",
        "model",
        "flowsheet",
        "unit",
        "control_volume",
        "properties",
        "reactions"
    ],
    "valid_logger_tags":[
        "framework",
        "model",
        "flowsheet",
        "unit",
        "control_volume",
        "properties",
        "reactions",
        "ui"
    ],
    "ipopt":{
        "options":{
            "nlp_scaling_method":"gradient-based",
            "tol":1e-6
        }
    },
    "logging":{
        "version":1,
        "disable_existing_loggers":false,
        "formatters":{
            "default_format":{
                "format": "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S"
            }
        },
        "handlers":{
            "console":{
                "class": "logging.StreamHandler",
                "formatter": "default_format",
                "stream": "ext://sys.stdout"
            }
        },
        "loggers":{
            "idaes":{
                "level": "INFO",
                "propagate": true,
                "handlers": ["console"]
            },
            "idaes.solve":{
                "propagate": false,
                "level": "INFO",
                "handlers": ["console"]
            },
            "idaes.init":{
                "propagate": false,
                "level": "INFO",
                "handlers": ["console"]
            },
            "idaes.model":{
                "propagate":false,
                "level": "INFO",
                "handlers": ["console"]
            }
        }
    }
}
"""

cfg = None # the idaes ConfigBlock once it's created by new_idaes_config_block()

class ConfigBlockJSONEncoder(json.JSONEncoder):
    """ This class handles non-serializable objects that may appear in the IDAES
    ConfigBlock. For now this is only set objects.
    """
    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        return obj

def _new_idaes_config_block():
    """The idaes configuration is stored in a Pyomo ConfigBlock created by this.
    This function is called when importing IDAES.  This function should only be
    called in ``__init__.py`` for ``idaes``.  Calling it anywhere else will cause
    the idaes configuration system to function improperly.
    """
    global cfg
    cfg = pyomo.common.config.ConfigBlock("idaes", implicit=True)
    cfg.declare(
        "logging",
        pyomo.common.config.ConfigBlock(
            implicit=True,
            description="Logging configuration dictionary",
            doc="This stores the logging configuration. See the Python "
            "logging.config.dictConfig() documentation for details.",
        ),
    )

    cfg.declare(
        "ipopt",
        pyomo.common.config.ConfigBlock(
            implicit=False,
            description="Default config for 'ipopt' solver",
            doc="Default config for 'ipopt-iades' solver"
        ),
    )

    cfg["ipopt"].declare(
        "options",
        pyomo.common.config.ConfigBlock(
            implicit=True,
            description="Default solver options for 'ipopt'",
            doc="Default solver options for 'ipopt' solver"
        ),
    )

    cfg["ipopt"]["options"].declare(
        "nlp_scaling_method",
        pyomo.common.config.ConfigValue(
            domain=str,
            default="gradient-based",
            description="Ipopt NLP scaling method",
            doc="Ipopt NLP scaling method"
        ),
    )

    cfg["ipopt"]["options"].declare(
        "tol",
        pyomo.common.config.ConfigValue(
            domain=float,
            default=1e-6,
            description="Ipopt tol option",
            doc="Ipopt tol option"
        ),
    )

    cfg.declare(
        "default_solver",
        pyomo.common.config.ConfigValue(
            default="ipopt",
            domain=str,
            description="Default solver.  See Pyomo's SolverFactory for detauls.",
            doc="Default solver.  See Pyomo's SolverFactory for detauls.",
        ),
    )

    cfg.declare(
        "use_idaes_solvers",
        pyomo.common.config.ConfigValue(
            default=True,
            domain=bool,
            description="If True, search the IDAES bin directory for executables"
                        " first; otherwise, use IDAES bin directory as last resort.",
            doc="If True the IDAES bin directory will be searched for executables "
                "first, which will result in the solvers installed by IDAES being "
                "used in preference to solvers installed on the machine.  If False, "
                "IDAES will only fall back on solvers installed in the IDAES bin "
                "directory if they are not otherwise available.",
        ),
    )

    cfg.declare(
        "valid_logger_tags",
        pyomo.common.config.ConfigValue(
            default=set(),
            domain=set,
            description="List of valid logger tags",
        ),
    )

    cfg.declare(
        "logger_tags",
        pyomo.common.config.ConfigValue(
            default=set(),
            domain=set,
            description="List of logger tags to allow",
        ),
    )

    cfg.declare(
        "logger_capture_solver",
        pyomo.common.config.ConfigValue(
            default=True,
            description="Solver output captured by logger?",
        ),
    )
    d = json.loads(default_config)
    cfg.set_value(d)
    logging.config.dictConfig(cfg["logging"])
    return cfg


def read_config(read_config):
    """Read either a JSON formatted config file or a configuration dictionary.
    Args:
        read_config: A config file path or dict
    Returns:
        None
    """
    global cfg
    config_file = None
    if read_config is None:
        return
    elif isinstance(read_config, (dict, pyomo.common.config.ConfigBlock)):
        pass  # don't worry this catches ConfigBlock too it seems
    else:
        config_file = read_config
        try:
            with open(config_file, "r") as f:
                read_config = json.load(f)
        except IOError:  # don't require config file
            _log.debug("Config file {} not found (this is okay)".format(read_config))
            return
    cfg.set_value(read_config)
    logging.config.dictConfig(cfg["logging"])
    if config_file is not None:
        _log.debug("Read config {}".format(config_file))


def write_config(path, default=False):
    if default:
        _cd = json.loads(default_config)
        with open(path, 'w') as f:
            json.dump(_cd, f, indent=4)
    else:
        with open(path, 'w') as f:
            json.dump(cfg.value(), f, cls=ConfigBlockJSONEncoder, indent=2)


def reconfig():
    read_config(cfg)
    setup_environment(bin_directory, cfg.use_idaes_solvers)


def create_dir(d):
    """Create a directory if it doesn't exist.

    Args:
        d(str): directory path to create

    Retruns:
        None
    """
    if os.path.exists(d):
        return
    else:
        os.mkdir(d)


def get_data_directory():
    """Return the standard data directory for idaes, based on the OS."""
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
    # Standard location for testing files
    if data_directory is not None:
        testing_directory = os.path.join(data_directory, "testing")
    else:
        testing_directory = None

    return data_directory, bin_directory, testing_directory


def setup_environment(bin_directory, use_idaes_solvers):
    """
    Set environment variables for the IDAES session.

    Args:
        bin_directory: directory to find idaes libraries and executables
        use_idaes_solvers: If true look first in the idaes bin directory for
                           executables if false look last in the idaes bin
                           directory for executables.

    Returns:
        None
    """
    oe = orig_environ
    if use_idaes_solvers:
        os.environ['PATH'] = os.pathsep.join([bin_directory, oe.get('PATH', '')])
    else:
        os.environ['PATH'] = os.pathsep.join([oe.get('PATH', ''), bin_directory])
    if os.name != 'nt':  # If not Windwos set lib search path, Windows uses PATH
        os.environ['LD_LIBRARY_PATH'] = os.pathsep.join(
            [oe.get('LD_LIBRARY_PATH', ''), bin_directory])
        # This is for OSX, but won't hurt other UNIX
        os.environ['DYLD_LIBRARY_PATH'] = os.pathsep.join(
            [oe.get('DYLD_LIBRARY_PATH', ''), bin_directory])
