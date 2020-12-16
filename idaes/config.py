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
import os
import importlib

_log = logging.getLogger(__name__)
default_binary_release = "2.2.2"
binary_platform_map = {
    "rhel6": "centos6",
    "rhel7": "centos7",
    "rhel8": "centos8",
    "linux": "centos7",
}
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

default_config = """
{
    "use_idaes_solvers":true,
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

def new_idaes_config_block():
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
        "use_idaes_solvers",
        pyomo.common.config.ConfigValue(
            default=True,
            domain=bool,
            description="Add the IDAES bin directory to the path.",
            doc="Add the IDAES bin directory to the path such that solvers provided "
            "by IDAES will be used in preference to previously installed solvers.",
        ),
    )

    _config.declare(
        "valid_logger_tags",
        pyomo.common.config.ConfigValue(
            default=set(),
            domain=set,
            description="List of valid logger tags",
        ),
    )

    _config.declare(
        "logger_tags",
        pyomo.common.config.ConfigValue(
            default=set(),
            domain=set,
            description="List of logger tags to allow",
        ),
    )

    _config.declare(
        "logger_capture_solver",
        pyomo.common.config.ConfigValue(
            default=True,
            description="Solver output captured by logger?",
        ),
    )

    d = json.loads(default_config)
    _config.set_value(d)
    logging.config.dictConfig(_config["logging"])
    return _config


def read_config(read_config, write_config):
    """Read either a TOML formatted config file or a configuration dictionary.
    Args:
        config: A config file path or dict
    Returns:
        None
    """
    config_file = None
    if read_config is None:
        return
    elif isinstance(read_config, dict):
        pass  # don't worry this catches ConfigBlock too it seems
    else:
        config_file = read_config
        try:
            with open(config_file, "r") as f:
                write_config = json.load(f)
        except IOError:  # don't require config file
            _log.debug("Config file {} not found (this is okay)".format(read_config))
            return
    write_config.set_value(read_config)
    logging.config.dictConfig(write_config["logging"])
    if config_file is not None:
        _log.debug("Read config {}".format(config_file))


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

    return data_directory, bin_directory


def setup_environment(bin_directory, use_idaes_solvers):
    if use_idaes_solvers:
        # Add IDAES stuff to the path unless you configure otherwise
        os.environ['PATH'] = os.pathsep.join([bin_directory, os.environ['PATH']])
        if os.name != 'nt':  # Windows (this is to find MinGW libs)
            os.environ['LD_LIBRARY_PATH'] = os.pathsep.join(
                [os.environ.get('LD_LIBRARY_PATH', ''), bin_directory]
            )
