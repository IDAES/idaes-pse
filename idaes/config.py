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
import pyomo.common.config
import logging.config
import json
import os

_log = logging.getLogger(__name__)
# Default release version if no options provided for get-extensions
default_binary_release = "3.1.0"
# Where to download releases from get-extensions
release_base_url = "https://github.com/IDAES/idaes-ext/releases/download"
# Where to get release checksums
release_checksum_url = (
    "https://raw.githubusercontent.com/IDAES/idaes-ext/main/releases/sha256sum_{}.txt"
)
# This is a list of platforms with builds
base_platforms = (
    "darwin-aarch64",
    "el7-x86_64",
    "el8-x86_64",
    "el8-aarch64",
    "ubuntu1804-x86_64",
    "ubuntu1804-aarch64",
    "ubuntu2004-x86_64",
    "ubuntu2004-aarch64",
    "ubuntu2204-x86_64",
    "ubuntu2204-aarch64",
    "windows-x86_64",
)
# Map some platform names to others for get-extensions
binary_distro_map = {
    "rhel7": "el7",
    "rhel8": "el8",
    "scientific7": "el7",
    "centos7": "el7",
    "centos8": "el8",
    "rocky8": "el8",
    "almalinux8": "el8",
    "debian9": "el7",
    "debian10": "el8",
    "debian11": "ubuntu2004",
    "linuxmint20": "ubuntu2004",
    "kubuntu1804": "ubuntu1804",
    "kubuntu2004": "ubuntu2004",
    "kubuntu2204": "ubuntu2204",
    "xubuntu1804": "ubuntu1804",
    "xubuntu2004": "ubuntu2004",
    "xubuntu2204": "ubuntu2204",
}
# Machine map
binary_arch_map = {
    "x64": "x86_64",
    "intel64": "x86_64",
    "amd64": "x86_64",
    "arm64": "aarch64",
}
# Set of extra binary packages and basic build platforom where available
extra_binaries = {
    "petsc": base_platforms,
}
# Store the original environment variable values so we can revert changes
orig_environ = {
    "PATH": os.environ.get("PATH", ""),
    "LD_LIBRARY_PATH": os.environ.get("LD_LIBRARY_PATH", ""),
    "DYLD_LIBRARY_PATH": os.environ.get("DYLD_LIBRARY_PATH", ""),
}
# Define set of default units of measurement for common derived quantities
default_uom = {
    "pressure": "Pa",
    "energy": "J",
    "energy_mol": "J/mol",
    "energy_mass": "J/kg",
    "entropy": "J/K",
    "entropy_mol": "J/mol/K",
    "entropy_mass": "J/kg/K",
    "power": "W",
}


def canonical_arch(arch):
    """Get the offical machine type in {x86_64, aarch64} if possible, otherwise
    just return arch.lower().

    Args:
        arch (str): machine type string usually from platform.machine()

    Returns (str):
        Canonical machine type used by the binary package names.
    """
    arch = arch.lower()
    return binary_arch_map.get(arch, arch)


def canonical_distro(dist):
    """Get the offical distro name if possible, otherwise just return
    dist.lower(). Distro is used loosely here and includes Windows, Darwin
    (macOS), and other OSs in addition to Linux.

    Args:
        arch (str): machine type string usually from platform.machine()

    Returns (str):
        Canonical machine type used by the binary package names.
    """
    dist = dist.lower()
    return binary_distro_map.get(dist, dist)


class ConfigBlockJSONEncoder(json.JSONEncoder):
    """This class handles non-serializable objects that may appear in the IDAES
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
    cfg = pyomo.common.config.ConfigBlock("idaes", implicit=True)
    cfg.declare(
        "warning_to_exception",
        pyomo.common.config.ConfigValue(
            domain=bool,
            default=False,
            description="Convert any logged warnings or errors to exceptions",
            doc="Convert any logged warnings or errors to exceptions",
        ),
    )
    cfg.declare(
        "deprecation_to_exception",
        pyomo.common.config.ConfigValue(
            domain=bool,
            default=False,
            description="Convert any logged deprecation warnings to exceptions",
            doc="Convert any logged deprecation warnings to exceptions",
        ),
    )
    cfg.declare(
        "logging",
        pyomo.common.config.ConfigBlock(
            implicit=True,
            description="Logging configuration dictionary",
            doc="This stores the logging configuration. See the Python "
            "logging.config.dictConfig() documentation for details.",
        ),
    )
    cfg["logging"].declare(
        "version", pyomo.common.config.ConfigValue(domain=int, default=1)
    )
    cfg["logging"].declare(
        "disable_existing_loggers",
        pyomo.common.config.ConfigValue(domain=bool, default=False),
    )
    cfg["logging"].declare("formatters", pyomo.common.config.ConfigBlock(implicit=True))
    cfg["logging"]["formatters"].declare(
        "default_format", pyomo.common.config.ConfigBlock(implicit=True)
    )
    cfg["logging"]["formatters"]["default_format"].declare(
        "format",
        pyomo.common.config.ConfigValue(
            domain=str, default="%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        ),
    )
    cfg["logging"]["formatters"]["default_format"].declare(
        "datefmt",
        pyomo.common.config.ConfigValue(domain=str, default="%Y-%m-%d %H:%M:%S"),
    )
    cfg["logging"].declare("handlers", pyomo.common.config.ConfigBlock(implicit=True))
    cfg["logging"]["handlers"].declare(
        "console", pyomo.common.config.ConfigBlock(implicit=True)
    )
    cfg["logging"]["handlers"]["console"].declare(
        "class",
        pyomo.common.config.ConfigValue(domain=str, default="logging.StreamHandler"),
    )
    cfg["logging"]["handlers"]["console"].declare(
        "formatter",
        pyomo.common.config.ConfigValue(domain=str, default="default_format"),
    )
    cfg["logging"]["handlers"]["console"].declare(
        "stream",
        pyomo.common.config.ConfigValue(domain=str, default="ext://sys.stdout"),
    )
    cfg["logging"].declare(
        "loggers",
        pyomo.common.config.ConfigValue(
            domain=dict,
            default={  # the period in the logger names is trouble for ConfigBlock
                "idaes": {"level": "INFO", "propagate": True, "handlers": ["console"]},
                "idaes.solve": {"propagate": False, "handlers": ["console"]},
                "idaes.init": {"propagate": False, "handlers": ["console"]},
                "idaes.model": {"propagate": False, "handlers": ["console"]},
            },
        ),
    )
    cfg.declare(
        "ipopt",
        pyomo.common.config.ConfigBlock(
            implicit=False,
            description="Default config for 'ipopt' solver",
            doc="Default config for 'ipopt' solver",
        ),
    )

    cfg["ipopt"].declare(
        "options",
        pyomo.common.config.ConfigBlock(
            implicit=True,
            description="Default solver options for 'ipopt'",
            doc="Default solver options for 'ipopt' solver",
        ),
    )

    cfg["ipopt"]["options"].declare(
        "nlp_scaling_method",
        pyomo.common.config.ConfigValue(
            domain=str,
            default="gradient-based",
            description="Ipopt NLP scaling method",
            doc="Ipopt NLP scaling method",
        ),
    )

    cfg["ipopt"]["options"].declare(
        "tol",
        pyomo.common.config.ConfigValue(
            domain=float,
            default=1e-6,
            description="Ipopt tol option",
            doc="Ipopt tol option",
        ),
    )

    cfg.declare(
        "petsc_ts",
        pyomo.common.config.ConfigBlock(
            implicit=False,
            description="Default config for 'petsc_ts' solver",
            doc="Default config for 'petsc_ts' solver",
        ),
    )

    cfg["petsc_ts"].declare(
        "options",
        pyomo.common.config.ConfigBlock(
            implicit=True,
            description="Default solver options for 'petsc_ts'",
            doc="Default solver options for 'petsc_ts' solver",
        ),
    )

    cfg["petsc_ts"]["options"].declare(
        "--ts_save_trajectory",
        pyomo.common.config.ConfigValue(
            domain=int,
            default=1,
            description="Save the trajectory data from PETSc",
            doc="Save the trajectory data from PETSc",
        ),
    )

    cfg["petsc_ts"]["options"].declare(
        "--ts_max_snes_failures",
        pyomo.common.config.ConfigValue(
            domain=int,
            default=200,
            description="Number of nonliner solver failures before giving up",
            doc="Number of nonliner solver failures before giving up",
        ),
    )

    cfg["petsc_ts"]["options"].declare(
        "--ts_max_reject",
        pyomo.common.config.ConfigValue(
            domain=int,
            default=20,
            description="Maximum number of time steps to reject",
            doc="Maximum number of time steps to reject",
        ),
    )

    cfg["petsc_ts"]["options"].declare(
        "--ts_type",
        pyomo.common.config.ConfigValue(
            domain=str,
            default="beuler",
            description="TS solver to use",
            doc="TS solver to use",
        ),
    )

    cfg["petsc_ts"]["options"].declare(
        "--ts_adapt_type",
        pyomo.common.config.ConfigValue(
            domain=str,
            default="basic",
            description="TS adaptive step size method to use",
            doc="TS adaptive step size method to use",
        ),
    )

    cfg["petsc_ts"]["options"].declare(
        "--ts_exact_final_time",
        pyomo.common.config.ConfigValue(
            domain=str,
            default="matchstep",
            description="How to handle the final time step",
            doc="How to handle the final time step",
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
        "use_idaes_solver_config",
        pyomo.common.config.ConfigValue(
            default=False,
            domain=bool,
            description="If True, use the configure IDAES solver default options.",
            doc="If True, use the configure IDAES solver default options.",
        ),
    )

    cfg.declare(
        "valid_logger_tags",
        pyomo.common.config.ConfigValue(
            default=set(
                [
                    "framework",
                    "model",
                    "flowsheet",
                    "unit",
                    "control_volume",
                    "properties",
                    "reactions",
                    "ui",
                ]
            ),
            domain=set,
            description="List of valid logger tags",
        ),
    )

    cfg.declare(
        "logger_tags",
        pyomo.common.config.ConfigValue(
            default=set(
                [
                    "framework",
                    "model",
                    "flowsheet",
                    "unit",
                    "control_volume",
                    "properties",
                    "reactions",
                ]
            ),
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

    cfg.declare(
        "reporting_units",
        pyomo.common.config.ConfigValue(
            default=default_uom,
            description="Units of measurement for reporting quantities",
        ),
    )
    return cfg


def read_config(val, cfg):
    """Read either a JSON formatted config file or a configuration dictionary.
    Args:
        val: dict, ConfigDict, or file to read
        cfg: config block
    Returns:
        None
    """
    config_file = None
    if val is None:
        return
    elif isinstance(val, (dict, pyomo.common.config.ConfigBlock)):
        pass  # don't worry this catches ConfigBlock too it seems
    else:
        config_file = val
        try:
            with open(config_file, "r") as f:
                val = json.load(f)
        except IOError:  # don't require config file
            _log.debug("Config file {} not found (this is okay)".format(config_file))
            return
    cfg.set_value(val)
    if config_file is not None:
        _log.debug("Read config {}".format(config_file))
    reconfig(cfg)


def write_config(path, cfg=None):
    if cfg is None:
        cfg = _new_idaes_config_block()
    with open(path, "w") as f:
        json.dump(cfg.value(), f, cls=ConfigBlockJSONEncoder, indent=4)


class _WarningToExceptionFilter(logging.Filter):
    """Filter applied to IDAES loggers returned by this module."""

    @staticmethod
    def filter(record):
        if record.levelno >= logging.WARNING:
            raise RuntimeError(f"Logged Warning converted to exception: {record.msg}")


class _DeprecationToExceptionFilter(logging.Filter):
    """Filter applied to IDAES loggers returned by this module."""

    @staticmethod
    def filter(record):
        if record.levelno >= logging.WARNING:
            if "deprecat" in record.msg.lower():
                raise RuntimeError(
                    f"Logged deprecation converted to exception: {record.msg}"
                )


def reconfig(cfg):
    logging.config.dictConfig(cfg.logging.value())
    _log = logging.getLogger("idaes")
    if cfg.deprecation_to_exception:
        _log.addFilter(_DeprecationToExceptionFilter)
    else:
        _log.removeFilter(_DeprecationToExceptionFilter)
    if cfg.warning_to_exception:
        _log.addFilter(_WarningToExceptionFilter)
    else:
        _log.removeFilter(_WarningToExceptionFilter)
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
    if "IDAES_DATA" in os.environ:
        data_directory = os.environ["IDAES_DATA"]
    else:
        try:
            if os.name == "nt":  # Windows
                data_directory = os.path.join(os.environ["LOCALAPPDATA"], "idaes")
            else:  # any other OS
                data_directory = os.path.join(os.environ["HOME"], ".idaes")
        except AttributeError:
            data_directory = None
    if data_directory is None or not os.path.isdir(os.path.dirname(data_directory)):
        _log.warning(
            "IDAES data directory location does not exist: "
            f"{os.path.dirname(data_directory)}"
        )
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
    if bin_directory is None:
        return
    oe = orig_environ
    if use_idaes_solvers:
        os.environ["PATH"] = os.pathsep.join([bin_directory, oe.get("PATH", "")])
    else:
        os.environ["PATH"] = os.pathsep.join([oe.get("PATH", ""), bin_directory])
    if os.name != "nt":  # If not Windows set lib search path, Windows uses PATH
        os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(
            [oe.get("LD_LIBRARY_PATH", ""), bin_directory]
        )
        # This is for macOS, but won't hurt other UNIX
        os.environ["DYLD_LIBRARY_PATH"] = os.pathsep.join(
            [oe.get("DYLD_LIBRARY_PATH", ""), bin_directory]
        )
