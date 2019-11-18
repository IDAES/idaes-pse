import pyomo.common.config
import logging.config

_log = logging.getLogger(__name__)

default_config = """
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
                write_config = toml.load(f)
        except IOError:  # don't require config file
            _log.debug("Config file {} not found (this is okay)".format(read_config))
            return
    write_config.set_value(read_config)
    logging.config.dictConfig(write_config["logging"])
    if config_file is not None:
        _log.debug("Read config {}".format(config_file))
