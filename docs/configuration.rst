Configuration
=============

Some behavior of IDAES, especially logging, is configurable through configuration
files. IDAES's configuration is obtained by first setting everything to internal
defaults; then loading a global config file, if it exists; then loading a config
file from the current working directory, if it exists.

Configuration file are in `TOML format <https://github.com/toml-lang/toml>`_. The
default configuration is shown below and can be used as a template to create new
configuration files. This is the configuration used by IDAES if nothing else is
provided.

.. code-block:: yaml

  default_binary_url = "https://github.com/IDAES/idaes-ext/releases/download/1.0.1/"
  use_idaes_solvers = true
  [plugins]
    required = ["idaes"]
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
      propagate = true
      handlers = ["console"]
    [logging.loggers."idaes.init"]
      level = "INFO"
      propagate = false
      handlers = ["console"]
    [logging.loggers."idaes.model"]
      level = "INFO"
      propagate = false
      handlers = ["console"]

Global Configuration
--------------------

IDAES configuration files are named idaes.conf. The easiest way to find where the
global configuration file should be placed is to run the command
``idaes data-directory``.  A global configuration file won't exist unless a user
creates one. The default configuration above can be used as a start.

Windows
~~~~~~~

On Windows the global configuration file is located at
``%LOCALAPPDATA%\idaes\idaes.conf``.

UNIX-Like
~~~~~~~~~

On Unix-like systems the global configuration files is located at
``$HOME/.idaes/idaes.conf``.

Other
~~~~~

On systems that have neither an ``%LOCALAPPDATA%`` or ``$HOME`` environment
variable, global config files are not currently supported.


Important Configuration Entries
-------------------------------

The configuration file has several fields, but they are not all important to
end-users. This section lists the commonly used entries.

logging
~~~~~~~

This section of the file configures IDAES loggers.  Once the configuration is
read, Python's standard ``logging.config.dictConfig()`` is used to set the logger
configuration.  See Python's logging documentation for more information.

IDAES has three main loggers defined in the standard configuration, although
additional loggers can be added if desired.  The standard loggers are:

  1. idaes, this is the root logger of most IDAES logging, unless otherwise noted.

  2. idaes.init, this is the root of IDAES initialization loggers.

  3. idaes.model, this is the root of model loggers.  Model loggers are usually used models written using the IDAES framework, but not part of the ``idaes`` package.

use_idaes_solvers
~~~~~~~~~~~~~~~~~

This option can be set to ``false`` to direct the IDAES framework not to use
solvers obtained with the ``idaes get-extensions`` command.  This can be used if
a user would prefer to use solver versions they have installed apart from IDAES.
