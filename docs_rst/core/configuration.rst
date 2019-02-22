IDAES Framework Configuration
=============================

The IDAES framework can be configured with configuration files in TOML format.
Supplying a configuration file is optional. Currently this file sets logging
configuration and modules that should be searched for plugins. The
configuration is done when first importing any idaes.* module. The IDAES
framework will first attempt to read a user-level configuration file at
``%LOCALAPPDATA%\idaes\idaes.conf`` on Windows or ``$HOME/.idaes/idaes.conf`` on
other operating systems (e.g. Linux or Mac).  Next if an idaes.conf file exists
in the working directory it will be read. Configuration files in the working
directory will override settings in the user-level configuration file.  The user
level configuration file will override default settings.  Not all setting need
to be set in a configuration file.

An example configuration file is given below with the default settings.

.. this is toml but ini is probably close enough for now.  toml to recognized
.. code-block:: ini

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
      stream = "ext://sys.stderr"
    [logging.loggers.idaes]
      level = "INFO"
      handlers = ["console"]

The Python dictConfig_ method is used to set up the logger.  The required and
optional elements under plugins are string lists of modules to search for Pyomo
style plugins. Any failure to import plugins in the required modules will raise
an exception, while any failure to import optional plugins will only result in
the exception being logged and execution continuing.

.. _dictConfig: https://docs.python.org/3/library/logging.config.html#logging.config.dictConfig
