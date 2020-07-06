Configuration
=============

Some behavior of IDAES, especially logging, is configurable through configuration
files. IDAES's configuration is obtained by first setting everything to internal
defaults; then loading a global config file, if it exists; then loading a config
file from the current working directory, if it exists.

Configuration file are in `JSON format <https://www.json.org/json-en.html>`_. The
default configuration is shown below and can be used as a template to create new
configuration files. This is the configuration used by IDAES if nothing else is
provided.

.. code-block:: json

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
          "reactions"
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


Global Configuration
--------------------

IDAES configuration files are named ``idaes.conf``. The easiest way to find where the
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

IDAES has four main loggers defined in the standard configuration, although
additional loggers can be added if desired.  The standard loggers are:

  1. idaes, this is the root logger of most IDAES logging, unless otherwise noted.

  2. idaes.init, this is the root of IDAES initialization loggers.

  3. idaes.solve, this is the root of IDAES solver loggers and solver information.

  4. idaes.model, this is the root of model loggers.  Model loggers are
     usually used by models written using the IDAES framework, but not
     part of the ``idaes`` package.

use_idaes_solvers
~~~~~~~~~~~~~~~~~

This option can be set to ``false`` to direct the IDAES framework not to use
solvers obtained with the ``idaes get-extensions`` command.  This can be used if
a user would prefer to use solver versions they have installed apart from IDAES.

logger_capture_solver
~~~~~~~~~~~~~~~~~~~~~

If a solver call is done from inside a solver logging context, this setting will
send the solver output to the logger if ``true``, and not capture the solver output
for the logger if ``false``.  If solver output is not captured it will be sent to
the screen, and not be logged.

logger_tags
~~~~~~~~~~~

Loggers created with the ``idaes.logging`` module can be assigned tags.  Output
from these loggers is recorded if the loggers tag is in the ``logger_tags`` set.
The default behavior can be configured in a configuration file. The tag set can
also be modified at any time via functions in the ``idaes.logging`` module.


valid_log_tags
~~~~~~~~~~~~~~

When setting logger tags for ``idaes.logging`` loggers they are compared against
a list of valid tags.  This is done to guard against spelling errors. If the default
set of defined tags is not sufficient tags can be added here, or later through
functions in the ``idaes.logging`` module.
