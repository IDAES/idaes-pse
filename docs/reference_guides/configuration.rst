Configuration
=============

.. module:: idaes.config
  :noindex:

Some behavior of IDAES is configurable through the IDAES global ConfigBlock.
IDAES's configuration is obtained by first setting everything to internal
defaults; then loading a global config file, if it exists; then loading a config
file from the current working directory, if it exists.  After the ``idaes``
module is imported, the ``idaes`` ConfigBlock can be accessed at ``idaes.cfg``.
Some configuration options can be changed after importing ``idaes`` by calling
``idaes.reconfig()``.

Configuration files are in `JSON format <https://www.json.org/json-en.html>`_.
The default configuration is shown below and can be used as a template to create
new configuration files. To get the IDAES default configuration the command
``idaes config-write --file idaes.conf --default`` can be used.

Global Configuration Files
--------------------------

IDAES configuration files are named ``idaes.conf``. The easiest way to find
where the global configuration file should be placed is to run the command
``idaes data-directory``.  A global configuration file won't exist unless a
user creates one. The default configuration above can be used as a start.

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

Local Configuration Files
-------------------------

Local configuration files are also named ``idaes.conf`` and can be placed in the
working directory, which is the directory you launch Python from.  You can also
use the Python command ``chdir()`` to change the working directory before
importing ``idaes``.

In addition to reading local configuration files when ``idaes`` is imported, you
can read a configuration file anytime by calling ``idaes.read_config(path)``.
Reading a configuration file will automatically apply any resulting
configuration changes.

Changing the Configuration in a Script or Interactive Session
-------------------------------------------------------------

The idaes configuration can be changed anytime after the ``idaes`` module is
imported.  The standard ConfigBlock options are described in detail below.  For
example to change whether you want to use the solvers provided by idaes or ones
you have installed elsewhere, you would first use the command
``idaes.cfg["use_idaes_solvers"] = False`` then to make the change take effect
use ``idaes.reconfig()``.  Not all option changes require ``idaes.reconfig()``,
so whether they do or don't is provided in the options descriptions below.

Important Configuration Entries
-------------------------------

The ConfigBlock has several options, but they are not all important to
end-users. This section lists the commonly used entries.

warning_to_exception
~~~~~~~~~~~~~~~~~~~~

If this option is True, any log messages at level warning or above will be
converted to a RuntimeError exception.  This can be used to ensure a model doesn't
generate warnings. It can also be used to generate a traceback for warnings, which
can provide additional debugging information.

Changes require ``idaes.reconfig()``.  The default setting is ``False``.

deprecation_to_exception
~~~~~~~~~~~~~~~~~~~~~~~~

If this option is True, any log messages at level warning or above that contain
"deprecation," "deprecate," or "deprecated" will be converted to a RuntimeError
exception.  This can be used to ensure a model doesn't use any deprecated models
or methods.

Changes require ``idaes.reconfig()``.  The default setting is ``False``.

use_idaes_solvers
~~~~~~~~~~~~~~~~~

This option can be set to ``False`` (``false`` in JSON) to direct the IDAES
framework not to use solvers obtained with the ``idaes get-extensions`` command
before using the solvers that may have been otherwise installed by the user.
This can be used if a user would prefer to use solver versions they have
installed apart from IDAES.

Changes require ``idaes.reconfig()``.  The default setting is ``True``.

logger_capture_solver
~~~~~~~~~~~~~~~~~~~~~

If a solver call is done from inside a solver logging context, this setting will
send the solver output to the logger if ``True``, and not capture the solver output
for the logger if ``False``.  If solver output is not captured it will be sent to
the screen, and not be logged.

Changes do not require ``idaes.reconfig()``.  The default setting is ``True``.

logger_tags
~~~~~~~~~~~

Loggers created with the ``idaes.logging`` module can be assigned tags.  Output
from these loggers is recorded if the loggers tag is in the ``logger_tags`` set.
The default behavior can be configured in a configuration file. The tag set can
also be modified at any time via functions in the ``idaes.logging`` module. This
is a subset of ``valid_log_tags``.

Changes do not require ``idaes.reconfig()``.  The default setting is:
``["framework", "model", "flowsheet", "unit", "control_volume", "properties", "reactions"]``.

valid_log_tags
~~~~~~~~~~~~~~

When setting logger tags for ``idaes.logging`` loggers they are compared against
a list of valid tags.  This is done to guard against spelling errors. If the
default set of defined tags is not sufficient tags can be added.

Changes do not require ``idaes.reconfig()``.  The default setting is:
``["framework", "model", "flowsheet", "unit", "control_volume", "properties", "reactions", "ui"]``.

ipopt
~~~~~
This is a config block that provides default configuration for the ``ipopt`` solver.
These options are used for ipopt solvers by default when the IDAES SolverFactory
wrapper is used. Currently only solver options can be configured in the ``options``
sub-ConfigBlock.

For example to set the default NLP scaling method for ipopt to use idaes-provided
scaling factors, use the command
``idaes.cfg["ipopt"]["options"]["nlp_scaling_method"] = "user-scaling"``

Any ipopt solver options that can be passed via command line argument to the ipopt
AMPL executable solver can be set under ``idaes.cfg["ipopt"]["options"]``
or equivalently in a configuration file.

Changes do not require ``idaes.reconfig()``.  The default options are:
``{"nlp_scaling_method": "gradient-based"}``.
