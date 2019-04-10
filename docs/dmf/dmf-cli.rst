
DMF Command-line Interface
==========================

.. contents::
    :depth: 2

.. program:: dmf

dmf
---

Data management framework command wrapper.

dmf options
^^^^^^^^^^^

.. option:: -v
.. option:: --verbose

Increase verbosity. Show warnings if given once, then info, and then
debugging messages.

.. option:: -q
.. option:: --quiet

Increase quietness. If given once, only show critical messages.
If given twice, show no messages.

.. program:: dmf-init

dmf usage
^^^^^^^^^

Run ``sub-command`` with logging at level "error":

.. code-block:: console

    $ dmf <sub-command>

Run ``sub-command`` and log warnings:

.. code-block:: console

    $ dmf <sub-command>

Run ``sub-command`` and log informational / warning messages:

.. code-block:: console

    $ dmf -vv <sub-command>

Run ``sub-command`` only logging fatal errors:

.. code-block:: console

    $ dmf -q <sub-command>

Run ``sub-command`` with no logging at all:

.. code-block:: console

    $ dmf -qq <sub-command>


dmf init
--------
Initialize the current workspace. Optionally, create a new workspace.

dmf init options
^^^^^^^^^^^^^^^^

.. option:: --path path

Use the provided ``path`` as the workspace path.

.. option:: --create

Create a new workspace at location provided by :option:`--path`. Use the
:option:`--name` and :option:`--desc` options to set the workspace name and
description, respectively. If these are not given, they will be prompted for
interactively.

.. option:: --name

Workspace name, used by :option:`--create`

.. option:: --desc

Workspace description, used by :option:`--create`

dmf init usage
^^^^^^^^^^^^^^
This command sets a value in the user-global configuration file
in ``.dmf``, in the user's home directory, so that all other dmf
commands know which workspace to use. With the :option:`--create` option,
a new empty workspace can be created.

.. note:: In the following examples, the current working directory is
          set to ``/home/myuser``.

Create new workspace in sub-directory ``ws``, with given name and description:

.. code-block:: console

    $ dmf init --path ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml

Create new workspace in sub-directory ``ws``, providing the name and
description interactively:

.. code-block:: console

    $ dmf init --path ws --create
    New workspace name: foo
    New workspace description: foo workspace description
    Configuration in '/home/myuser/ws/config.yaml

Switch to workspace ``ws2``:

.. code-block:: console

    $ dmf init --path ws2

If you try to switch to a non-existent workspace, you will get an error message:

.. code-block:: console

    $ dmf init --path doesnotexist
    Existing workspace not found at path='doesnotexist'
    Add --create flag to create a workspace.

If the workspace exists, you cannot create it:

.. code-block:: console

    $ dmf init --path ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml
    $ dmf init --path ws --create
    Cannot create workspace: path 'ws' already exists

dmf status
----------
This command shows basic information about the current active workspace
and, optionally, some additional details. It does not (yet) give any way
to modify the workspace configuration. To do that, you need to edit the
``config.yaml`` file in the workspace root directory.

dmf status options
^^^^^^^^^^^^^^^^^^

.. option:: --color

Allow (if terminal supports it) colored terminal output. This is the default.

.. option:: --no-color

Disallow, even if terminal supports it, colored terminal output.
UNIX output streams to pipes should be detected and have color disabled,
but this option can force that behavior if detection is failing.

.. option:: -s,--show info

Show one of the following types of information:

files
    Count and total size of files in workspace
htmldocs
    Configured paths to the HTML documentation (for "%dmf help" magic in the
    Jupyter Notebook)
logging
    Configuration for logging
all
    Show all items above

.. option:: -a,--all

This option is just an alias for "--show all".

dmf status usage
^^^^^^^^^^^^^^^^

In the following examples, the current working directory is
set to ``/home/myuser`` and the workspace is named ``ws``.

Also note that the output shown below is plain (black) text. This is due to our
limited understanding of how to do colored text in our documentation tool
(Sphinx). In a color-capable terminal, the output will be more colorful.


Show basic workspace status:

.. code-block:: console

    $ dmf status
    settings:
      workspace: /home/myuser/ws
    workspace:
      location: /home/myuser/ws
      name: myws
      description: my workspace
      created: 2019-04-09 12:46:40
      modified: 2019-04-09 12:46:40

Add the file information:

.. code-block:: console

    $ dmf status --show files
    settings:
      workspace: /home/myuser/ws
    workspace:
      location: /home/myuser/ws
      name: myws
      description: my workspace
      created: 2019-04-09 12:52:49
      modified: 2019-04-09 12:52:49
      files:
        count: 3
        total_size: 1.3 MB

You can repeat the :option:`-s,--show` option to add more things:

.. code-block:: console

    $ dmf status --show files --show htmldocs
    settings:
      workspace: /home/myuser/ws
    workspace:
      location: /home/myuser/ws
      name: myws
      description: my workspace
      created: 2019-04-09 12:54:10
      modified: 2019-04-09 12:54:10
      files:
        count: 3
        total_size: 1.3 MB
      html_documentation_paths:
        -: /home/myuser/idaes/docs/build

However, showing everything is less typing, and not overwhelming:

.. code-block:: console

    $ dmf status -a
    settings:
      workspace: /home/myuser/ws
    workspace:
      location: /home/myuser/ws
      name: myws
      description: my workspace
      created: 2019-04-09 12:55:05
      modified: 2019-04-09 12:55:05
      files:
        count: 3
        total_size: 1.3 MB
      html_documentation_paths:
        -: /home/myuser/idaes/docs/build
      logging:
        not configured

dmf ls
------
This command lists resources in the current workspace.

dmf ls options
^^^^^^^^^^^^^^

.. option:: -s,--show

Pick field to show in output table. This field can be repeated to show
any known subset of fields. Default fields, if this option is not
specified at all, are "type", "desc", and "modified". The resource identifier
field is always shown first.

codes
    List name of code(s) in resource. May be shortened with ellipses.
created
    Date created.
desc
    Description of resource.
files
    List names of file(s) in resource. May be shortened with ellipses.
modified
    Date modified.
type
    Name of the type of resource.
version
    Resource version.

.. option:: -S,--sort

Sort by given field; if repeated, combine to make a compound sort key. These
fields are a subset of those in :option:`-s,--show`, with the addition of
``id`` for sorting by the identifier: "id", "type", "desc", "created", "modified",
and/or "version".

.. option:: --no-prefix

By default, shown identifier is the shortest unique prefix, but if you don't
want the identifier shortened, this option will force showing it in full.

.. option:: -r,--reverse

Reverse the order of the sorting given by (or implied by absence of) the
:option:`-S,--sort` option.

dmf ls usage
^^^^^^^^^^^^
In the following examples, the current working directory is
set to ``/home/myuser`` and the workspace is named ``ws``.

Without arguments, show the resources in an arbitrary (though consistent)
order:

.. code-block:: console

    $ dmf ls
    id   type desc  modified
    0b62 data foo13 2019-03-23 17:49:35
    1e41 data foo10 2019-03-23 17:47:53
    6c9a data foo14 2019-03-23 17:51:59
    d3d5 data bar1  2019-03-26 13:07:02
    e780 data foo11 2019-03-23 17:48:11
    eb60 data foo12 2019-03-23 17:49:08

Add a sort key to sort by, e.g. modified date

.. code-block:: console

    $ dmf ls -S modified
    id   type desc  modified
    1e41 data foo10 2019-03-23 17:47:53
    e780 data foo11 2019-03-23 17:48:11
    eb60 data foo12 2019-03-23 17:49:08
    0b62 data foo13 2019-03-23 17:49:35
    6c9a data foo14 2019-03-23 17:51:59
    d3d5 data bar1  2019-03-26 13:07:02


Especially for resources of type "data", showing the first (possibly only) file
that is referred to by the resource is useful:

.. code-block:: console

    $ dmf ls -S modified -s type -s modified -s files
    id   type modified            files
    1e41 data 2019-03-23 17:47:53 foo10
    e780 data 2019-03-23 17:48:11 foo11
    eb60 data 2019-03-23 17:49:08 foo12
    0b62 data 2019-03-23 17:49:35 foo13
    6c9a data 2019-03-23 17:51:59 foo14
    d3d5 data 2019-03-26 13:07:02 bar1

Note that you don't actually have to show a field to sort by it (compare sort
order with results from command above):

.. code-block:: console

    $ dmf ls -S modified -s type -s files
    id   type files
    1e41 data foo10
    e780 data foo11
    eb60 data foo12
    0b62 data foo13
    6c9a data foo14
    d3d5 data bar1

Add ``--no-prefix`` to show the full identifier:

.. code-block:: console

    $ dmf ls -S modified -s type -s files --no-prefix
    id                               type files
    1e41e6ae882b4622ba9043f4135f2143 data foo10
    e7809d25b390453487998e1f1ef0e937 data foo11
    eb606172dde74aa79eea027e7eb6a1b6 data foo12
    0b62d999f0c44b678980d6a5e4f5d37d data foo13
    6c9a85629cb24e9796a2d123e9b03601 data foo14
    d3d5981106ce4d9d8cccd4e86c2cd184 data bar1


.. include:: ../global.rst
