
DMF Command-line Interface
==========================

.. contents::
    :depth: 2

.. program:: dmf

dmf
---

Data management framework command wrapper. This base command has
some options for verbosity that can be applied to any sub-command.

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

.. program:: dmf

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

dmf subcommands
^^^^^^^^^^^^^^^
The subcommands are detailed below. For each, keep in mind that any unique
prefix of that command will be accepted. For example, for ``dmf init``, the
user may also type ``dmf ini``. However, ``dmf in`` will not work because that
would also be a valid prefix for ``dmf info``.

In addition, there are some aliases for some of the sub-commands:

- status => describe
- register => add
- info => resource, show
- ls => list

.. program:: dmf-init

dmf init
--------
Initialize the current workspace. Optionally, create a new workspace.

dmf init options
^^^^^^^^^^^^^^^^

.. option:: path

Use the provided ``path`` as the workspace path. This is required.

.. option:: --create

Create a new workspace at location provided by :option:`path`. Use the
:option:`--name` and :option:`--desc` options to set the workspace name and
description, respectively. If these are not given, they will be prompted for
interactively.

.. option:: --name

Workspace name, used by :option:`--create`

.. option:: --desc

Workspace description, used by :option:`--create`

dmf init usage
^^^^^^^^^^^^^^
In the following examples, the current working directory is
set to ``/home/myuser``.

This command sets a value in the user-global configuration file
in ``.dmf``, in the user's home directory, so that all other dmf
commands know which workspace to use. With the :option:`--create` option,
a new empty workspace can be created.


Create new workspace in sub-directory ``ws``, with given name and description:

.. code-block:: console

    $ dmf init ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml

Create new workspace in sub-directory ``ws``, providing the name and
description interactively:

.. code-block:: console

    $ dmf init  ws --create
    New workspace name: foo
    New workspace description: foo workspace description
    Configuration in '/home/myuser/ws/config.yaml

Switch to workspace ``ws2``:

.. code-block:: console

    $ dmf init  ws2

If you try to switch to a non-existent workspace, you will get an error message:

.. code-block:: console

    $ dmf init doesnotexist
    Existing workspace not found at path='doesnotexist'
    Add --create flag to create a workspace.

If the workspace exists, you cannot create it:

.. code-block:: console

    $ dmf init ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml
    $ dmf init ws --create
    Cannot create workspace: path 'ws' already exists

.. program:: dmf-status

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

.. program:: dmf-ls

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

.. program:: dmf-info

dmf info
--------
Show detailed information about a resource.
This command may also be referred to as ``dmf show``.

dmf info options
^^^^^^^^^^^^^^^^

.. option:: identifier

Identifier, or unique prefix thereof, of the resource.
Any unique prefix of the identifier will work, but if that prefix
matches multiple identifiers, you need to add :option:`--multiple`
to allow multiple records in the output.

.. option:: --multiple

Allow multiple records in the output (see :option:`identifier`)

.. option:: -f,--format value

Output format. Accepts the following values:

term
    Terminal output (colored, if the terminal supports it), with values
    that are empty left out and some values simplified for easy reading.
json
    Raw JSON value for the resource, with newlines and indents for readability.
jsonc
    Raw JSON value for the resource, "compact" version with no extra whitespace
    added.

dmf info usage
^^^^^^^^^^^^^^

The default is to show, with some terminal colors, a summary of the resource:


.. code-block:: console

                              Resource 0b62d999f0c44b678980d6a5e4f5d37d
      created
         '2019-03-23 17:49:35'
      creator
         name: dang
      datafiles
         - desc: foo13
           is_copy: true
           path: foo13
           sha1: feee44ad365b6b1ec75c5621a0ad067371102854
      datafiles_dir
         /home/dang/src/idaes/dangunter/idaes-dev/ws2/files/71d101327d224302aa8875802ed2af52
      desc
         foo13
      doc_id
         4
      id_
         0b62d999f0c44b678980d6a5e4f5d37d
      modified
         '2019-03-23 17:49:35'
      relations
         - 1e41e6ae882b4622ba9043f4135f2143 --[derived]--> ME
      type
         data
      version
         0.0.0 @ 2019-03-23 17:49:35


The same resource in JSON format:

.. code-block:: console

        $ dmf info --format json 0b62
        {
          "id_": "0b62d999f0c44b678980d6a5e4f5d37d",
          "type": "data",
          "aliases": [],
          "codes": [],
          "collaborators": [],
          "created": 1553363375.817961,
          "modified": 1553363375.817961,
          "creator": {
            "name": "dang"
          },
          "data": {},
          "datafiles": [
            {
              "desc": "foo13",
              "path": "foo13",
              "sha1": "feee44ad365b6b1ec75c5621a0ad067371102854",
              "is_copy": true
            }
          ],
          "datafiles_dir": "/home/dang/src/idaes/dangunter/idaes-dev/ws2/files/71d101327d224302aa8875802ed2af52",
          "desc": "foo13",
          "relations": [
            {
              "predicate": "derived",
              "identifier": "1e41e6ae882b4622ba9043f4135f2143",
              "role": "object"
            }
          ],
          "sources": [],
          "tags": [],
          "version_info": {
            "created": 1553363375.817961,
            "version": [
              0,
              0,
              0,
              ""
            ],
            "name": ""
          },
          "doc_id": 4
        }

And one more time, in "compact" JSON:

.. code-block:: console

        $ dmf info --format jsonc 0b62
        {"id_": "0b62d999f0c44b678980d6a5e4f5d37d", "type": "data", "aliases": [], "codes": [], "collaborators": [], "created": 1553363375.817961, "modified": 1553363375.817961, "creator": {"name": "dang"}, "data": {}, "datafiles": [{"desc": "foo13", "path": "foo13", "sha1": "feee44ad365b6b1ec75c5621a0ad067371102854", "is_copy": true}], "datafiles_dir": "/home/dang/src/idaes/dangunter/idaes-dev/ws2/files/71d101327d224302aa8875802ed2af52", "desc": "foo13", "relations": [{"predicate": "derived", "identifier": "1e41e6ae882b4622ba9043f4135f2143", "role": "object"}], "sources": [], "tags": [], "version_info": {"created": 1553363375.817961, "version": [0, 0, 0, ""], "name": ""}, "doc_id": 4}

dmf register
------------
Register a new resource with the DMF, using a file as an input.
An alias for this command is ``dmf add``.

dmf register options
^^^^^^^^^^^^^^^^^^^^

.. option:: --no-copy

Do not copy the file, instead remember path to current location.
Default is to copy the file under the workspace directory.

.. option:: -t,--resource-type

Explicitly specify the type of resource. If this is not given, then
try to infer the resource type from the file. The default will be 'data'.
The full list of resource types is in :py:data:`idaes.dmf.resource.RESOURCE_TYPES`


.. option:: --strict

If inferring the type fails, report an error. With ``--no-strict``, or no option,
if inferring the type fails, fall back to importing as a generic file.

.. option:: --no-unique

Allow duplicate files. The default is ``--unique``, which will
stop and print an error if another resource has a file matching this
file's name and contents.

.. option:: --contained resource

Add a 'contained in' relation to the given resource.

.. option:: --derived resource

Add a 'derived from' relation to the given resource.

.. option:: --used resource

Add a 'used by' relation to the given resource.

.. option:: --prev resource

Add a 'version of previous' relation to the given resource.


dmf register usage
^^^^^^^^^^^^^^^^^^

Register a new file, which is a CSV data file, and use the ``--info``
option to show the created resource.

.. code-block:: console

    $ printf "index,time,value\n1,0.1,1.0\n2,0.2,1.3\n" > file.csv
    $ dmf reg file.csv --info
                          Resource 7fbd197c58374e4abcb6bb0334102ca0
    created
     '2019-04-11 03:51:06'
    creator
     name: dang
    datafiles
     - desc: file.csv
       do_copy: false
       is_copy: false
       path: /home/dang/src/idaes/dangunter/idaes-dev/docs/file.csv
       sha1: f1171a6442bd6ce22a718a0e6127866740c9b52c
    datafiles_dir
     /home/dang/src/idaes/dangunter/idaes-dev/ws2/files/9fa4f1eabe8d457aa7b6a19a23e61f53
    desc
     file.csv
    doc_id
     1
    id_
     7fbd197c58374e4abcb6bb0334102ca0
    modified
     '2019-04-11 03:51:06'
    type
     data
    version
     0.0.0 @ 2019-04-11 03:51:06

    eaa0d83e1ff04200a5359150790fb319



.. include:: ../global.rst

