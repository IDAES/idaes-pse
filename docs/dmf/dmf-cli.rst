
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

.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: ../_images/blue-white-band.png
    :width: 100%
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
.. note:: In the following examples, the current working directory is
          set to ``/home/myuser``.

This command sets a value in the user-global configuration file
in ``.dmf``, in the user's home directory, so that all other dmf
commands know which workspace to use. With the :option:`--create` option,
a new empty workspace can be created.

.. testsetup:: dmf-init

    import pathlib
    import os
    from click.testing import CliRunner
    from idaes.dmf.cli import init
    from idaes.dmf.workspace import Workspace
    from idaes.dmf.dmfbase import DMFConfig
    runner = CliRunner()

Create new workspace in sub-directory ``ws``, with given name and description:

.. code-block:: console

    $ dmf init ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml

.. testcode:: dmf-init
    :hide:

    with runner.isolated_filesystem():
        DMFConfig._filename = str(pathlib.Path('.dmf').absolute())
        result = runner.invoke(init, ['ws', '--create', '--name', 'foo',
            '--desc', 'foo workspace description'])
        assert result.exit_code == 0
        assert (pathlib.Path('ws') / Workspace.WORKSPACE_CONFIG).exists()

Create new workspace in sub-directory ``ws``, providing the name and
description interactively:

.. code-block:: console

    $ dmf init  ws --create
    New workspace name: foo
    New workspace description: foo workspace description
    Configuration in '/home/myuser/ws/config.yaml

.. testcode:: dmf-init
    :hide:

    with runner.isolated_filesystem():
        DMFConfig._filename = str(pathlib.Path('.dmf').absolute())
        result = runner.invoke(init, ['ws', '--create'], input='foo\nfoo desc\n')
        assert result.exit_code == 0
        assert (pathlib.Path('ws') / Workspace.WORKSPACE_CONFIG).exists()

Switch to workspace ``ws2``:

.. code-block:: console

    $ dmf init  ws2

If you try to switch to a non-existent workspace, you will get an error message:

.. code-block:: console

    $ dmf init doesnotexist
    Existing workspace not found at path='doesnotexist'
    Add --create flag to create a workspace.
    $ mkdir some_random_directory
    $ dmf init some_random_directory
    Workspace configuration not found at path='some_random_directory/'

.. testcode:: dmf-init
    :hide:

    with runner.isolated_filesystem():
        DMFConfig._filename = str(pathlib.Path('.dmf').absolute())
        runner.invoke(init, ['ws', '--create'], input='foo\nfoo desc\n')
        result = runner.invoke(init, ['doesnotexist'])
        assert result.exit_code != 0
        assert 'not found' in result.output
        os.mkdir('some_random_directory')
        result = runner.invoke(init, ['some_random_directory'])
        assert result.exit_code != 0

If the workspace exists, you cannot create it:

.. code-block:: console

    $ dmf init ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml
    $ dmf init ws --create
    Cannot create workspace: path 'ws' already exists

.. testcode:: dmf-init
    :hide:

    with runner.isolated_filesystem():
        DMFConfig._filename = str(pathlib.Path('.dmf').absolute())
        runner.invoke(init, ['ws', '--create'], input='foo\nfoo desc\n')
        result = runner.invoke(init, ['ws', '--create'])
        assert result.exit_code != 0
        assert 'exists' in result.output

And, of course, you can't create workspaces anywhere you don't
have permissions to creat directories:

.. code-block:: console

    $ mkdir forbidden
    $ chmod 000 forbidden
    $ dmf init forbidden/ws --create
    Cannot create workspace: path 'forbidden/ws' not accessible

.. testcode:: dmf-init
    :hide:

    with runner.isolated_filesystem():
        DMFConfig._filename = str(pathlib.Path('.dmf').absolute())
        os.mkdir('forbidden')
        os.chmod('forbidden', 0)
        result = runner.invoke(init, ['forbidden/ws', '--create'])
        assert result.exit_code != 0
        os.chmod('forbidden', 0o700)

.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: ../_images/blue-white-band.png
    :width: 100%
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
.. note:: In the following examples, the current working directory is
          set to ``/home/myuser`` and the workspace is named ``ws``.

Also note that the output shown below is plain (black) text. This is due to our
limited understanding of how to do colored text in our documentation tool
(Sphinx). In a color-capable terminal, the output will be more colorful.

.. testsetup:: dmf-status

    from pathlib import Path
    from click.testing import CliRunner
    from idaes.dmf.cli import init, status
    from idaes.dmf.dmfbase import DMFConfig
    runner = CliRunner()

    fsctx = runner.isolated_filesystem()
    fsctx.__enter__()
    DMFConfig._filename = str(Path('.dmf').absolute())
    runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])

.. testcleanup:: dmf-status

    fsctx.__exit__(None, None, None)
    DMFConfig._filename = str(Path('~/.dmf').expanduser())

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

.. testcode:: dmf-status
    :hide:

    result = runner.invoke(status, ['--no-color'])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output

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

.. testcode:: dmf-status
    :hide:

    result = runner.invoke(status, ['--no-color', '--show', 'files'])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    assert "files:" in result.output

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

.. testcode:: dmf-status
    :hide:

    result = runner.invoke(status, ['--no-color', '--show', 'files',
        '--show', 'htmldocs'])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    assert "html" in result.output

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

.. testcode:: dmf-status
    :hide:

    result = runner.invoke(status, ['--no-color', '-a'])
    assert result.exit_code == 0
    assert "settings" in result.output
    assert "name: foo" in result.output
    assert "html" in result.output
    assert "logging:" in result.output

.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: ../_images/blue-white-band.png
    :width: 100%
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. program:: dmf-ls

dmf ls
------
This command lists resources in the current workspace.

dmf ls options
^^^^^^^^^^^^^^

.. option:: --color

Allow (if terminal supports it) colored terminal output. This is the default.

.. option:: --no-color

Disallow, even if terminal supports it, colored terminal output.

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
.. note:: In the following examples, the current working directory is
          set to ``/home/myuser`` and the workspace is named ``ws``.

.. testsetup:: dmf-ls

    from pathlib import Path
    from click.testing import CliRunner
    from idaes.dmf.cli import init, ls, register
    from idaes.dmf.dmfbase import DMFConfig
    runner = CliRunner()

    fsctx = runner.isolated_filesystem()
    fsctx.__enter__()
    DMFConfig._filename = str(Path('.dmf').absolute())
    runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
    files = [f"foo1{n}" for n in range(5)]
    files.append("bar1")
    for f in files:
        with open(f, 'w') as fp:
            fp.write("This is some sample data")
        runner.invoke(register, [f])  # add file to DMF

.. testcleanup:: dmf-ls

    fsctx.__exit__(None, None, None)
    DMFConfig._filename = str(Path('~/.dmf').expanduser())


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

.. testcode:: dmf-ls
    :hide:

    result = runner.invoke(ls, ['--no-color'])
    assert result.exit_code == 0
    output1 = result.output
    result = runner.invoke(ls, ['--no-color'])
    assert result.output == output1

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


.. testcode:: dmf-ls
    :hide:

    result = runner.invoke(ls, ['--no-color', '-S', 'modified'])
    assert result.exit_code == 0
    output1 = result.output
    result = runner.invoke(ls, ['--no-color', '--sort', 'modified'])
    assert result.output == output1


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

.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: ../_images/blue-white-band.png
    :width: 100%
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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


.. testsetup:: dmf-info

    from pathlib import Path
    from click.testing import CliRunner
    from idaes.dmf.cli import init, register, info
    from idaes.dmf.dmfbase import DMFConfig
    runner = CliRunner()

    fsctx = runner.isolated_filesystem()
    fsctx.__enter__()
    DMFConfig._filename = str(Path('.dmf').absolute())
    runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
    filename = "foo.txt"
    with open(filename, 'w') as fp:
        fp.write("This is some sample data")
    result = runner.invoke(register, [filename])
    id_all = result.output.strip()
    id_4 = id_all[:4]


.. testcleanup:: dmf-info

    fsctx.__exit__(None, None, None)
    DMFConfig._filename = str(Path('~/.dmf').expanduser())

.. code-block:: console

    $ dmf info 0b62

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


.. testcode:: dmf-info
    :hide:

    result = runner.invoke(info, ["--no-color", id_4], catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output

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

.. testcode:: dmf-info
    :hide:

    result = runner.invoke(info, ["--no-color", "--format", "json", id_4],
                           catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    out = result.output.strip()
    assert out.startswith("{") and out.endswith("}")
    assert '"relations"' in out

And one more time, in "compact" JSON:

.. code-block:: console

        $ dmf info --format jsonc 0b62
        {"id_": "0b62d999f0c44b678980d6a5e4f5d37d", "type": "data", "aliases": [], "codes": [], "collaborators": [], "created": 1553363375.817961, "modified": 1553363375.817961, "creator": {"name": "dang"}, "data": {}, "datafiles": [{"desc": "foo13", "path": "foo13", "sha1": "feee44ad365b6b1ec75c5621a0ad067371102854", "is_copy": true}], "datafiles_dir": "/home/dang/src/idaes/dangunter/idaes-dev/ws2/files/71d101327d224302aa8875802ed2af52", "desc": "foo13", "relations": [{"predicate": "derived", "identifier": "1e41e6ae882b4622ba9043f4135f2143", "role": "object"}], "sources": [], "tags": [], "version_info": {"created": 1553363375.817961, "version": [0, 0, 0, ""], "name": ""}, "doc_id": 4}

.. testcode:: dmf-info
    :hide:

    result = runner.invoke(info, ["--no-color", "--format", "jsonc", id_4],
                           catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    out = result.output.strip()
    import json
    j = json.loads(out)
    assert len(j['datafiles']) == 1

.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: ../_images/blue-white-band.png
    :width: 100%
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. program:: dmf-register

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

.. option:: --is-subject

If given, reverse the sense of any relation(s) added to the resource so that the
newly created resource is the subject and the existing resource is the object.
Otherwise, the new resource is the object of the relation.

dmf register usage
^^^^^^^^^^^^^^^^^^
.. note:: In the following examples, the current working directory is
          set to ``/home/myuser`` and the workspace is named ``ws``.

Register a new file, which is a CSV data file, and use the ``--info``
option to show the created resource.

.. testsetup:: dmf-register

    from pathlib import Path
    import re, json
    from click.testing import CliRunner
    from idaes.dmf.cli import init, register, info
    from idaes.dmf.dmfbase import DMFConfig
    runner = CliRunner()

    fsctx = runner.isolated_filesystem()
    fsctx.__enter__()
    DMFConfig._filename = str(Path('.dmf').absolute())
    runner.invoke(init, ['ws', '--create', '--name', 'foo', '--desc', 'foo desc'])
    filename = "file.csv"
    with open(filename, 'w') as fp:
        fp.write("index,time,value\n1,0.1,1.0\n2,0.2,1.3\n")


.. testcleanup:: dmf-register

    fsctx.__exit__(None, None, None)
    DMFConfig._filename = str(Path('~/.dmf').expanduser())


.. code-block:: console

    $ printf "index,time,value\n1,0.1,1.0\n2,0.2,1.3\n" > file.csv
    $ dmf reg file.csv --info
                          Resource 117a42287aec4c5ca333e0ff3ac89639
  created
     '2019-04-11 03:58:52'
  creator
     name: dang
  datafiles
     - desc: file.csv
       is_copy: true
       path: file.csv
       sha1: f1171a6442bd6ce22a718a0e6127866740c9b52c
  datafiles_dir
     /home/myuser/ws/files/4db42d92baf3431ab31d4f91ab1a673b
  desc
     file.csv
  doc_id
     1
  id_
     117a42287aec4c5ca333e0ff3ac89639
  modified
     '2019-04-11 03:58:52'
  type
     data
  version
     0.0.0 @ 2019-04-11 03:58:52

.. testcode:: dmf-register
    :hide:

    result = runner.invoke(register, ["file.csv", "--info"], catch_exceptions=False)
    assert result.exit_code == 0
    assert filename in result.output
    assert "version" in result.output

If you try to register (add) the same file twice, it will be an error by default.
You need to add the :option:`--no-unique` option to allow it.

.. code-block:: console

    $ printf "index,time,value\n1,0.1,1.0\n2,0.2,1.3\n" > timeseries.csv
    $ dmf add timeseries.csv
    2315bea239c147e4bc6d2e1838e4101f
    $ dmf add timeseries.csv
    This file is already in 1 resource(s): 2315bea239c147e4bc6d2e1838e4101f
    $ dmf add --no-unique timeseries.csv
    3f95851e4931491b995726f410998491

.. testcode:: dmf-register
    :hide:

    result = runner.invoke(register, ["file.csv",], catch_exceptions=False)
    assert result.exit_code != 0
    result = runner.invoke(register, ["file.csv", "--no-unique"], catch_exceptions=False)
    assert result.exit_code == 0

If you register a file ending in ".json", it will be parsed (unless it is
over 1MB) and, if it passes, registered as type JSON. If the parse fails, it
will be registerd as a generic file *unless* the :option:`--strict` option is
given (with this option, failure to parse will be an error):

.. code-block:: console

    $ echo "totally bogus" > notreally.json
    $ dmf reg notreally.json
    2019-04-12 06:06:47,003 [WARNING] idaes.dmf.resource: File ending in '.json' is not valid JSON: treating as generic file
    d22727c678a1499ab2c5224e2d83d9df
    $ dmf reg --strict notreally.json
    Failed to infer resource: File ending in '.json' is not valid JSON

.. testcode:: dmf-register
    :hide:

    not_json = "notreally.json"
    with open(not_json, "w") as fp:
        fp.write("totally bogus\n")
    result = runner.invoke(register, [not_json], catch_exceptions=False)
    assert result.exit_code == 0
    result = runner.invoke(register, [not_json, "--strict", "--no-unique"], catch_exceptions=False)
    assert result.exit_code != 0

You can explicitly specify the type of the resource with the
:option:`-t,--resource-type` option. In that case, any failure
to validate will be an error. For example, if you say the resource is a Jupyter
Notebook file, and it is not, it will fail. But the same file with type "data"
will be fine:

.. code-block:: console

    $ echo "Ceci n'est pas une notebook" > my.ipynb
    $ dmf reg -t notebook my.ipynb
    Failed to load resource: resource type 'notebook': not valid JSON
    $ dmf reg -t data my.ipynb
    0197a82abab44ecf980d6e42e299b258

.. testcode:: dmf-register
    :hide:

    not_nb = "my.ipynb"
    with open(not_nb, "w") as fp:
        fp.write("foo\n")
    result = runner.invoke(register, [not_nb, '-t', 'notebook'])
    assert result.exit_code != 0
    result = runner.invoke(register, [not_nb, '-t', 'data'])
    assert result.exit_code == 0


You can add links to existing resources with the options :option:`--contained`,
:option:`--derived`, :option:`--used`, and :option:`--prev`. For all of these,
the new resource being registered is the target of the relation and the
option argument is the identifier of an existing resource that is the subject of the
relation.

For example, here we add a "shoebox" resource and then some "shoes" that are contained
in it:

.. code-block:: console

    $ touch shoebox.txt shoes.txt closet.txt
    $ dmf add shoebox.txt
    755374b6503a47a09870dfbdc572e561
    $ dmf add shoes.txt --contained 755374b6503a47a09870dfbdc572e561
    dba0a5dc7d194040ac646bf18ab5eb50
    $ dmf info 7553  # the "shoebox" contains the "shoes"
                                Resource 755374b6503a47a09870dfbdc572e561
      created
         '2019-04-11 20:16:50'
      creator
         name: dang
      datafiles
         - desc: shoebox.txt
           is_copy: true
           path: shoebox.txt
           sha1: da39a3ee5e6b4b0d3255bfef95601890afd80709
      datafiles_dir
         /home/dang/src/idaes/dangunter/idaes-dev/docs/ws/files/7f3ff820676b41689bb32bc325fd2d1b
      desc
         shoebox.txt
      doc_id
         9
      id_
         755374b6503a47a09870dfbdc572e561
      modified
         '2019-04-11 20:16:50'
      relations
         - dba0a5dc7d194040ac646bf18ab5eb50 <--[contains]-- ME
      type
         data
      version
         0.0.0 @ 2019-04-11 20:16:50

    $ dmf info dba0  # the "shoes" are in the "shoebox"
                                Resource dba0a5dc7d194040ac646bf18ab5eb50
      created
         '2019-04-11 20:17:28'
      creator
         name: dang
      datafiles
         - desc: shoes.txt
           is_copy: true
           path: shoes.txt
           sha1: da39a3ee5e6b4b0d3255bfef95601890afd80709
      datafiles_dir
         /home/dang/src/idaes/dangunter/idaes-dev/docs/ws/files/a27f98c24d1848eaba1b26e5ef87be88
      desc
         shoes.txt
      doc_id
         10
      id_
         dba0a5dc7d194040ac646bf18ab5eb50
      modified
         '2019-04-11 20:17:28'
      relations
         - 755374b6503a47a09870dfbdc572e561 --[contains]--> ME
      type
         data
      version
         0.0.0 @ 2019-04-11 20:17:28

.. testcode:: dmf-register
    :hide:

    for text_file in "shoebox", "shoes", "closet":
        open(f"{text_file}.txt", "w")
    result = runner.invoke(register, ["shoebox.txt"], catch_exceptions=False)
    assert result.exit_code == 0
    shoebox_id = result.output.strip()
    result = runner.invoke(register, ["shoes.txt", "--contained", shoebox_id], catch_exceptions=False)
    assert result.exit_code == 0
    shoe_id = result.output.strip()
    result = runner.invoke(info, [shoe_id, "--format", "jsonc"])
    assert result.exit_code == 0

To reverse the sense of the relation, add the :option:`--is-subject` flag.
For example, we now add a "closet" resource that contains the existing "shoebox".
This means the shoebox now has two different "contains" type of relations.

.. code-block:: console

    $ dmf add closet.txt --is-subject --contained 755374b6503a47a09870dfbdc572e561
    22ace0f8ed914fa3ac3e7582748924e4
    $ dmf info 7553
                                Resource 755374b6503a47a09870dfbdc572e561
      created
         '2019-04-11 20:16:50'
      creator
         name: dang
      datafiles
         - desc: shoebox.txt
           is_copy: true
           path: shoebox.txt
           sha1: da39a3ee5e6b4b0d3255bfef95601890afd80709
      datafiles_dir
         /home/dang/src/idaes/dangunter/idaes-dev/docs/ws/files/7f3ff820676b41689bb32bc325fd2d1b
      desc
         shoebox.txt
      doc_id
         9
      id_
         755374b6503a47a09870dfbdc572e561
      modified
         '2019-04-11 20:16:50'
      relations
         - dba0a5dc7d194040ac646bf18ab5eb50 <--[contains]-- ME
         - 22ace0f8ed914fa3ac3e7582748924e4 --[contains]--> ME
      type
         data
      version
         0.0.0 @ 2019-04-11 20:16:50

.. testcode:: dmf-register
    :hide:

    for text_file in "shoebox", "shoes", "closet":
        open(f"{text_file}.txt", "w")
    result = runner.invoke(register, ["closet.txt", "--is-subject",
                           "--contained", shoebox_id], catch_exceptions=False)
    assert result.exit_code == 0
    closet_id = result.output.strip()
    result = runner.invoke(info, [shoebox_id, "--format", "jsonc"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data["relations"]) == 2
    for rel in data["relations"]:
        if rel["role"] == "subject":
            assert rel["identifier"] == shoe_id
        else:
            assert rel["identifier"] == closet_id



.. include:: ../global.rst

