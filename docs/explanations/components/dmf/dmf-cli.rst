.. _dmf-cli:

DMF Command-line Interface
==========================

This page lists the commands and options for the DMF command-line interface,
which is a Python program called `dmf`. There are several usage examples for each
sub-command. These examples assume the UNIX `bash` shell.

.. contents::
    :local:
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
The subcommands are listed alphabetically below. For each, keep in mind that any unique
prefix of that command will be accepted. For example, for ``dmf init``, the
user may also type ``dmf ini``. However, ``dmf in`` will not work because that
would also be a valid prefix for ``dmf info``.

In addition, there are some aliases for some of the sub-commands:

- ``dmf info`` => `dmf resource` or `dmf show`
- ``dmf ls`` => `dmf list`
- ``dmf register`` => `dmf add`
- ``dmf related`` => `dmf graph`
- ``dmf rm`` => `dmf delete`
- ``dmf status`` => `dmf describe`

usage overview
^^^^^^^^^^^^^^
To give a feel for the context in which you might actually run these
commands, below is a simple example that uses each command:

.. code-block:: console

    # create a new workspace
    $ dmf init ws --name workspace --desc "my workspace" --create
    Configuration in '/home/dang/src/idaes/dangunter/idaes-dev/docs/ws/config.yaml

    # view status of the workspace
    $ dmf status
    settings:
      workspace: /home/myuser/ws
    workspace:
      location: /home/myuser/ws
      name: workspace
      description: my workspace
      created: 2019-04-20 08:32:59
      modified: 2019-04-20 08:32:59

    # add some resources from files
    $ echo "one" > oldfile ; echo "two" > newfile
    $ dmf register oldfile --version 0.0.1
    2792c0ceb0734ed4b302c44884f2d404
    $ dmf register newfile --version 0.0.2 --prev 2792c0ceb0734ed4b302c44884f2d404
    6ddee9bb2bb3420ab10aaf4c74d186f6

    # list the current workspace contents
    $ dmf ls
    id   type desc    modified
    2792 data oldfile 2019-04-20 15:33:11
    6dde data newfile 2019-04-20 15:33:23

    # look at one one resource (newfile)
    $ dmf info 6dde
                                Resource 6ddee9bb2bb3420ab10aaf4c74d186f6
      created
         '2019-04-20 15:33:23'
      creator
         name: dang
      datafiles
         - desc: newfile
           is_copy: true
           path: newfile
           sha1: 7bbef45b3bc70855010e02460717643125c3beca
      datafiles_dir
         /home/myuser/ws/files/8027bf92628f41a0b146a5167d147e9d
      desc
         newfile
      doc_id
         2
      id_
         6ddee9bb2bb3420ab10aaf4c74d186f6
      modified
         '2019-04-20 15:33:23'
      relations
         - 2792c0ceb0734ed4b302c44884f2d404 --[version]--> ME
      type
         data
      version
         0.0.2 @ 2019-04-20 15:33:23

    # see relations
    $ dmf related 2792
    2792 data
        │
        └──┤version├─▶ 6dde data -

    # remove the "old" file
    $ dmf rm 2792
    id                               type desc    modified
    2792c0ceb0734ed4b302c44884f2d404 data oldfile 2019-04-20 15:33:11
    Remove this resource [y/N]? y
    resource removed

    $ dmf ls
    id   type desc    modified
    6dde data newfile 2019-04-20 15:33:23


.. program:: dmf-find

dmf find
--------
Search for resources by a combination of their fields.
Several convenient fields are provided. At this time, a comprehensive
capability to search on any field is not available.

dmf find options
^^^^^^^^^^^^^^^^

In addition to the options below, this command also accepts all the
`dmf ls options`_, although the ``--color/--no-color`` option is
ignored for JSON output.

.. option:: --output value

Output style/format. Possible values:

list
    (Default) Show results as a listing, as from the `ls` subcommand.
info
    Show results as individual records, as from the `info` subcommand.
json
    Show results are JSON objects

.. option:: --by value

Look for "value" in the value of the `creator.name` field.

.. option:: --created value

Use "value" as a date or date range and filter on records that
have a `created` date in that range. Dates should be in the form::

    YYYY-MM-DD[*HH[:MM[:SS[.fff[fff]]]][+HH:MM[:SS[.ffffff]]]]

To indicate a date range, separate two dates with a "..".

* ``2012-03-19``: On March 19, 2012
* ``2012-03-19..2012-03-22``: From March 19 to March 22, 2012
* ``2012-03-19..``: After March 19, 2012
* ``..2012-03-19``: Before March 19, 2012

Note that times may also be part of the date strings.

.. option:: --file value

Look for "value" in the value of the `desc` field in one of the `datafiles`.

.. option:: --modified value

Use "value" as a date or date range and filter on records that
have a `modified` date in that range. See :option:`--created` for
details on the date format.

.. option:: --name value

Look for "value" as one of the values of the `alias` field.

.. option:: --type value

Look for "value" as the value of the `type` field.

dmf find usage
^^^^^^^^^^^^^^

By default, find will essentially provide a filtered listing of
resources. If used without options, it is basically an alias for
`ls`.

.. code-block:: console

    $ dmf ls
    id   type desc      modified
    2517 data file1.txt 2019-04-29 17:29:00
    344c data file2.txt 2019-04-29 17:29:01
    5d98 data A         2019-04-29 17:28:41
    602a data B         2019-04-29 17:28:56
    8c55 data C         2019-04-29 17:28:58
    9cbe data D         2019-04-29 17:28:59
    $ dmf find
    id   type desc      modified
    2517 data file1.txt 2019-04-29 17:29:00
    344c data file2.txt 2019-04-29 17:29:01
    5d98 data A         2019-04-29 17:28:41
    602a data B         2019-04-29 17:28:56
    8c55 data C         2019-04-29 17:28:58
    9cbe data D         2019-04-29 17:28:59

The find-specific options add filters. In the example below, the find
filters for files that were modified after the given date and time.

.. code-block:: console

    $ dmf  find --modified 2019-04-29T17:29:00..
    id   type desc      modified
    2517 data file1.txt 2019-04-29 17:29:00
    344c data file2.txt 2019-04-29 17:29:01

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
    $ mkdir some_random_directory
    $ dmf init some_random_directory
    Workspace configuration not found at path='some_random_directory/'

If the workspace exists, you cannot create it:

.. code-block:: console

    $ dmf init ws --create --name "foo" --desc "foo workspace description"
    Configuration in '/home/myuser/ws/config.yaml
    $ dmf init ws --create
    Cannot create workspace: path 'ws' already exists

And, of course, you can't create workspaces anywhere you don't
have permissions to create directories:

.. code-block:: console

    $ mkdir forbidden
    $ chmod 000 forbidden
    $ dmf init forbidden/ws --create
    Cannot create workspace: path 'forbidden/ws' not accessible

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

Pick field to show in output table. This option can be repeated to show
any known subset of fields. Also the option value can have commas
in it to hold multiple fields. Default fields, if this option is not
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

You can specify other fields from the schema, as long as they are not
arrays of objects, i.e. you can say ``--show tags`` or ``--show version_info.version``,
but ``--show sources`` is too complicated for a tabular listing. To
see detailed values in a record use the `dmf info`_ command.

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

.. option:: -t,--type

Explicitly specify the type of resource. If this is not given, then
try to infer the resource type from the file. The default will be 'data'.
The full list of resource types is in :py:data:`idaes.core.dmf.resource.RESOURCE_TYPES`


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

.. option:: --version

Set the semantic version of the resource.
From 1 to 4 part semantic versions are allowed, e.g.

* `1`
* `1.0`
* `1.0.1`
* `1.0.1-alpha`

See http://semver.org and the function :func:`idaes.core.dmf.resource.version_list` for more details.

dmf register usage
^^^^^^^^^^^^^^^^^^
.. note:: In the following examples, the current working directory is
          set to ``/home/myuser`` and the workspace is named ``ws``.

Register a new file, which is a CSV data file, and use the ``--info``
option to show the created resource.

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

If you register a file ending in ".json", it will be parsed (unless it is
over 1MB) and, if it passes, registered as type JSON. If the parse fails, it
will be registered as a generic file *unless* the :option:`--strict` option is
given (with this option, failure to parse will be an error):

.. code-block:: console

    $ echo "totally bogus" > notreally.json
    $ dmf reg notreally.json
    2019-04-12 06:06:47,003 [WARNING] idaes.core.dmf.resource: File ending in '.json' is not valid JSON: treating as generic file
    d22727c678a1499ab2c5224e2d83d9df
    $ dmf reg --strict notreally.json
    Failed to infer resource: File ending in '.json' is not valid JSON

You can explicitly specify the type of the resource with the
:option:`-t,--type` option. In that case, any failure
to validate will be an error. For example, if you say the resource is a Jupyter
Notebook file, and it is not, it will fail. But the same file with type "data"
will be fine:

.. code-block:: console

    $ echo "Ceci n'est pas une notebook" > my.ipynb
    $ dmf reg -t notebook my.ipynb
    Failed to load resource: resource type 'notebook': not valid JSON
    $ dmf reg -t data my.ipynb
    0197a82abab44ecf980d6e42e299b258

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

You can give your new resource a version with the :option:`--version` option.
You can use this together with the :option:`--prev` option to link
between multiple versions of the same underlying data:

.. code-block:: console

    # note: following command stores the output of "dmf reg", which is the
    #       id of the new resource, in the shell variable "oldid"
    $ oldid=$( dmf reg oldfile.py --type code --version 0.0.1 )
    $ dmf reg newfile.py --type code --version 0.0.2 --prev $oldid
    ef2d801ca29a4a0a8c6f79ee71d3fe07
    $ dmf ls --show type --show version --show codes --sort version
    id   type version codes
    44e7 code 0.0.1   oldfile.py
    ef2d code 0.0.2   newfile.py
    $ dmf related $oldid
    44e7 code
        │
        └──┤version├─▶ ef2d code -


.. program:: dmf-related

dmf related
------------
This command shows resources related to a given resource.

dmf related options
^^^^^^^^^^^^^^^^^^^^

.. option:: -d,--direction

Direction of relationships to show / follow. The possible values are:

in
    Show incoming connection/relationship edges. Since all relations have a
    bi-directional counterpart, this effectively only shows the immediate neighbors
    of the root resource. For example, if the root resource is "A", and "A"
    `contains` "B" and "B" `contains` "C", then this option shows the incoming edge
    from "B" to "A" but not the edge from "C" to "B".

out
    (Default) Show the outgoing connection/relationship edges. This will continue
    until there are no more connections to show, avoiding cycles.
    For example, if the root resource is "A", and "A"
    `contains` "B" and "B" `contains` "C", then this option shows the outgoing edge
    from "A" to "B" and also from "B" to "C".

The default value is ``out``.

.. option:: --color

Allow (if terminal supports it) colored terminal output. This is the default.

.. option:: --no-color

Disallow, even if terminal supports it, colored terminal output.

.. option:: --unicode

Allow unicode drawing characters in the output. This is the default.

.. option:: --no-unicode

Use only ASCII characters in the output.

dmf related usage
^^^^^^^^^^^^^^^^^

In the following examples, we work with 4 resources arranged as a fully
connected square (A, B, C, D). This is not currently possible just with the
command-line, but the following Python code does the job:

.. code-block:: python

    from idaes.core.dmf import DMF, resource
    dmf = DMF()
    rlist = [resource.Resource(value={"desc": ltr, "aliases": [ltr],
                               "tags": ["graph"]})
             for ltr in "ABCD"]
    relation = resource.PR_USES
    for r in rlist:
        for r2 in rlist:
            if r is r2:
                continue
            resource.create_relation_args(r, relation, r2)
    for r in rlist:
        dmf.add(r)

If you save that script as `r4.py`, then the following command-line
actions will run it and verify that everything is created.

.. code-block:: console

    $ python r4.py
    $ dmf ls
    id   type  desc modified
    1e7f other B    2019-04-20 15:43:49
    3bc5 other D    2019-04-20 15:43:49
    ba67 other A    2019-04-20 15:43:49
    f7e9 other C    2019-04-20 15:43:49

You can then see the connections by looking at any one of the
four resource (e.g., `A`):

.. code-block:: console

    $ dmf rel ba67
    ba67 other A
        │
        ├──┤uses├─▶ 3bc5 other D
        ┆  │
        ┆  ├──┤uses├─▶ f7e9 other C
        ┆  │
        ┆  ├──┤uses├─▶ 1e7f other B
        ┆  │
        ┆  └──┤uses├─▶ ba67 other A
        │
        ├──┤uses├─▶ f7e9 other C
        ┆  │
        ┆  ├──┤uses├─▶ 3bc5 other D
        ┆  │
        ┆  ├──┤uses├─▶ 1e7f other B
        ┆  │
        ┆  └──┤uses├─▶ ba67 other A
        │
        └──┤uses├─▶ 1e7f other B
           │
           ├──┤uses├─▶ 3bc5 other D
           │
           ├──┤uses├─▶ f7e9 other C
           │
           └──┤uses├─▶ ba67 other A

If you change the direction of relations, you will get much the same
result, but with the arrows reversed.

.. program:: dmf-rm

dmf rm
-------
Remove one or more resources. This also removes relations (links) to other resources.


dmf rm options
^^^^^^^^^^^^^^

.. option:: identifier

The identifier, or identifier prefix, of the resource(s) to remove

.. option:: --list,--no-list

With the `--list` option, which is the default, the resources to remove, 
or removed, will be listed as if by the ``dmf ls`` command. With 
`--no-list`, then do not produce this output.  

.. option:: -y,--yes

If given, do not confirm removal of the resource(s) with a prompt.
This is useful for scripts that do not want to bother with input,
or people with lots of confidence.

.. option:: --multiple

If given, allow multiple resources to be selected by an identifier prefix. Otherwise,
if the given identifier matches more than one resource, the program will print a message and stop.

dmf rm usage
^^^^^^^^^^^^
.. note:: In the following examples, there are 5 text files named "file1.txt", "file2.txt", .., "file5.txt", in the workspace.
          The identifiers for these files may be different in each example.

Remove one resource, by its full identifier:

.. code-block:: console

    $ dmf ls --no-prefix
    id                               type desc      modified           
    096aa2491e234c4b941f32b537dd3017 data file5.txt 2019-04-16 02:51:30
    821fc8f8e54e4c65b481f483be7f5a2d data file4.txt 2019-04-16 02:51:29
    c20f3a6e338a40ee8a3a4972544adb74 data file1.txt 2019-04-16 02:51:25
    c8f2b5cb80824e649008c414db5287f7 data file3.txt 2019-04-16 02:51:28
    cd62e3bcb9a4459c9f2f5405ca442961 data file2.txt 2019-04-16 02:51:26
    $ dmf rm c20f3a6e338a40ee8a3a4972544adb74
    id                               type desc      modified           
    c20f3a6e338a40ee8a3a4972544adb74 data file1.txt 2019-04-16 02:51:25
    Remove this resource [y/N]? y
    resource removed
    [dmfcli-167 !?]idaes-dev$ dmf ls --no-prefix
    id                               type desc      modified           
    096aa2491e234c4b941f32b537dd3017 data file5.txt 2019-04-16 02:51:30
    821fc8f8e54e4c65b481f483be7f5a2d data file4.txt 2019-04-16 02:51:29
    c8f2b5cb80824e649008c414db5287f7 data file3.txt 2019-04-16 02:51:28
    cd62e3bcb9a4459c9f2f5405ca442961 data file2.txt 2019-04-16 02:51:26

Remove a single resource by its prefix:

.. code-block:: console

    $ dmf ls
    id   type desc      modified           
    6dd5 data file2.txt 2019-04-16 18:51:10
    7953 data file3.txt 2019-04-16 18:51:12
    7a06 data file4.txt 2019-04-16 18:51:13
    e5d7 data file1.txt 2019-04-16 18:51:08
    fe0c data file5.txt 2019-04-16 18:51:15
    $ dmf rm 6d
    id                               type desc      modified           
    6dd57ecc50a24efb824a66109dda0956 data file2.txt 2019-04-16 18:51:10
    Remove this resource [y/N]? y
    resource removed
    $ dmf ls
    id   type desc      modified           
    7953 data file3.txt 2019-04-16 18:51:12
    7a06 data file4.txt 2019-04-16 18:51:13
    e5d7 data file1.txt 2019-04-16 18:51:08
    fe0c data file5.txt 2019-04-16 18:51:15

Remove multiple resources that share a common prefix. In this case, use the
:option:`-y,--yes` option to remove without prompting.

.. code-block:: console

    $ dmf ls
    id   type desc      modified           
    7953 data file3.txt 2019-04-16 18:51:12
    7a06 data file4.txt 2019-04-16 18:51:13
    e5d7 data file1.txt 2019-04-16 18:51:08
    fe0c data file5.txt 2019-04-16 18:51:15
    $ dmf rm --multiple --yes 7
    id                               type desc      modified           
    7953e67db4a543419b9988c52c820b68 data file3.txt 2019-04-16 18:51:12
    7a06435c39b54890a3d01a9eab114314 data file4.txt 2019-04-16 18:51:13
    2 resources removed
    $ dmf ls
    id   type desc      modified           
    e5d7 data file1.txt 2019-04-16 18:51:08
    fe0c data file5.txt 2019-04-16 18:51:15

.. note this is harder to test since we need to force a non-unique
.. prefix. Is it worth it??

.. program:: dmf-status

dmf status
----------
This command shows basic information about the current active workspace
and, optionally, some additional details. It does not (yet) give any way
to modify the workspace configuration. To do that, you need to edit the
``config.yaml`` file in the workspace root directory.
See :ref:`dmf-config`.

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

