.. _dmf-overview:

DMF Overview
============
The Data Management Framework (DMF) is used to manage all the data needed by the
IDAES framework, including flowsheets, models, and results. It stores
metadata and data in persistent storage. It does not require that the user
run a server or connect to a remote service. The DMF can be accessed through its
:ref:`API <dmf-api>` or :ref:`command-line interfaces <dmf-cli>`.

.. contents::
    :local:
    :depth: 1

.. _dmf-concepts:

Concepts
--------
The DMF is designed to allow multiple separate threads of work. These are
organized in ``workspaces``. Inside a given workspace, all the information is
represented by containers called ``resources``. A resource describes some
data in the system in a standard way, so it can be searched and manipulated
by the rest of the IDAES framework.
Resources can be connected to each other with ``relations`` such as
"derived", "contains", "uses", and "version".

Below is an illustration of these components.

.. image:: /images/dmf-workspace-resource.png
    :width: 600px


.. _dmf-config:

Configuration
-------------
The DMF is configured with an optional global configuration file and a
required per-workspace configuration file. By default the global file is
looked for as ``.dmf`` in the user's home directory. Its main function at the
moment is to set the default workspace directory with the ``workspace``
keyword. For example::

  # global DMF configuration
  workspace: ~/data/workspaces/workspace1

The per-workspace configuration has more options. See the documentation
in the :class:`Workspace <idaes.dmf.workspace.Workspace>` class for details.
The configuration file is in YAML (or JSON) format. Here is an example file, with some
description in comments:

.. code-block:: yaml

    settings:                               # Global settings
      workspace: /home/myuser/ws            # Path to current workspace
    workspace:                              # Per-workspace settings
      location: /home/myuser/ws             # Path to this workspace
      name: myws                            # Name of this workspace
      description: my workspace             # Description (if any) of this workspace
      created: 2019-04-09 12:55:05          # Date workspace was created
      modified: 2019-04-09 12:55:05         # Date workspace was modified
      files:                                # Basic information about data files
        count: 3                            # How many files
        total_size: 1.3 MB                  # Total size of the files
      html_documentation_paths:             # List of paths for HTML documentation
        -: /home/myuser/idaes/docs/build
      logging:                              # Logging configuration
        idaes.dmf:                          # Name of the logger
            level: DEBUG                    # Log level (Python logging constant)
            output: /tmp/debug.log          # File path or "_stdout_" or "_stderr_"

This configuration file is used whether you use the DMF from the command-line,
Jupyter notebook, or in a Python program. For details see the
:mod:`DMF package <idaes.dmf>` documentation.

Jupyter notebook usage
----------------------
In the Jupyter Notebook, there are some "magics" defined that make
initializing the DMF pretty easy. For example::

    from idaes.dmf import magics
    %dmf init path/to/workspace

The code above loads the "%dmf" *line magic* in the first line, then uses it
to initialize the DMF with the workspace at "path/to/workspace".

From there, other "line magics" will operate in the context of that DMF
workspace.

* ``%dmf help`` - Provide help on IDAES objects and classes. See `dmf-help`_.
* ``%dmf info`` - Provide information about DMF current state for whatever 'topics' are provided
* ``%dmf list`` - List resources in the current workspace
* ``%dmf workspaces`` - List DMF workspaces; you can do this *before* `%dmf init`

.. index::
    pair: dmf;Help

.. _dmf-help:

DMF help
^^^^^^^^

The IDAES Python interfaces are documented with `Sphinx`_. This includes
automatic translation of the comments and structure of the code into
formatted and hyperlinked HTML pages. The ``%dmf help`` command lets you easily
pull up this documentation for an IDAES module, class, or
object. Below are a couple of examples::

    # Initialize the DMF first
    from idaes.dmf import magics
    %dmf init path/to/workspace create

    # Get help on a module (imported)
    from idaes.core import control_volume1d
    %dmf help control_volume1d

    # Get help on a module (by name, no import)
    %dmf help idaes.core.control_volume0d

    # Get help on a class
    from idaes.core.control_volume1d import ControlVolume1DBlock
    %dmf help ControlVolume1DBlock

    # Get help on a class (by name, no import)
    %dmf help idaes.core.control_volume1d.ControlVolume1DBlock

    # Get help on an object (will show help for the object's class)
    # This will end up showing the same help as the previous two examples
    obj = control_volume1d.ControlVolume1DBlock()
    %dmf help obj

The help pages will open in a new window. The location of the built
documentation that they use is configured in the per-workspace DMF
configuration under the ``htmldocs`` keyword (a default value is filled in
when the DMF is first initialized).


.. _Sphinx: https://www.sphinx-doc.org

Sharing
-------

The contents of a DMF workspace can be shared quite simply because
the data is all contained within a directory in the local file system.
So, some ways to share (with one or many people) include:

* Put the workspace directory in a cloud/shared drive like `Dropbox`_ ,
  `Box`_ , `Google Drive`_ , or `OneDrive`_ .
* Put the workspace directory under version control like `Git`_ and
  share that versioned data using Git commands and a service like `Github`_ ,
  `BitBucket`_ or `Gitlab`_.
* Package up the directory with a standard archiving utility like "zip"
  or "tar" and share it like any other file (e.g. attach it to an email).

.. _Box: https://www.box.com/
.. _Dropbox: https://www.dropbox.com/
.. _Google Drive: https://google.com/drive/
.. _OneDrive: https://onedrive.live.com/about/en-us/
.. _Git: https://git-scm.com/
.. _Github: https://github.com/
.. _BitBucket: https://bitbucket.org/
.. _GitLab: https://gitlab.com/

.. note:: These modes of sharing allow users to see the same data, but are not
   designed for real-time collaboration (reading and writing) of the same
   data. That mode of operation requires a proper database server to mediate
   operations on the same data. This is in the roadmap for the DMF, but
   not currently implemented.

Reference
---------
See the :mod:`idaes.dmf` package documentation that is generated
automatically from the source code.
