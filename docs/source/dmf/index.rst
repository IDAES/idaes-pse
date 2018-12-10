Data Management Framework
=========================

.. contents::


Overview
--------
The Data Management Framework (DMF) is used to manage all the data needed by the
IDAES framework, including flowsheets, models, and results. It stores
metadata and data in files on disk. The DMF can be accessed through its
Python :term:`API` or command-line interfaces. There is work in progress on adding
graphical interfaces for Jupyter Notebooks and stand-alone desktop apps.

The DMF is designed to allow multiple separate threads of work. These are
organized in ``workspaces``. Inside a given workspace, all the information is
represented by generic containers called ``resources``. A resource could be
a file, a set of files, Python code, or really anything the user wants it
to be. A number of metadata fields are provided, but the user can add
arbitrary other information in a "data" field, and as many files as
they wish in a "datafiles" section.
Resources can be connected to each other with ``relations``. There are
four standard relations defined: derived, contains, uses, version. These may
be expanded in the future.

Here is a very simple graphical illustration of these concepts::

    +-------------------------+  +--------------------+
    | WORKSPACE  1            |  | WORKSPACE  2       |
    |                         |  |                    |
    |  (Resource)  (Resource) |  |                    |   ...
    |     ^            ^      |  |       ...          |
    |     | <relation> |      |  |                    |
    |     v            /      |  |                    |
    |   (Resource) <--^       |  |                    |
    +-------------------------+  +--------------------+


For most users the goal of the DMF is to
provide high-level APIs or interface that don't require knowing the
internal resource representation.

Configuration
-------------

The DMF is configured with an optional global configuration file and a
required per-workspace configuration file. By default the global file is
looked for as ``.dmf`` in the user's home directory. Its main function at the
moment is to set the default workspace directory with the ``workspace``
keyword. For example::

  # global DMF configuration
  workspace: ~/data/workspaces/workspace1

The per-workspace configuration has more options. MORE HERE

Usage
-----

You can use the DMF programmatically by instantiating the Python classes.
For details see the :ref:`idaes.dmf package` documentation.

Jupyter notebook usage
^^^^^^^^^^^^^^^^^^^^^^
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

.. _dmf-help:

DMF help
~~~~~~~~

The IDAES Python interfaces are documented with `Sphinx`_. This includes
automatic translation of the comments and structure of the code into
formatted and hyperlinked HTML pages. The ``%dmf help`` command lets you easily
pull up this documentation for an IDAES module, class, or
object. Below are a couple of examples::

  # Need to initialize the DMF first
  from idaes.dmf import magics
  %dmf init path/to/workspace
  # Get help on a module
  Examples go HERE
  # Get help on a class
  EXAMPLE
  # Get help on an object (will show help for the object's class)
  EXAMPLE

The help pages will open in a new window. The location of the built
documentation that they use is configured in the per-workspace DMF
configuration under the ``htmldocs`` keyword (a default value is filled in
when the DMF is first initialized).


.. _Sphinx: https://www.sphinx-doc.org

Reference
---------
See the :ref:`idaes.dmf package` documentation that is generated
automatically from the source code.
