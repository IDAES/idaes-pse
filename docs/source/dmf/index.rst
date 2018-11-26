Data Management Framework
=========================

DMF Contents
------------

.. toctree::
    :maxdepth: 1

    resource

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

For most users the goal of the DMF is to
provide high-level APIs or interface that don't require knowing how the
internal resource representation works. For example, see the
property data indexing API.
