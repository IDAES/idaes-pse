Getting Started
===============

.. _IDAES Installation:

Installation
------------
To install the IDAES PSE framework, follow the set of instructions below that are
appropriate for your needs and operating system. If you get stuck, please contact
`idaes-support@idaes.org <idaes-support@idaes.org>`_.

After installing and testing IDAES, it is strongly recommended to do the IDAES tutorials
located on the |examples-site|.

If you expect to develop custom models, we recommend following the
:ref:`advanced user installation<tutorials/advanced_install/index:Advanced User Installation>`.

The OS specific instructions provide information about optionally installing
Miniconda. If you already have a Python installation you prefer, you can skip
the Miniconda install section.

.. note:: IDAES supports Python 3.7 and above.

.. toctree::
    :glob:
    :hidden:

    linux
    windows
    mac_osx
    binaries
    opt_dependencies

.. list-table::
   :header-rows: 1

   * - System
     - Section
   * - Linux
     - :ref:`tutorials/getting_started/linux:Linux Installation Guide`
   * - Windows
     - :ref:`tutorials/getting_started/windows:Windows Installation Guide`
   * - Mac OSX
     - :ref:`Mac Installation<tutorials/getting_started/mac_osx:Mac/OSX Installation Guide>`

.. warning:: If you are using Python for other complex projects, you may want to
            consider using environments of some sort to avoid conflicting
            dependencies.  There are several good options including conda
            environments if you use Anaconda.

.. _min_updating_install:

Updating an existing installation
---------------------------------

When a new version is released, an IDAES installation can be updated without having to remove and reinstall it from scratch.

The following steps describe how to upgrade an existing installation in-place,
assuming that the installation was done using one of the methods described earlier in this section.

.. warning:: If IDAES was installed in a dedicated environment (e.g. a Conda environment, or Python virtual environment), activate the environment before running any of these commands.

1. Open a terminal and verify the currently installed version of IDAES::

    idaes --version

2. Install the upgraded version of the ``idaes-pse`` package using ``pip install``::

    pip install --upgrade idaes-pse

..

    If a newer version of the ``idaes-pse`` package is available, the currently installed version will be removed and replaced by the newest available version.
    Check again the IDAES version to verify that the upgrade was successful::

        idaes --version

3. Run the ``idaes get-extension`` command to install compiled binaries compatible with the newly upgraded IDAES version.  These binaries include solvers and function libraries.  See :ref:`Binary Packages <tutorials/getting_started/binaries:Binary Packages>` for more details.::

    idaes get-extensions

4. Finally, use the ``idaes get-examples`` command to install the most recent version of the IDAES examples compatible with the upgraded IDAES version.

    .. warning:: If the examples target installation directory is not empty, its contents, including examples installed with a previous IDAES version and other files, **will be overwritten without warning**.
        To avoid losing data, **it is strongly recommended that you make a backup copy of any existing examples directory** before proceeding.

..

    After creating a backup copy of the existing examples directory, run::

        idaes get-examples
