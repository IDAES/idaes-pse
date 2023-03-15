Getting Started
===============

.. _IDAES Installation:

.. note:: IDAES supports Python |python-min| to |python-max|.  Newer versions of Python may work, but are untested. 

Installation
------------
To install the IDAES PSE framework, follow the set of instructions below that are appropriate for 
your needs. The OS specific instructions provide optional steps for installing Miniconda, which can be 
skipped. If you are an IDAES developer or expect to change IDAES code, we recommend following the
:ref:`advanced user installation<tutorials/advanced_install/index:Advanced User Installation>`.
Please contact `idaes-support@idaes.org <idaes-support@idaes.org>`_, if you have difficulty installing 
IDAES.

After installing and testing IDAES, it is recommended that you do IDAES tutorials located on 
the |examples-site|.

OS Specific Instructions
~~~~~~~~~~~~~~~~~~~~~~~~

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

.. warning:: If you are using Python for other projects or installing multiple versions of IDAES, 
            you may want to consider using environments to avoid conflicting dependencies.  There are
            several good options including 
            `conda environments <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ 
            or `venv <https://docs.python.org/3/library/venv.html>`_.

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

3. Run the ``idaes get-extension`` command to install compiled binaries compatible with the newly upgraded IDAES version. The ``--extra petsc`` argument installs the optional PETSc solver. These binaries include solvers and function libraries.  See :ref:`Binary Packages <tutorials/getting_started/binaries:Binary Packages>` for more details.::

    idaes get-extensions --extra petsc

4. Install the IDAES example Jupyter notebooks.

.. include:: install_templates/examples.txt

