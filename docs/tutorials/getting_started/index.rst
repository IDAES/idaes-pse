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

+------------------+-----------------------------+
| System           | Section                     |
+==================+=============================+
| Linux            | :ref:`min_install_linux`    |
+------------------+-----------------------------+
| Windows          | :ref:`min_install_windows`  |
+------------------+-----------------------------+
| Mac OSX          | :ref:`min_install_osx`      |
+------------------+-----------------------------+
| Generic          | :ref:`min_install_generic`  |
+------------------+-----------------------------+

.. warning:: If you are using Python for other complex projects, you may want to
            consider using environments of some sort to avoid conflicting
            dependencies.  There are several good options including conda
            environments if you use Anaconda.


.. _min_install_windows:

Windows
-------

**Install Miniconda (optional)**

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
2. Install anaconda from the downloaded file in (1).
3. Open the Anaconda Prompt (Start -> "Anaconda Prompt").
4. In the Anaconda Prompt, follow the :ref:`min_install_generic` instructions.

.. _min_install_linux:

Linux
-----

**Install  Miniconda (optional)**

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
2. Open a terminal window
3. Run the script you downloaded in (1).

**Install Dependencies**

1. The IPOPT solver depends on the GNU FORTRAN, GOMP, Blas, and Lapack libraries,
   If these libraries are not already installed on your Linux system, you or your
   system administrator can use the sample commands below to install them. If you
   have a Linux distribution that is not listed, IPOPT should still work, but
   the commands to install the required libraries may differ. If these libraries
   are already installed, you can skip this and proceed with the next step.

   .. note:: Depending on your distribution, you may need to prepend ``sudo`` to
            these commands or switch to the "root" user.

   Ubuntu 18.04 and 19.10 and distributions based on them::

      sudo apt-get install libgfortran4 libgomp1 liblapack3 libblas3

   Ubuntu 20.04 and distributions based on it ::

      sudo apt-get install libgfortran5 libgomp1 liblapack3 libblas3

   Current RedHat based distributions, including CentOS::

      yum install lapack blas libgfortran libgomp

**Complete Generic Install**

Follow the :ref:`min_install_generic` instructions.


.. _min_install_osx:

Mac/OSX
-------

**Install  Miniconda (optional)**

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
2. For the next steps, open a terminal window
3. Run the script you downloaded in (1).

**Complete Generic Install**

Follow the :ref:`min_install_generic` instructions.


.. _min_install_generic:

Generic Install
---------------

The remaining steps performed in either the Linux or OSX Terminal or Powershell.
If you installed Miniconda on Windows use the Anaconda Prompt or Anaconda
Powershell Prompt.  Regardless of OS and shell, the following steps are the same.


**Install IDAES**

1. Install IDAES with pip by one of the following methods

  a. To get the latest release::

      pip install idaes-pse

  b. To get a specific release, for example 1.7::

      pip install idaes-pse==1.7

  c. To get the latest version from the GitHub main branch::

      pip install 'idaes-pse[prerelease] @ https://github.com/IDAES/idaes-pse/archive/main.zip'

  d. To get a specific fork or branch, for example myfork (of idaes-pse) and mybranch::

      pip install 'idaes-pse[prerelease] @ https://github.com/myfork/idaes-pse/archive/mybranch.zip'

  e. For developers: follow the :ref:`advanced user installation<tutorials/advanced_install/index:Advanced User Installation>`.

2. Run the :doc:`idaes get-extensions command<../../../reference_guides/commands/get_extensions>`
   to install the compiled binaries. These binaries include solvers and function libraries.
   See :ref:`Binary Packages <tutorials/getting_started/binaries:Binary Packages>`
   for more details.

    idaes get-extensions

..

    .. note:: If you are not able to successfully run the ``idaes get-extensions``
              command due to network security settings or another reason, you can
              download binary release files from
              https://github.com/IDAES/idaes-ext/releases, and extract them in the
              directory indicated by the ``idaes bin-directory`` command. You will
              need both the ``idaes-lib-*`` and ``idaes-solvers-*`` files
              appropriate for your operating system.

..

   .. warning:: The IDAES binary extensions are not yet supported on Mac/OSX.

                As a fallback (assuming you are using a conda env) you can install
                the generic ipopt solver with the command ``conda install -c
                conda-forge ipopt`` though this will not have all the features
                of our extensions package.

3. Run the :doc:`idaes get-examples command <../../../reference_guides/commands/get_examples>` to download
   and install the example files::

    idaes get-examples

..

    By default this will install in a folder "examples" in the current directory.
    The command has many options, but an important
    one is `--dir`, which specifies the folder in which to install.

    for Mac and Linux users this would look like::

        idaes get-examples --dir ~/idaes/examples

    or, for Windows users, it would look like::

        idaes get-examples --dir C:\Users\MyName\IDAES\Examples

    Refer to the full :doc:`idaes get-examples command documentation <../../../reference_guides/commands/get_examples>`
    for more information.

4. Run tests::

    pytest --pyargs idaes -W ignore

5. You should see the tests run and all should pass to ensure the installation worked. You
   may see some "Error" level log messages, but they are okay, and produced by tests for
   error handling. The number of tests that failed and succeeded is reported at the end of the pytest
   output. You can report problems on the |github-issues|
   (Please try to be specific about the command and the offending output.)

**Install IDAES using Conda**

As an alternative to the ``pip install`` method described above, IDAES can also be installed using the Conda package manager.

1. Create a new Conda environment with a name of your choice (in this example, ``my-idaes-env``)::

    conda create --yes --name my-idaes-env python=3.8

..

    .. note:: The ``--yes`` optional flag can be used to perform the installation without having it pausing to ask for confirmation.

2. Activate the ``my-idaes-env`` Conda environment::

    conda activate my-idaes-env

..

    .. note:: This step is needed when starting a new session or a new console tab/window.

3. Install the IDAES Conda package using the ``conda install`` subcommand::

    conda install --yes -c IDAES-PSE -c conda-forge idaes-pse

..

    .. note:: The most recent stable release will be selected by default.
              To instead select a particular version, specify the version tag using an ``=`` after the package name, e.g. ``idaes-pse=1.9.0rc0``.

4. To complete the installation, follow the instructions described in the previous section from Step 2 ("Run the ``idaes get-extensions`` command...") onward.

Optional Dependencies
---------------------
Some tools in IDAES may require additional dependencies. Instructions for installing these dependencies
are located :ref:`here<tutorials/getting_started/opt_dependencies:Optional Dependencies>`.

.. toctree::
    :glob:
    :hidden:

    *

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
