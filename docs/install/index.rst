.. _idaes_installation:

Installation
============

.. toctree::
    :hidden:

    docker

To install the IDAES PSE framework, follow the set of instructions below that are
appropriate for your needs and operating system. If you get stuck, please contact
`idaes-support@idaes.org <idaes-support@idaes.org>`_.

The Docker_ installation works on any platform that supports Docker, but
of course requires installation of, and some understanding of, Docker itself
to operate.

.. _Docker: https://www.docker.com/

The OS specific instructions provide information about installing Miniconda. If
you already have a Python installation you prefer, you can skip to the generic
install procedure.

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
| Docker-based     | :ref:`install_docker`       |
+------------------+-----------------------------+

.. _min_install_windows:

Windows
-------

**Install Miniconda**

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
2. Install anaconda from the downloaded file in (1).
3. Open the Anaconda powershell (Start -> "Anaconda Powershell Prompt").
4. In the Anaconda Powershell, follow the :ref:`min_install_generic` instructions.

.. _min_install_linux:

Linux
-----

**Install  Miniconda**

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
2. Open a terminal window
3. Run the script you downloaded in (1).
4. Follow the :ref:`min_install_generic` instructions.
5. The IPOPT solver depends on the GNU FORTRAN libraries, which are not bundled
   with it. Unless you know that these are already installed on your system,
   you should manually install them using "apt-get", "yum" or other
   appropriate package manager. Depending on how your system is set up you either
   need to use the "sudo" command or install packages as the "root" user.

   apt-get (Debian or Ubuntu based distributions)::

      apt-get install libgfortran3

   yum (RedHat based distributions)::

      yum install libgfortran

.. _min_install_osx:

Mac/OSX
-------

**Install  Miniconda**

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
2. For the next steps, open a terminal window
3. Run the script you downloaded in (1).
4. Follow the :ref:`min_install_generic` instructions.


.. _min_install_generic:

Generic install
---------------

Once you have Conda installed, the remaining steps, performed in either the
Anaconda Powershell (Prompt) or a Linux terminal, are the same.

If you are familiar with Python/Conda environments, you will probably
want to create a new environment for your IDAES installation before
starting to install Python and/or Conda packages,
*e.g.*, ``conda create -n <env> python=3.7`` then ``conda activate <env>``.
If you are not familiar with these commands, don't worry, this is
an optional step.

**Install IDAES**

1. Install IDAES with pip::

    pip install idaes-pse

2. Run the :doc:`idaes get-extensions command <../commands/get_extensions>`
to install the compiled binaries::

    idaes get-extensions

   .. note:: The IDAES binary extensions are not yet supported on Mac/OSX

3. Run the :doc:`idaes get-examples command <../commands/get_examples>` to download
   and install the example files::

    idaes get-examples


4. Run tests::

    pytest --pyargs idaes -W ignore

5. You should see the tests run and all should pass to ensure the installation worked. You
   may see some "Error" level log messages, but they are okay, and produced by tests for
   error handling. The number of failed and succeed test is reported at the end of the pytest
   output. You can report problems on the `Github issues page <https://github.com/IDAES/idaes-pse/issues>`_
   (Please try to be specific about the command and the offending output.)
