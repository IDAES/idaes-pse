.. _idaes_installation:

Installation
============

.. toctree::
    :hidden:


To install the IDAES PSE framework, follow the set of instructions below that are
appropriate for your needs and operating system. If you get stuck, please contact
`idaes-support@idaes.org <idaes-support@idaes.org>`_.


The OS specific instructions provide information about installing Miniconda. If
you already have a Python installation you prefer, you can skip to the generic
install procedure. Note that IDAES only supports Python 3.6 and above. 

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
4. The IPOPT solver depends on the GNU FORTRAN libraries, which are not bundled
   with it. Unless you know that these are already installed on your system,
   you should manually install them using "apt-get", "yum" or other
   appropriate package manager. Depending on how your system is set up you either
   need to use the "sudo" command or install packages as the "root" user. If these libraries are already installed, you can skip this and proceed with the next step.
   
   apt-get (Debian or Ubuntu based distributions)::

      apt-get install libgfortran3

   yum (RedHat based distributions)::

      yum install libgfortran

5. Follow the :ref:`min_install_generic` instructions.

  

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


**Install IDAES**

1. Install IDAES with pip::

    pip install idaes-pse

2. Run the idaes command to install the compiled binaries::

    idaes get-extensions

   .. note:: The IDAES binary extensions are not yet supported on Mac/OSX

3. Run tests::

    pytest --pyargs idaes -W ignore

4. You should see the tests run and all should pass to ensure the installation worked. You
   may see some "Error" level log messages, but they are okay, and produced by tests for
   error handling. The number of tests that failed and succeeded is reported at the end of the pytest
   output. You can report problems on the `Github issues page <https://github.com/IDAES/idaes-pse/issues>`_
   (Please try to be specific about the command and the offending output.)
