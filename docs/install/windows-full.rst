.. _full_install_windows:

Windows Installation
====================
.. note:: Windows is not officially supported at this time.

This is a complete guide to installing the IDAES framework on Windows. 
The :ref:`windows_extras` includes additional information which may be useful.
This guide includes compiling C++ components.  In the future precompiled versions of these 
libraries will be made available, simplifying the installation process.

Tools
-----
Before installing the IDAES software there are a few development tools that need to be installed.
There are alternatives, but an attempt was made to provide the easiest path here.

1. Install a good text editor (Atom, notepad++, spyder, ... whatever you prefer).
2. Install a *git* client from https://git-scm.com/download/win.
   A git client is not necessary for all users, but
   if you are a developer or advanced user, you will likely want it.
3. Install MSYS2. MSYS2 provides a shell which will allow use of Linux style build tools.
   It also provides a convenient package manager (pacman) which allows for easy
   installation of build tools.

    a. Go to https://www.msys2.org/
    #. Download the x86_64 installer
    #. Run the installer (the default options should be okay)
    #. Open the MSYS2 MinGW 64-bit terminal (go to: start menu/MSYS2 64Bit/MSYS2 MinGW 64Bit).
    #. Update the MSYS2 software::

        pacman -Syu

    #. Repeat the previous step until there are no more updates.
    #. Install the build tools and libraries::

        pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-boost unzip patch make

    #. While MinGW does produce Windows native binaries, depending on linking options,
       some DLLs may be required.  Add the MinWG/MSYS2 DLLs to your path.  For example if MSYS2
       was installed in the default location you would probably want to add ``C:\msys64\mingw64\bin``.
       See :ref:`modify_path_env`.

.. note:: In the MSYS2 terminal the directory structure looks different than the
          regular Windows directory structure.
          The Windows C: drive is located at ``/c``.

Install Miniconda
------------------
1. Download Miniconda (https://docs.conda.io/en/latest/miniconda.html)
2. Run the Miniconda installer (default options should be fine)

Get IDAES
---------
The two main options for getting IDAES are to download the files or to clone the repository.
Cloning the repository requires a git client. For core IDAES developers or users who
need to track the latest developments **and** have access to the idaes-dev repo,
replace "idaes-pse" with "idaes-dev."

Option 1: Download from Github
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Most users can download the release files from https://github.com/IDAES/idaes-pse/releases.
The latest development version can be downloaded by  going to https://github.com/IDAES/idaes-pse
and clicking the "Clone or Download" button then clicking on "Download Zip." Unzip the files to a convenient location.

Option 2: Fork and Clone the Repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
For people who are not IDAES core developers but potentially would like to make
contributions to the IDAES project or closely follow IDAES development, the best way
to get the IDAES files is to fork the IDAES repo on Github, then clone the new fork.
To fork the repository sign into your Github account, and go to https://github.com/IDAES/idaes-pse.
Then, click the "Fork" button in the upper righthand corner of the page.

To clone a repository:

1. Open a command window.
2. Go to the directory where you want to create the local repo.
3. Enter the command (replace "Github_Account" with the Github account of the
   fork you wish to clone)::

    git clone https://github.com/Githhub_Account/idaes-pse

4. The clone command should create a new idaes-pse subdirectory with a local repository.

IDAES Location
^^^^^^^^^^^^^^
In the instructions that follow ``idaes_dir`` will refer to the directory containing the IDAES files.

Compiling ASL
-------------
The AMPL Solver Library (ASL) is required to compile some user-defined functions used
in parts of the IDAES framework (mainly some property packages).

1. Open the MSYS2 MinGW 64-bit terminal (go to: start menu/MSYS2 64Bit/MSYS2 MinGW 64Bit).
2. Create a directory for complied source code in a convenient location, which will be
   referred to as ``src`` in these instructions.  For example (obviously change the
   user name and ``/c`` is the location of the C: drive in Windows) ``mkdir /c/Users/jeslick/src``.
3. Go to the source directory (again replace src with the actual directory)::

    cd src

4. Download the ASL and compile the ASL::

    wget https://ampl.com/netlib/ampl/solvers.tgz
    tar -zxvf solvers.tgz
    cd solvers
    ./configure
      make

Compiling IDAES AMPL Function Extensions
----------------------------------------

IDAES uses some additional user defined AMPL functions for various purposes, but
mainly for physical properties.  Before installing IDAES these functions must be
compiled.

1. Open the MSYS2 MinGW 64-bit terminal.

2. Set the ASL_BUILD environment variable (the directory may differ depending on the
   architecture and replace ``.../src`` with the actual location of your src directory)::

    export ASL_BUILD=C:/.../src/solvers/sys.x86_64.MINGW64_NT-10.0

3. Go to the IDAES directory (replace ``/c/idaes_dir`` with the location
   of the IDAES files)::

    cd /c/idaes_dir/idaes_pse/

4. Run: ``make``

If the compile finishes without errors you can proceed to installing IDAES.

Install IDAES
-------------

1. Open the Anaconda Command prompt
2. Create an ``idaes`` environment and activate it (optional)::

    conda create -n idaes "python>=3.6" pip
    conda activate idaes

.. note::
  If you are using a version of conda older than 4.4 the command on Windows to
  activate a conda environment (for example idaes) is ``activate idaes``.

3. Install requirements::

    pip install -r requirements.txt

4. Install IDAES::

    python setup.py develop

5. (Optional) Install IPOPT::

    conda install -c conda-forge ipopt

.. _windows_extras:

Extras
------

Building Documentation
^^^^^^^^^^^^^^^^^^^^^^

Most users do not need to build this documentation, but if necessary you can.  The instructions here use the ``make`` from the MSYS2 installed above.

  1. Open the Anaconda Command prompt, and activate the IDAES environment
  2. Go to the IDAES directory
  3. Go to the docs subdirectory
  4. Add the MSYS2 bin directory to your path temporarily.
     For example, if MSYS2 is installed in the default location::

        set Path=%Path%;C:\msys64\usr\bin

  5. Run make (from MSYS2)::

        make html

The HTML documentation will be in the "build" subdirectory.

Compiling IPOPT
^^^^^^^^^^^^^^^

It's not required to compile IPOPT yourself, and these are pretty much the standard
IPOPT compile instructions.  If you have set up MSYS2 as above, you should be able to
follow these instructions to compile IPOPT for Windows.

1. Download IPOPT from https://www.coin-or.org/download/source/Ipopt/, and put the zip file in the ``src`` directory created above.
2. Open the MSYS2 MinGW 64-bit terminal (go to: start menu/MSYS2 64Bit/MSYS2 MinGW 64Bit).
3. Unzip Ipopt (the ``*`` here represents the portion of the file names with the Ipopt
   version information)::

    unzip Ipopt*.zip
    cd Ipopt*

4. Get third party libraries::

    cd ThirdParty/ASL
    ./get.ASL
    cd ../Blas
    ./get.Blas
    # and so on for all the other subdirectories except HSL.
    # If you have an HSL license follow the instructions in the HSL directory

5. Go to the IPOPT directory (replace $IPOPT_DIR with the IPOPT directory)::

    cd $IPOPT_DIR
    ./configure
    make

6. The IPOPT AMPL executable will be in ./Ipopt/src/Apps/AmplSolver/ipopt.exe, you
   can move the executable to a location in the path (environment variable).
   See :ref:`modify_path_env`.

.. _modify_path_env:

Modifying the Path Environment Variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Windows ``Path`` environment variable provides a search path for executable code
and dynamically linked libraries (DLLs).  You can temporarily modify the path in a
command window session or permanently modify it for the whole system.

**Changing Path Via the Control Panel**

This method will modify the path for the whole system.  Running programs especially
open command windows will need to be restarted for this change to take effect.

A. Any version of Windows

    1. Press the "Windows Key."
    2. Start to type "Control Panel"
    3. Click on "Control Panel" in the start menu.
    4. Click "System and Security."
    5. Click "System."
    6. Click "Advanced system settings."
    7. Click "Environment Variables."

B. In Windows 10

    1. Press the "Windows Key."
    2. Start to type "Environment"
    3. Click on "Edit the system environment" in the start menu.
    4. Click "Environment Variables."

**Temporary Change in Command Window**

This method temporarily changes the path in just the active command window.
Once the command window is closed the change will be lost.

Just set the Path variable to include any additional directories you want to add to
the path.  Replace "added_directory" with the directory you want to add::

    set Path=%Path%;added_directory

