.. _full_install_linux:

Linux installation
==================

This section has the instructions for a "full" Linux installation. If you want to just try a few
examples and find these instructions difficult to follow, you may try the :ref:`min_install_linux`.

System Requirements
-------------------
The IDAES toolkit can be installed on Linux, Windows, or MacOSX. **The officially supported
platform, and the one we use for our automated testing, is Linux.** Therefore it is recommended
that for maximum stability you use this platform. However we realize many users have
Windows or Mac OSX environments. We include best-effort instructions, that we have gotten
to work for us, for those platforms as well.

    * Linux operating system
    * Python 3.6 or above (Python 2 is no longer supported)
    * Basic GNU/C compilation tools: make, gcc/g++
    * ``wget`` (for downloading software)
    * ``git`` (for getting the IDAES source code)
    * Access to the Internet

Things you must know how to do:

    * Get root permissions via `sudo`.
    * Install packages using the package manager.

Installation steps
------------------

.. code-block:: sh

    sudo apt-get install gcc g++ make libboost-dev

We use a Python packaging system called Conda_.
Below are instructions for installing a minimal version of Conda, called Miniconda_.
The full version installs a large number of scientific analysis and visualization libraries
that are not required by the IDAES framework.

.. _Conda: https://conda.io/
.. _Miniconda: https://conda.io/en/latest/miniconda.html

.. code-block:: sh

    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Next, obtain the source code for IDAES from GitHub:

.. code-block:: sh

    git clone https://github.com/IDAES/idaes-pse.git

Download and compile the AMPL Solver Library (ASL) and compile external property
functions; this is required for steam properties and cubic equations of state.
This step is optional, but highly recommended.

.. code-block:: sh

    cd <Location to keep the ASL>
    wget https://ampl.com/netlib/ampl/solvers.tgz
    tar -xf solvers.tgz
    cd solvers
    ./configure
    make
    export ASL_BUILD=`pwd`/sys.`uname -m`.`uname -s`
    cd <IDAES source main directory>
    make

.. note:: If you get an error about ``funcadd.h`` not being found, either ``ASL_BUILD`` is not set correctly or the ASL did not compile properly.

If you are familiar with Python/Conda environments, you will probably
want to create a new environment for your IDAES installation before
starting to install Python and/or Conda packages,
*e.g.*, ``conda create -n <env>`` then ``conda activate <env>``.
If you are not familiar with these commands, don't worry, this is
an optional step.

Install the required Python packages:

.. code-block:: sh

    pip install -r requirements.txt
    python setup.py develop  # or "install"

Install ipopt.  If you have an HSL license, you may prefer to compile ipopt with
HSL support.  Please see the ipopt `documentation <https://projects.coin-or.org/Ipopt>`_
in that case.  Otherwise ipopt can be installed with conda.

.. code-block:: sh

    conda install -c conda-forge ipopt


At this point, you should be able to launch the Jupyter Notebook server and successfully `run examples <examples.html>`_ from the ``examples`` folder:

.. code-block:: sh

    jupyter notebook

Solvers
-------
Some of the model code depends on external solvers. The installation instructions
above include the free IPOPT_ solver. Most of the examples can run with this solver,
but a significant number of more advanced problems will not be handled well. Some
other solvers you can install that may improve (or make possible) solutions for
these models are:

    * CPLEX: a linear optimization package from `IBM <https://www.ibm.com/analytics/cplex-optimizer>`_.
    * Gurobi: LP/MILP/MIQP, etc., solvers from `Gurobi <http://www.gurobi.com>`_.

.. _IPOPT: https://www.coin-or.org/Ipopt/documentation/documentation.html

ASL and AMPL
^^^^^^^^^^^^
In some cases, IDAES uses AMPL user-defined functions written in C for property
models.  Compiling these functions is optional, but some models may not work
without them.

The AMPL solver library (ASL) is required, and can be downloaded from
from https://ampl.com/netlib/ampl/solvers.tgz.  Documentation is available at
https://ampl.com/resources/hooking-your-solver-to-ampl/.
