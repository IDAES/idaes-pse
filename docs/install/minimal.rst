Minimal installation
====================

To make it easier to use basic functionality and try the IDAES PSE Toolkit,
we have compiled these  "minimal" instructions, that only allow one to use the free 
IPOPT_ solver with `MUMPS`_. This will not be appropriate for some models.
We are working on an easy installer with better
solvers, but for now you will need to use the full install instructions in
the next sections if this is not sufficient for your needs.

.. _IPOPT: https://www.coin-or.org/Ipopt/documentation/documentation.html

.. _MUMPS: http://mumps.enseeiht.fr/

.. _min_install_windows:

Minimal install with IPOPT/MUMPS for Windows
--------------------------------------------

**Install Miniconda** [2]_

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe
2. Install anaconda from the downloaded file in (1).
3. Open the Anaconda powershell (Start -> "Anaconda Powershell Prompt").
4. In the Anaconda Powershell, follow the :ref:`min_install_generic` instructions.

.. _min_install_linux:

Minimal install with IPOPT/MUMPS for Linux
------------------------------------------

**Install  Miniconda** [2]_

1. Download: https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
2. For the next steps, open a terminal window
3. Run the script you downloaded in (1).
4. Follow the :ref:`min_install_generic` instructions.

.. [2] Miniconda is a product from `Anaconda <https://anaconda.com>`_ that contains
       their package manager, "Conda" (and not much else). This is the package manager we
       will use here for setting up the software development environment
       and installing IDAES' software (package) dependencies.


.. _min_install_generic:

Generic minimal install with IPOPT/MUMPS
----------------------------------------

Once you have Conda installed, the remaining steps, performed in either the
Anaconda Powershell (Prompt) or a Linux terminal, are the same.

If you are familiar with Python/Conda environments, you will probably
want to create a new environment for your IDAES installation before
starting. If you do not know what that means, don't worry, this is
an optional step.

**Install a git client**

1. Install the git client::

    conda install -c anaconda git

**Install IPOPT**

2. Install IPOPT from "conda-forge"::

    conda install -c conda-forge ipopt

3. Check if the installation worked by checking for the ipopt version::

    ipopt -v

**Download IDAES source code and install required packages**

4. Go to the idaes-pse releases page, https://github.com/IDAES/idaes-pse/releases/, and
   look at the most recent release. Under the
   section labeled "Assets" there will be a zip file. Download that file and
   extract the contents in any location of your choice.
5. In the Linux terminal or Anaconda Powershell, navigate to the folder you created
   in the previous step.
6. Install the packages required for IDAES using the following command::

    pip install -r requirements.txt

**Install IDAES**

7. In the folder where the idaes source code was downloaded, run the *setup.py* file::

    python setup.py develop

8. Run tests on unit models::

     pytest idaes/unit_models

9. You should see the tests run and all should pass to ensure the installation worked.
    You can report problems on the `Github issues page <https://github.com/IDAES/idaes-pse/issues>`_
    (Please try to be specific about the command and the offending output.)
10. Launch the Jupyter Notebook

    a. Navigate to `examples` and run Jupyter notebook::

            cd idaes/examples
            jupyter notebook

    b. Open a web browser to the URL that is printed from the previous command.

