Linux Installation Guide
========================

Quickstart
----------
The following commands should be sufficient to get you started with installing and using IDAES. For more information on these steps, see below.

.. sidebar:: Learn more about the Installation
   :class: quickstart-toc

   * `Install Prerequisites`_
   * `Install IDAES-PSE`_
   * `Install IDAES Extensions`_
   * `Install IDAES Examples`_
   * `Run IDAES Tests`_

.. Import quick start guide, including OS specific lines & skipping non-OS lines & comments
.. literalinclude:: install_templates/quickstart.txt
   :language: bash

------------------------------------------------

Installing IDAES
----------------
To get IDAES fully set up on your machine, we'll go through the steps to get idaes-pse package installed as well as setting up the IDAES extensions, which includes some extra solvers and function libraries, the IDAES example files, and the IDAES tests.

.. include:: install_templates/idaes_data.txt

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install  Miniconda**

1. Download `Miniconda <https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh>`_
2. Open a terminal window & run the downloaded script.

**Install Dependencies**

1. The IPOPT solver depends on the GNU FORTRAN, GOMP, Blas, and Lapack libraries.
   If these libraries are not already installed on your Linux system, you or your
   system administrator can use the sample commands below to install them. If you
   have a Linux distribution that is not listed, IPOPT should still work, but
   the commands to install the required libraries may differ. If these libraries
   are already installed, you can skip this and proceed with the next step.

   .. note:: Depending on your distribution, you may need to prepend ``sudo`` to
            these commands or switch to the "root" user.

   .. container:: collapsible

      .. container:: header

         Ubuntu 18.04 and 19.10 and distributions based on them

      .. code-block:: console

         sudo apt-get install libgfortran4 libgomp1 liblapack3 libblas3

   .. container:: collapsible

      .. container:: header

         Ubuntu 20.04 and distributions based on it

      .. code-block:: console

         sudo apt-get install libgfortran5 libgomp1 liblapack3 libblas3

   .. container:: collapsible

      .. container:: header

         Current RedHat based distributions, including CentOS

      .. code-block:: console

         yum install lapack blas libgfortran libgomp

Install IDAES-PSE
^^^^^^^^^^^^^^^^^^

.. include:: install_templates/conda_idaes_pse.txt

Install IDAES Extensions
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/extensions.txt

Install IDAES Examples
^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/examples.txt

Run IDAES Tests
^^^^^^^^^^^^^^^

.. include:: install_templates/run_tests.txt

