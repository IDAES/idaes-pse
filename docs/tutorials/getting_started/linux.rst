Linux Installation Guide
========================

Prerequisites
-------------

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

.. |os_specific_fpath| replace:: `~\/idaes/examples`

.. include:: generic_install.md