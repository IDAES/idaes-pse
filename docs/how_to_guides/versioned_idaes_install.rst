Installing Specific IDAES Versions
==================================

This guide walks users through installing a version of IDAES other than the current stable release. 

.. warning:: The IDAES binary extensions are not yet supported on Mac/OSX.
             
             By installing IDAES outside of a conda environment, you will not be able
             to install the open source ipopt solver. 

We recommend using Conda to manage your environment & modules.

Setting up Conda Environment
----------------------------
Creating a New Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^
We strongly recommend using a separate environment for your versioned installation of IDAES.

.. code-block:: bash
   
   # Create a new Conda environment (Here named my-versioned-idaes-env, but this is up to you)
   conda create --yes --name my-versioned-idaes-env python=3.10
   conda activate my-versioned-idaes-env

Updating an Existing Install
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
While we recommend creating a fresh environment for your versioned install, but if you've already installed IDAES using the `Getting Started <../tutorials/getting_started/index>`_ guide using Conda for the package installation and would like to use that same environment, you'll need to uninstall the existing version of IDAES before you can install a different version using pip.

.. code-block:: bash

   # Update an existing Conda environment & remove idaes-pse (Here named my-idaes-env, but this is up to you)
   conda remove -n my-idaes-env idaes-pse
   activate my-idaes-env


Installing IDAES
----------------
.. _updating_install:

* To get a previous `IDAES release <https://github.com/IDAES/idaes-pse/releases>`_, for example v1.13 ::

   pip install idaes-pse==1.13

* To get the latest version from the GitHub main branch ::

   pip install 'idaes-pse[prerelease] @ https://github.com/IDAES/idaes-pse/archive/main.zip'

* To get a specific fork or branch, for example myfork (of idaes-pse) and mybranch ::

   pip install 'idaes-pse[prerelease] @ https://github.com/myfork/idaes-pse/archive/mybranch.zip'

* **For IDAES Contributors**: follow the :ref:`advanced user installation<tutorials/advanced_install/index:Advanced User Installation>`.