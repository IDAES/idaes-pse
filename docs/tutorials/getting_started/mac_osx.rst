Mac/OSX Installation Guide
==========================

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
.. literalinclude:: install_templates/quickstart.md
   :language: bash
   :linenos:
   :start-after: modification warning
   :lines: 1-6,10,11-

------------------------------------------------

Installing IDAES
----------------
To get IDAES fully set up on your machine, we'll go through the steps to get idaes-pse package installed as well as setting up the IDAES extensions, which includes some extra solvers and function libraries, the IDAES example files, and the IDAES tests.

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install  Miniconda**

1. Download `Miniconda <https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh>`_
2. Open a terminal window & run the downloaded script.

Install IDAES-PSE
^^^^^^^^^^^^^^^^^^

.. include:: install_templates/conda_idaes_pse.md


Install IDAES Extensions
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/extensions.md

.. warning:: 

   The IDAES binary extensions are not yet supported on Mac/OSX.

   .. container:: collapsible

      .. container:: header

         **Fallback solution**

      As a fallback (assuming you are using a conda env) you can install
      the open source ipopt solver with the command ``conda install -c
      conda-forge ipopt`` though this will not have all the features
      of our extensions package.


Install IDAES Examples
^^^^^^^^^^^^^^^^^^^^^^

.. |os_specific_fpath| replace:: `~\/idaes/examples`
.. include:: install_templates/examples.md

Run IDAES Tests
^^^^^^^^^^^^^^^

.. include:: install_templates/run_tests.md
