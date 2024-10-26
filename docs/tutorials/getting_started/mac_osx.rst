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

.. literalinclude:: install_templates/quickstart.txt
   :language: bash

.. warning:: Currently, the macOS x86_64 binaries do not include HSL linear solvers, k_aug, or dot_sens, so
   some linear solvers in Ipopt and some uncertainty propagation features may not be available. Some tests
   may also fail due to missing features.

------------------------------------------------

Installing IDAES
----------------
To get IDAES fully set up on your machine, we'll go through the steps to get idaes-pse package installed as well as setting up the IDAES extensions, which includes some extra solvers and function libraries, the IDAES example files, and the IDAES tests.

.. include:: install_templates/idaes_data.txt

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install  Miniforge**
-- If your device features apple silicon (M1/2/3/X):
1. Download `Miniforge <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh>`_
2. Open a terminal window & run the downloaded script.

-- If your device has an Intel Chip (x86_64):
1. Download `Miniforge <https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh>`_
2. Open a terminal window & run the downloaded script.

To determine your chip type:
** Click the Apple icon on the top right of your screen --> `About This Mac` **

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