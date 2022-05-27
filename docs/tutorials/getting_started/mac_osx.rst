Mac/OSX Installation Guide
==========================

Quickstart
----------

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
   :lines: 1-5,7,10-

------------------------------------------------

Installing IDAES
----------------

Install Prerequisites
^^^^^^^^^^^^^^^^^^^^^

**Install  Miniconda (optional)**

1. Download `Anaconda <https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh>`_
2. Open a terminal window & run the downloaded script.

Install IDAES-PSE
^^^^^^^^^^^^^^^^^^

.. include:: install_templates/conda_idaes_pse.md

.. container:: collapsible

   .. container:: header

      Installing other versions
   
   .. warning:: The IDAES binary extensions are not yet supported on Mac/OSX.

                By installing IDAES outside of a conda environment, you will be able
                to install the generic ipopt solver. 

   .. include:: install_templates/other_releases.md


Install IDAES Extensions
^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: install_templates/extensions.md

.. warning:: The IDAES binary extensions are not yet supported on Mac/OSX.

   .. raw:: html

      <details>
       <summary> Fallback solution </summary>

   As a fallback (assuming you are using a conda env) you can install
   the generic ipopt solver with the command ``conda install -c
   conda-forge ipopt`` though this will not have all the features
   of our extensions package.

   .. raw:: html

      </details>


Install IDAES Examples
^^^^^^^^^^^^^^^^^^^^^^

.. |os_specific_fpath| replace:: `~\/idaes/examples`
.. include:: install_templates/examples.md

Run IDAES Tests
^^^^^^^^^^^^^^^

.. include:: install_templates/run_tests.md