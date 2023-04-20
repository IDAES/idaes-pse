Advanced User Installation
==========================

Advanced users who plan to develop their own models or tools are encouraged to install IDAES using Git and GitHub as described in this section, rather than using the instructions in the :ref:`getting started<tutorials/getting_started/index:Getting Started>` section. These advanced users will greatly benefit from improved version control and code integration capabilities.

.. contents:: :local:

Git and GitHub Basics
---------------------

Git is a distributed version control system that keeps track of changes in a set of files, while GitHub is a hosting service for Git repositories that adds many other features that are useful for collaborative software development.

Both Git and GitHub are widely used and there are excellent tutorials and resources for each. See `Atlassian Github tutorials <https://www.atlassian.com/git/tutorials>`_ , `GitHub help <https://help.github.com/>`_, and `Git documentation <https://git-scm.com/doc>`_.

A limited reference for Git and GitHub terminology and commands is provided :ref:`here<tutorials/advanced_install/terms_commands:Terminology and Commands>`, users that are new to Git and GitHub are strongly encouraged to use the more detailed resources above.

.. toctree::
    :hidden:

    terms_commands

Installation with GitHub
------------------------

The main IDAES GitHub repository is idaes-pse. This repository includes the core framework, model libraries, and integrated tools. It contains all of the release versions of IDAES and is frequently updated with new features.

The following instructions describe how to install and update the idaes-pse repository.


Github Setup
^^^^^^^^^^^^^

In order to use GitHub, you need to create a login on `GitHub <https://github.com>`_.

Fork the Repository
^^^^^^^^^^^^^^^^^^^
You use a "fork" of a repository (or "repo" for short) to create a space where you have complete control and can make changes without directly affecting the main repository.

.. figure:: /images/github-fork-repo_pse.png
    :align: right
    :width: 500px

    Figure 1. Screenshot showing where to click to fork the Github repo

You should first visit the idaes-pse repo on Github at https://github.com/IDAES/idaes-pse/. Then you should click on the fork icon in the top right and click on your username. These steps will have created your own fork of the repo with the same name under your username.

Clone Your Fork
^^^^^^^^^^^^^^^
A "clone" is a copy of a Github repository on your local machine. This is what you need to do in order to actually edit and change the files. To make a clone of the fork you created in the previous step, change to a directory where you want to put the source code and run the command::

    git clone https://github.com/MYNAME/idaes-pse.git
    cd idaes-pse

Of course, replace MYNAME with your username. This will download all the files in the latest version of the repository onto your local disk.

.. note:: After the ``git clone``, subsequent Git commands should be performed from the "idaes-pse" directory.

Add Upstream Remote
^^^^^^^^^^^^^^^^^^^
In order to guarantee that your fork can be synchronized with the "main" idaes-pse repo in the GitHub IDAES organization, you need to add a pointer to that repository as a *remote*. This repository will be called *upstream* and linked with the following command::

    git remote add upstream https://github.com/IDAES/idaes-pse.git

To check to see if you added the remote correctly use the following command::

     git remote -v

You should see that there are two remotes, origin and upstream. Both have two lines showing the remote name, the url, and the access (fetch or push). Origin is the pointer to your fork and was automatically added with the clone command, while upstream is the pointer to the main idaes-pse repo that you just added.

Create the Python Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Once you have the repo cloned, you can change into that directory (by default, it will be called "idaes-pse" like the repo) and install the Python packages.

But before you do that, you need to get the Python package manager fully up and running. We use a Python packaging system called Conda_ and we specifically use its minimal version Miniconda_. If you do not already have Conda, please follow the installation instructions for your operating system in :ref:`getting started<tutorials/getting_started/index:Getting Started>`.

.. _Conda: https://conda.io/
.. _Miniconda: https://conda.io/en/latest/miniconda.html

After Miniconda is installed, we recommend creating a separate conda environment for IDAES. If you are unfamiliar with environments, a good starting guide is `here <https://towardsdatascience.com/getting-started-with-python-environments-using-conda-32e9f2779307>`__. Create and activate a conda environment for the new IDAES installation with the following commands (we officially support python |python-min| through |python-max|, with |python-default| as a default):

.. code-block:: sh

    conda create -n idaes python=|python-default|
    conda activate idaes

.. note::
    When setting up a conda environment like this, you **must** ``conda activate idaes`` whenever you open a fresh terminal window and wish to use IDAES.

Finish the Installation
^^^^^^^^^^^^^^^^^^^^^^^

Now that conda and pip are installed, and you are in the "idaes" conda environment, you can run the following commands to install the requirements for IDAES and get the extensions (e.g. binaries for the solver IPOPT and other external function calls):

.. code-block:: sh

    pip install -r requirements-dev.txt
    idaes get-extensions

.. warning::
    The IDAES binary extensions are not yet supported on Mac/OSX

.. note::
    This ``pip install`` command would override any package within the conda environment,
    so if you would like to use a specific version of a package (e.g. a local clone of the Pyomo git repository), you should look at the
    ``requirements-dev.txt`` file and use it as a reference to either install the individual packages manually, or create a separate requirements file customized to your development use case.

You can test that everything is installed properly by running the tests with
Pytest_:

.. code-block:: sh

    pytest -m "not integration"

The not integration tag skips some tests that are slow. If you like, you can run all of the tests with just ``pytest``.

.. _Pytest: https://pytest.org/

Update IDAES
^^^^^^^^^^^^^^^^^^^^^^^
The main branch of idaes-pse is frequently updated and a new IDAES release occurs quarterly. It is recommended that you update your fork and local repositories and conda environment periodically.

.. code-block:: sh

    pip install -r requirements-dev.txt
    idaes get-extensions  # IDAES extensions should also be updated
