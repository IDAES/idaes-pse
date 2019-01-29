Installation Instructions
=========================

.. contents:: Contents

Installation using Docker
-------------------------
The simplest way to install the IDAES PSE Framework is by using
the pre-built Docker_ container.

Install Docker on your system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#. Install the community edition (CE) of Docker_ (website: https://docker.io).
#. Start the Docker daemon. Depending on your system, this will vary and you need to follow through with the install instructions (linked in step 1) for your specific operating system until you reach the step that starts the docker daemon. Here are some options for common operating systems:

      OS X
         Docker should have been installed to your Applications directory. Browse to it and click on it from there.
         You will see a small icon in your toolbar that indicates if the daemon's running successfully.

      Linux
        - Ubuntu/Debian: The Docker daemon will start automatically once Docker is installed.
        - CentOS: Start Docker manually, e.g., run
          ``sudo systemctl start docker``.

    Windows
        - Docker will be started automatically

.. _Docker: https://docker.io/

Get the IDAES Docker image
^^^^^^^^^^^^^^^^^^^^^^^^^^
Next you need to get the pre-built Docker "image" containing the source
code and solvers for the IDAES PSE framework. This image is stored in
an online service called "Docker Hub". Full details on how to do
this are available in the `Docker Hub documentation`_

.. _Docker Hub documentation: https://docs.docker.com/docker-hub/

Starting a new container with only the Docker image available:
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#. Run the following command which will pull the latest IDAES image from DockerHub:

   .. code-block:: sh

     docker pull idaes/idaes_jupyterhub:latest

#. Run the tests directly on the docker container by using the following command. If everything went well, all tests should pass.

   .. code-block:: sh

     docker run -it idaes/idaes_jupyterhub /bin/bash -c "cd /home/idaes && pytest"

#. There are then two basic ways to use the image:

   #. Start a docker container and interact with it directly:

    .. code-block:: sh

      $ docker run -it idaes/idaes_jupyterhub /bin/bash
      jovyan@10c11ca29008:~$ ls /home/
      idaes  jovyan
      jovyan@10c11ca29008:~$ cd idaes/
      jovyan@10c11ca29008:~/idaes$ pytest
      ...

   #. Start a docker container and use it to run Jupyter notebooks:

    .. code-block:: sh

      $ docker run -p 8888:8888 -it idaes/idaes_jupyterhub
      Container must be run with group "root" to update passwd file
      Executing the command: jupyter notebook
      [I 07:54:20.117 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
      [I 07:54:20.414 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.6/site-packages/jupyterlab
      [I 07:54:20.414 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
      [I 07:54:20.424 NotebookApp] Serving notebooks from local directory: /home
      [I 07:54:20.424 NotebookApp] The Jupyter Notebook is running at:
      [I 07:54:20.424 NotebookApp] http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254
      [I 07:54:20.424 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
      [C 07:54:20.424 NotebookApp]

        Copy/paste this URL into your browser when you connect for the first time,
        to login with a token:
            http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254

   Browse to the URL provided in the output message (in the example above this is `http://127.0.0.1:8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f`) and then start a new notebook from New -> Python 3 or browse to the IDAES example notebook under idaes/examples/heat_exchange_simple/simple_hx_flowsheet_01.ipynb. To shutdown the notebook server click "{Ctrl,Command} + c" in your terminal.

The IDAES toolkit is written in Python. It should run under versions of Python 2.7 and 3.6, and above. The toolkit uses [Pyomo](https://www.pyomo.org), a Python-based optimization language. See the Pyomo website for details.

.. note:: Although Python can run on most operating systems, *we are currently only
    supporting installation of the IDAES PSE framework on Linux*. This is due largely
    to complications of installing third-party solvers, not inherent properties
    of the PSE framework itself, and we plan to support Windows and Mac OSX
    installation in the not-too-distant future.

Simplified Linux Installation Instructions
------------------------------------------

The following instructions assume that

    * You have sudo privilege on your system.
    * You have ``apt`` or may install it, or know how to adapt the instructions to a different package manager or how to install the packages directly. ``apt`` is default on Debian-based Linux distributions, including Ubuntu.
    * You are on a computer+network that is allowed (by your sysadmins) to access and download the various packages and tools, including the solvers from third-party sources.

Install `prerequisite system applications <#system-prerequisites>`_:

.. code-block:: sh

    sudo apt-get install gcc g++ make

Download and install `miniconda <https://conda.io/docs/user-guide/install/linux.html>`_ (follow the prompts in the installer; you may want to restart your Terminal window afterwards to ensure that new environment variables are set correctly):

.. code-block:: sh
    
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

Create and activate a conda environment (along with its own copy of ``pip``) for the new IDAES installation **(you will need to** ``conda activate idaes`` **when you open a fresh terminal window and wish to use IDAES)**:

.. code-block:: sh
    
    conda create -n idaes pip
    conda activate idaes

Obtain the source code for IDAES from GitHub:

.. code-block:: sh

    git clone https://github.com/IDAES/idaes.git

Install the python dependencies:

.. code-block:: sh

    cd idaes
    pip install -r requirements.txt

`Install the IDAES framework itself <#install-idaes>`_:

.. code-block:: sh

    python setup.py develop

Install the `main solver dependencies <#other-dependencies>`_:

.. code-block:: sh

    sudo apt-get update && sudo apt-get install -y libboost-dev
    wget https://ampl.com/netlib/ampl/solvers.tgz
    tar -xf solvers.tgz
    ( cd solvers && ./configure && make )
    ( export ASL_BUILD=`pwd`/solvers/sys.x86_64.Linux && cd idaes/property_models/iapws95 && make )
    wget https://ampl.com/dl/open/ipopt/ipopt-linux64.zip
    unzip ipopt-linux64.zip
    sudo cp ipopt /usr/local/bin/

At this point, you should be able to launch the Jupyter Notebook server and successfully `run examples <examples.html>`_ from the ``examples`` folder:

.. code-block:: sh

    jupyter notebook

Installation on Linux/Unix
--------------------------

System Prerequisites
^^^^^^^^^^^^^^^^^^^^

The following commonly-used programs must be installed:

 - make
 - gcc
 - g++

GCC and G++ are necessary if you wish to compile and use the solver libraries. The following command installs all three, and assumes you have ``apt`` installed, which is default on Debian-based systems.

.. code-block:: sh

    sudo apt-get install gcc g++ make

Additionally, for full functionality you may wish to consult the `Other Dependencies`_.


Install IDAES
^^^^^^^^^^^^^^

* The installation instructions assume a Python packaging system called `Conda <https://conda.io/docs/>`_ is available. Please first consult the `Conda documentation <https://conda.io/docs/user-guide/>`_ to install this on your system. You can use either Anaconda or Miniconda.

* Conda allows you to to create separate environments containing files, packages and their dependencies that will not interact with other environments.

**Create/switch to your preferred Python environment**

.. code-block:: sh

  conda create -n idaes python=3 pyqt pip
  conda activate idaes

You can replace idaes with any name you like.  PyQt is used for some IDAES
graphical user interface elements. ``pip`` is already installed with conda itself,
but a copy needs to exist within the environment in order to cleanly encapsulate
all of the requirements and IDAES itself.

**Install the master branch of IDAES from GitHub:**

.. code-block:: sh

  git clone https://github.com/IDAES/idaes.git
  cd idaes

**Install the requirements**

.. code-block:: sh

  pip install -r requirements.txt

**Install the IDAES Framework**

  To compile C functions for some property models, the location of the compiled ASL is required
  for the commands below a location of :code:`$HOME/local/src/solvers/sys.x86_64.Linux`;
  however, this location will depend on your system and where you put the files.

  The BOOST_HEADER environment variable can be set optionally if the the build
  fails due to not finding BOOST. This allows more flexibility for alternative
  locations.  Setting BOOST_HEADER is usually not needed.

  If make fails or you do not want to compile, you can skip to the last line, but
  some property packages may not work.

.. code-block:: sh

  export ASL_BUILD=$HOME/local/src/solvers/sys.x86_64.Linux
  make
  python setup.py develop

**OR**

.. code-block:: sh

  export ASL_BUILD=$HOME/local/src/solvers/sys.x86_64.Linux
  make
  python setup.py install

IDAES on Docker Containers:
---------------------------

JupyterHub instance on Amazon EC2 cluster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The JupyterHub instance is currently available for demo purposes and is started only when needed. It will be made available to IDAES users in the near future.

Single-user image for development and testing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Dockerfile in the top-level can be used to build a docker image that includes the IDAES package and its dependencies. The latest image is also maintained and **can be used for development and testing purposes**. 

Using the latest Docker image from DockerHub:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In our Jupyterhub deployment, this image serves as the single-user image that we use to spin up new containers for users to run Jupyter notebooks on. To pull the latest version of this image for development or testing, follow the steps outlined below. 


Starting a new container with the repo cloned:
""""""""""""""""""""""""""""""""""""""""""""""

Use the script `idaes-docker` in the top-level directory of the repo as follows:

#. As a sanity check, run the unit tests on the image (which will have the latest IDAES master baked into it):

   .. code-block:: sh

     hamdys-mbp:idaes helgammal$ ./idaes-docker test
     Running tests in container...
     =========================================================================================== test session starts ============================================================================================
     platform linux -- Python 3.6.7, pytest-4.0.2, py-1.7.0, pluggy-0.8.0
     rootdir: /home/idaes, inifile: pytest.ini
     plugins: cov-2.5.0
     collected 648 items  
     ...
     ================================================================================= 633 passed, 15 skipped in 89.09 seconds ==================================================================================


#. Run a Jupyter notebook from inside the container:

   .. code-block:: sh

     hamdys-mbp:idaes helgammal$ ./idaes-docker notebook
     Starting Jupyter...
     Container must be run with group "root" to update passwd file
     Executing the command: jupyter notebook
     [I 19:10:08.161 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
     [I 19:10:08.452 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.6/site-packages/jupyterlab
     [I 19:10:08.452 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
     [I 19:10:08.470 NotebookApp] Serving notebooks from local directory: /home
     [I 19:10:08.470 NotebookApp] The Jupyter Notebook is running at:
     [I 19:10:08.470 NotebookApp] http://(a9e555672b1c or 127.0.0.1):8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f
     [I 19:10:08.471 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
     [C 19:10:08.472 NotebookApp]
     
         Copy/paste this URL into your browser when you connect for the first time,
         to login with a token:
             http://(a9e555672b1c or 127.0.0.1):8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f


   Browse to the URL provided in the output message (in the example above this is `http://127.0.0.1:8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f`) and then start a new notebook from New -> Python 3 or browse to the IDAES example notebook under idaes/examples/heat_exchange_simple/simple_hx_flowsheet_01.ipynb. To shutdown the notebook server click "{Ctrl,Command} + c" in your terminal.

#. Refresh your IDAES docker image to the latest version from DockerHub:

   .. code-block:: sh

     hamdys-mbp:idaes helgammal$ ./idaes-docker refresh
     Refreshing IDAES image from DockerHub...
     latest: Pulling from idaes/idaes_jupyterhub
     Digest: sha256:17e2c1d5d184cde71cd67477cac467af7d2da798e9f9a0a297f5c2f94bdeb1ac
     Status: Image is up to date for idaes/idaes_jupyterhub:latest

Starting a new container with only the Docker image available: 
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

#. Run the following command which will pull the latest IDAES image from DockerHub:

   .. code-block:: sh

     docker pull idaes/idaes_jupyterhub:latest

#. Run the tests directly on the docker container by using the following command. If everything went well, all tests should pass.  

   .. code-block:: sh

     docker run -it idaes/idaes_jupyterhub /bin/bash -c "cd /home/idaes && pytest"

#. There are then 2 ways to use the image: 

   #. Start a docker container and interact with it directly:
 
    .. code-block:: sh

      $ docker run -it idaes/idaes_jupyterhub /bin/bash
      jovyan@10c11ca29008:~$ ls /home/
      idaes  jovyan
      jovyan@10c11ca29008:~$ cd idaes/
      jovyan@10c11ca29008:~/idaes$ pytest
      ...
  
   #. Start a docker container and use it to run Jupyter notebooks:

    .. code-block:: sh

      $ docker run -p 8888:8888 -it idaes/idaes_jupyterhub
      Container must be run with group "root" to update passwd file
      Executing the command: jupyter notebook
      [I 07:54:20.117 NotebookApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
      [I 07:54:20.414 NotebookApp] JupyterLab extension loaded from /opt/conda/lib/python3.6/site-packages/jupyterlab
      [I 07:54:20.414 NotebookApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
      [I 07:54:20.424 NotebookApp] Serving notebooks from local directory: /home
      [I 07:54:20.424 NotebookApp] The Jupyter Notebook is running at:
      [I 07:54:20.424 NotebookApp] http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254
      [I 07:54:20.424 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
      [C 07:54:20.424 NotebookApp]

        Copy/paste this URL into your browser when you connect for the first time,
        to login with a token:
            http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254

   Browse to the URL provided in the output message (in the example above this is `http://127.0.0.1:8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f`) and then start a new notebook from New -> Python 3 or browse to the IDAES example notebook under idaes/examples/heat_exchange_simple/simple_hx_flowsheet_01.ipynb. To shutdown the notebook server click "{Ctrl,Command} + c" in your terminal.

Build new image from Dockerfile:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Run the build command from the root of the IDAES repo. This will take some time to execute: 

  .. code-block:: sh

    docker build .

* Tag the image. You can get IMAGE_NAME from the very last line in the previous step's output, for e.g: `Successfully built 88528d8e1f11` indicates the image name is `88528d8e1f11`.

  .. code-block:: sh

    docker tag IMAGE_NAME idaes/idaes_jupyterhub:version_info_here

* You can then run a container as described in steps 3 and after in the previous section.

Other Dependencies
------------------

Solvers
^^^^^^^

Some of the model code depends on external solvers. All of the solvers are optional to some extent, however IPOPT is used extensively.

**CPLEX**

* `Getting CPLEX <https://www.ibm.com/developerworks/community/blogs/jfp/entry/CPLEX_Is_Free_For_Students?lang=en>`_
* `Setting up CPLEX Python <http://www.ibm.com/support/knowledgecenter/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html>`_

**Gurobi**

* `Gurobi license <https://user.gurobi.com/download/licenses/free-academic>`_
* `Gurobi solver <http://www.gurobi.com/downloads/gurobi-optimizer>`_
* `Gurobi Python setup <http://www.gurobi.com/documentation/6.5/quickstart_mac/the_gurobi_python_interfac.html>`_

**IPOPT**

* Installing `IPOPT <https://www.coin-or.org/Ipopt/documentation/node10.html>`_

Function Dependencies
^^^^^^^^^^^^^^^^^^^^^

In some cases, IDAES uses AMPL user-defined functions written in C for property
models.  Compiling these functions is optional, but some models may not work
without them.

**ASL**

The AMPL solver library (ASL) is required, and can be downloaded from
from https://ampl.com/netlib/ampl/solvers.tgz.  Documentation is available at
https://ampl.com/resources/hooking-your-solver-to-ampl/. Typically to build the
ASL the files can be extracted, then in the directory with the ASL file run the
commands below.

.. code-block:: sh

  ./configure
  make

**Boost**

The C++ Boost libraries should be available. One possibility is to use conda to
install boost, but the best option depends on your system.

Installation on Windows
-----------------------

.. note:: We are NOT supporting Windows at this time. Some developers on the team have had success with the following instructions, but we do not promise that they will work for all users, nor will we prioritize helping debug problems.

Python Distribution
^^^^^^^^^^^^^^^^^^^

* Install `Anaconda for Windows <https://www.anaconda.com/download/#download>`_

* Add Anaconda and Anaconda scripts to the path "c:\users\<user>\Anaconda2\" and "c:\users\<user>\Anaconda2\Scripts\". To do this, search for "Edit system variables" in Windows search.  Click on "Edit system environment variables". Click on "Environment Variables". Under "System   Variables", search for the variable "Path" and click "Edit"

	.. image:: _static/install_windows_system_properties.png
	   :align: center
	   :scale: 75%



	1. For Windows 10:

	      1. In the new dialog box, click on "New" and add the path where you find the python.exe file. If you installed Anaconda2, this should be in “c:\users\<user>\Anaconda2\”. Copy the address and paste it here.


	      2. Repeat for "c:\users\<user>\Anaconda2\Scripts\".

   	2. For earlier versions:

	      1. Add path to the existing list, use semicolon as separator

	      2. Type "c:\users\<user>\Anaconda2\;c:\users\<user>\Anaconda2\Scripts\"

* Restart the command prompt and type `python`. If the path variable was added correctly, then you should be able to see the python interpreter as shown below.

.. image:: _static/install_windows_cmd_python.png
   :align: center
   :scale: 75%

Pyomo
^^^^^
* See `instructions <http://www.pyomo.org/installation/>`_ for pyomo installation. As mentioned, you can either use the pip or the conda install methods which come included with the Anaconda distribution but conda may be preferable if you installed Anaconda.

* To install pyomo using python’s **pip** package, follow these steps:


    1. Launch the "Anaconda prompt". You can find this in the start menu under Anaconda.

    2. Navigate to the "Scripts" folder in Anaconda. Or simply type, `where pip` in the prompt. This should return 1 paths and this should be in the scripts folder.

    3. Pip install pyomo from trunk (we recommend installing the IDAES branch of pyomo)

        1. Install the master branch of PyUtilib from GitHub using pip:

           `pip.exe install git+https://github.com/PyUtilib/pyutilib`

        2. Install the master branch of Pyomo from GitHub using pip:

           `pip.exe install git+https://github.com/Pyomo/pyomo@IDAES`

* To install using python’s **conda** package, follow the following steps:


    1. Launch the "Anaconda prompt". You can find this in the start menu under Anaconda.

    2. Navigate to the "Scripts" folder in Anaconda. Or simply type, `where conda` in the prompt. This should return 2 paths and one of these should be in the scripts folder.

    3. In the scripts folder run the following commands:

        `conda.exe install -c conda-forge pyomo`

        `conda.exe install -c conda-forge pyomo.extras`
* If the installation was successful, you should see the pyomo executable listed in the Scripts folder. You can check this using the `where pyomo` command.

IDAES
^^^^^

Option 1: Download zip file
"""""""""""""""""""""""""""
* From the `IDAES <https://github.com/IDAES/idaes>`_ repository on GitHub, click on "Clone or download" on the right in green. Click on “Download zip”.

* Extract the contents in the desired directory you want IDAES in.

* Open command prompt and navigate to the folder where you extracted the contents of the IDAES repository (`cd <user>/.../<desired directory>/IDAES/`).

    1. Run: `python setup.py develop`

Option 2: Using Git
"""""""""""""""""""

* Install `git <https://git-scm.com/download/win>`_ for Windows.

* If cloning the repository from the command line, move to a directory where you want to install the IDAES repository. Then run the following command:

	1. `git clone https://github.com/IDAES/idaes.git`

* Enter your github user id and password. The git installation in 1 should have added the git executable to your system path and you should be able to execute git commands from the command line.

* Open command prompt and navigate to the folder where you extracted the contents of the IDAES repository (`cd <user>/.../<desired directory>/IDAES/`).

   1. Run: `python setup.py develop`



