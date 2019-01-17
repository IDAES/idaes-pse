<!-- BEGIN Status badges -->
[![CircleCI](https://circleci.com/gh/dangunter/idaes-tmp.svg?style=svg&circle-token=1c18c546ab1e9e63002a6c37593f1f9efa03d587)](https://circleci.com/gh/dangunter/idaes-tmp)
[![codecov](https://codecov.io/gh/IDAES/idaes-tmp/branch/master/graph/badge.svg?token=5F5EzQws6o)](https://codecov.io/gh/IDAES/idaes-tmp)
<!-- END Status badges -->

# IDAES Toolkit

The IDAES Toolkit aims to provide multi-scale, simulation-based, open source
computational tools and models to support the design, analysis, optimization,
scale-up, operation and troubleshooting of innovative, advanced energy systems.

## Contributing

**By contributing to this repository, you are agreeing to all the terms set out
in the LICENSE.txt and COPYRIGHT.txt files in this directory.**

## Getting Started
For installation instructions, please refer to the documentation.

A pre-built PDF version of the IDAES documentation is in the file
`docs/IDAES.pdf`.

The documentation for IDAES is built using [Sphinx](http://www.sphinx-doc.org/). To generate the HTML version of the documentation, first make sure Sphinx is installed for your version of Python, then go to the "docs/" subdirectory and run the _Makefile_:

```
cd docs
make html
```

To view the documentation you just built, open the file
`doc/build/html/index.html` in a web browser.


## Running tests

After you install, you can run tests to make sure everything is working. We use [pytest](https://pytest.org/) for testing and generating code coverage reports.  The `pytest` command should be available in the conda environment created by running the `install.sh` script as described in the installation instructions.

To run tests against the core modules and DMF (not `idaes/contrib`), and generate a coverage report:

```
$ source activate <idaes_conda_env>  # If you used "install.sh <idaes_conda_env>" to install
$ pytest
```

To run tests in `idaes/contrib` just add that to the pytest command:

```
$ pytest ideas/contrib  # These are not guarenteed to all succeed...
```

If there are errors, or you are having trouble, you can use our [issue tracker on Github](https://github.com/IDAES/idaes/issues) to look for other users experiencing similar problems, or to report a new bug.


## Running a notebook

There are example [Jupyter](https://jupyter.org) notebook(s) in the `examples/` directory. To run them, you should invoke Jupyter on a Notebook file (these end in the extension `.ipynb`).

```
jupyter notebook examples/run-mea-model.ipynb
```

This should start up a notebook server and then pop up a tab or window in your default web browser showing the Notebook. For more information on how to use Jupyter, see the "Help" menu in the Notebook window itself, and the extensive documentation on the [Jupyter website](https://jupyter.org).

## IDAES on Docker Containers:

### JupyterHub instance on Amazon EC2 cluster

The JupyterHub instance is currently available for demo purposes and is started only when needed. It will be made available to IDAES users in the near future.

### Single-user image for development and testing

The Dockerfile in the top-level can be used to build a docker image that includes the IDAES package and its dependencies. The latest image is also maintained and **can be used for development and testing purposes**. 

#### Using the latest Docker image from DockerHub:

In our Jupyterhub deployment, this image serves as the single-user image that we use to spin up new containers for users to run Jupyter notebooks on. To pull the latest version of this image for development or testing, follow the steps outlined below. 

#### Docker installation: 

1. Install the community edition (CE) of [docker](https://docs.docker.com/install/). 

2. Start the docker daemon. Depending on your system, this will vary and you need to follow through with the install instructions (linked in step 1) for your specific operating system until you reach the step that starts the docker daemon. Here are some options for common operating systems:
   
      a. **OS X** : Docker should have been installed to your Applications directory. Browse to it and click on it from there. 
         You will see a small icon in your toolbar that indicates if the daemon's running successfully.
   
      b. **Ubuntu/Debian** : The Docker daemon will start automatically once Docker is installed.
   
      c. **CentOS** : Run `sudo systemctl start docker`.

Based on whether or not you have this repository cloned on your host machine, you should follow one or the other of the set of steps outlined under "Starting a new container with the repo cloned" or "Starting a new container with only the Docker image available."

##### Starting a new container with the repo cloned:

Use the script `idaes-docker` in the top-level directory of the repo as follows:

1. As a sanity check, run the unit tests on the image (which will have the latest IDAES master baked into it):

   ```
   hamdys-mbp:idaes helgammal$ ./idaes-docker test
   Running tests in container...
   =========================================================================================== test session starts ============================================================================================
   platform linux -- Python 3.6.7, pytest-4.0.2, py-1.7.0, pluggy-0.8.0
   rootdir: /home/idaes, inifile: pytest.ini
   plugins: cov-2.5.0
   collected 648 items  
   ...
   ================================================================================= 633 passed, 15 skipped in 89.09 seconds ==================================================================================

   ```

2. Run a Jupyter notebook from inside the container:

   ```
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
   ```
   Browse to the URL provided in the output message (in the example above this is `http://127.0.0.1:8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f`) and then start a new notebook from New -> Python 3 or browse to the IDAES example notebook under idaes/examples/heat_exchange_simple/simple_hx_flowsheet_01.ipynb. To shutdown the notebook server click "{Ctrl,Command} + c" in your terminal.

3. Refresh your IDAES docker image to the latest version from DockerHub:

   ```
   hamdys-mbp:idaes helgammal$ ./idaes-docker refresh
   Refreshing IDAES image from DockerHub...
   latest: Pulling from idaes/idaes_jupyterhub
   Digest: sha256:17e2c1d5d184cde71cd67477cac467af7d2da798e9f9a0a297f5c2f94bdeb1ac
   Status: Image is up to date for idaes/idaes_jupyterhub:latest
   ```

##### Starting a new container with only the Docker image available: 

1. Run the following command which will pull the latest IDAES image from DockerHub:

   ```
   docker pull idaes/idaes_jupyterhub:latest
   ```

2. Run the tests directly on the docker container by using the following command. If everything went well, all tests should pass.  

   ```
   docker run -it idaes/idaes_jupyterhub /bin/bash -c "cd /home/idaes && pytest"
   ```

3. There are then 2 ways to use the image: 

   a. Start a docker container and interact with it directly:
 
      ```
	  $ docker run -it idaes/idaes_jupyterhub /bin/bash
	  jovyan@10c11ca29008:~$ ls /home/
	  idaes  jovyan
	  jovyan@10c11ca29008:~$ cd idaes/
	  jovyan@10c11ca29008:~/idaes$ pytest
	  ...

      ```	  
  
   b. Start a docker container and use it to run Jupyter notebooks:

      ```
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
      ```
	 Browse to the URL provided in the output message (in the example above this is `http://127.0.0.1:8888/?token=348184135dacb8e7bd80f1bdcff5b34fff9012a9d79ecd0f`) and then start a new notebook from New -> Python 3 or browse to the IDAES example notebook under idaes/examples/heat_exchange_simple/simple_hx_flowsheet_01.ipynb. To shutdown the notebook server click "{Ctrl,Command} + c" in your terminal.

#### Build new image from Dockerfile:

- Run the build command. This will take some time to execute: 

  ```
  docker build .
  ```

- Tag the image. You can get IMAGE_NAME from the very last line in the previous step's output, for e.g: `Successfully built 88528d8e1f11` indicates the image name is `88528d8e1f11`.

  ```
  docker tag IMAGE_NAME idaes/idaes_jupyterhub:version_info_here
  ```
- You can then run a container as described in steps 3 and after in the previous section. 

## Contacts and more information

Please see the
[IDAES main website](https://www.idaes.org) for general information
and people to contact.
