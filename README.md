<!-- BEGIN Status badges -->
[![CircleCI](https://circleci.com/gh/dangunter/idaes-tmp.svg?style=svg&circle-token=1c18c546ab1e9e63002a6c37593f1f9efa03d587)](https://circleci.com/gh/dangunter/idaes-tmp)
[![codecov](https://codecov.io/gh/IDAES/idaes-tmp/branch/master/graph/badge.svg?token=5F5EzQws6o)](https://codecov.io/gh/IDAES/idaes-tmp)
<!-- END Status badges -->

# idaes-tmp

The core of the IDAES framework.

# Jupyterhub deployment 

## Amazon EC2 cluster

TODO

## Single-user image

The Dockerfile in the top-level can be used to build a docker image that includes the IDAES package and its dependencies. 

### Pull latest image:
In our Jupyterhub deployment, this image serves as the single-user image that we use to spin up new containers for users to run notebooks on. To pull the latest version of this image for development or testing, follow the steps below: 

1. Install [docker](https://docs.docker.com/install/). 

2. Start docker and run the following command: 

   ```
   docker pull idaes/idaes_jupyterhub:latest
   ```

3. Run the tests directly on the docker container by using the following command. If everything went well, all tests should pass.  

   ```
   docker run -it idaes/idaes_jupyterhub /bin/bash -c "cd /home/jovyan/idaes && pytest"
   ```

4. There are then 2 ways to use the image: 

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
	  [I 07:54:20.424 NotebookApp] Serving notebooks from local directory: /home/jovyan
	  [I 07:54:20.424 NotebookApp] The Jupyter Notebook is running at:
	  [I 07:54:20.424 NotebookApp] http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254
	  [I 07:54:20.424 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
	  [C 07:54:20.424 NotebookApp]

	    Copy/paste this URL into your browser when you connect for the first time,
	    to login with a token:
	        http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254
      ```
	 Copy the URL to your browser and start a Jupyter notebook from "New -> Python3". You should also be able to browse to (and run) an example notebook from the ideas repo in the top-level of your Jupyter setup. Go to idaes -> examples -> heat_exchange_simple -> simple_hx_flowsheet_01.ipynb

### Build new image from Dockerfile:

- Run the build command. This will take some time to execute: 

  ```
  docker build .
  ```

- Tag the image. You can get IMAGE_NAME from the very last line in the previous step's output, for e.g: `Successfully built 88528d8e1f11` indicates the image name is `88528d8e1f11`.

  ```
  docker tag IMAGE_NAME idaes/idaes_jupyterhub:version_info_here
  ```
- You can then run the image as described in steps 3 and after in the previous section. 