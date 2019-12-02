.. _install_docker:

Installation using Docker
=========================
One way to install the IDAES PSE Framework is by using
the pre-built Docker_ image.

.. _Docker: https://www.docker.com/

A Docker image is essentially an embedded
instance of Linux (even if you are using Windows or Mac OSX)
that has all the code for the IDAES PSE framework
pre-installed. You can run commands and Jupyter Notebooks in that
image. This section describes how to set up your system, get the
Docker image, and interact with it.

Install Docker on your system
-----------------------------
#. Install the community edition (CE) of Docker_ (website: https://docker.io).
#. Start the Docker daemon. How to do this will depend on your operating system.

      OS X
         You should install `Docker Desktop for Mac`_.
         Docker should have been installed to your Applications directory. Browse to it and click on it from there.
         You will see a small icon in your toolbar that indicates
         that the daemon is running.

      Linux
         Install Docker using the package manager for your OS. Then
         start the daemon. If you are using Ubuntu or a Debian-based Linux distro,
         the Docker daemon will start automatically once Docker is installed.
         For CentOS, start Docker manually, e.g., run ``sudo systemctl start docker``.

      Windows
        You should install `Docker Desktop for Windows`_.
        Docker will be started automatically.

.. _Docker Desktop for Mac: https://docs.docker.com/docker-for-mac/install/
.. _Docker Desktop for Windows: https://docs.docker.com/docker-for-windows/install/

Get the IDAES Docker image
--------------------------
You need to get the ready made Docker image containing the source
code and solvers for the IDAES PSE framework. This image is available
for download from DockerHub (an online portal where Docker images are stored). Images
on DockerHub are versioned according to the release version.
See the Releases_ page on GitHub
for information about what is different about each version.

If you want the latest version, simply use the tag "latest" as the version number.
Thus, **running the following in a terminal will download the latest version**:

.. code-block:: shell
    
    docker pull idaes/jupyterhub:latest

-.. _Releases: https://github.com/IDAES/idaes-pse/releases

Run the IDAES Docker image
--------------------------

To start the Docker image, use a graphical user interface or a console or shell
command-line interface.

From the command-line, if you want to start up the Jupyter Notebook server, e.g.
to view and run the examples and tutorials, then run this command:

.. code-block:: console

      $ docker run -p 8888:8888 -it idaes/jupyterhub
      ... <debugging output from Jupyter>
      ...
      Copy/paste this URL into your browser when you connect for the first time,
      to login with a token:
          http://(305491ce063a or 127.0.0.1):8888/?token=812a290619211bef9177b0e8c0fd7e4d1f673d29909ac254

Copy and paste the URL provided at the end of the output into a browser window
and you should get a working Jupyter Notebook. You can browse to the examples
directory under ``/home/idaes/examples`` and click on the Jupyter Notebooks to
open them.

To interact with the image directly from the command-line (console), you can run the
following command:

.. code-block:: console

      $ docker run -p 8888:8888 -it idaes/jupyterhub /bin/bash
      jovyan@10c11ca29008:~$ cd /home/idaes
      ...

