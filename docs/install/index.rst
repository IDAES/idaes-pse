Installation
============

.. toctree::
    :hidden:
    
    minimal
    linux-full
    windows-full
    docker

To install the IDAES PSE framework, follow the
set of instructions below that are appropriate for your needs and operating system.

If you get stuck, please contact `idaes-support@idaes.org <idaes-support@idaes.org>`_.

The minimal installation only installs IDAES and the free IPOPT_ solver with MUMPS_.
The full installation is recommended for access to more advanced solvers.
The Docker_ installation works on any platform that supports Docker, but
of course requires installation of, and some understanding of, Docker itself
to operate.

.. _Docker: https://www.docker.com/

.. _IPOPT: https://www.coin-or.org/Ipopt/documentation/documentation.html

.. _MUMPS: http://mumps.enseeiht.fr/

+-----------------------+------------------+-----------------------------+
| Type of installation  | Operating System | Section                     |
+=======================+==================+=============================+
| Minimal IPOPT/MUMPS   | Linux            | :ref:`min_install_linux`    |
+-----------------------+------------------+-----------------------------+
|                       | Windows          | :ref:`min_install_windows`  |
+-----------------------+------------------+-----------------------------+
|                       | Mac OSX          | not supported [1]_          |
+-----------------------+------------------+-----------------------------+
| Full                  | Linux            | :ref:`full_install_linux`   |
+-----------------------+------------------+-----------------------------+
|                       | Windows          | :ref:`full_install_windows` |
+-----------------------+------------------+-----------------------------+
|                       | Mac OSX          | not supported [1]_          |
+-----------------------+------------------+-----------------------------+
| Docker-based          | Windows, Linux   | :ref:`install_docker`       |
|                       | OSX              |                             |
+-----------------------+------------------+-----------------------------+

.. [1] For advanced users, Mac OSX installation may be performed with some small changes
       to the Linux installation instructions.

.. note::

    These installation procedures are only fully tested on Debian-based Linux
    distributions.



