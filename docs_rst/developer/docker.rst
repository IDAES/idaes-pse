Docker container
----------------
This page documents information needed by developers for working with the IDAES
docker container.

As is expected by Docker, the main file for creating the Docker
image is the "Dockerfile" in the top-level directory.

docker-idaes script
^^^^^^^^^^^^^^^^^^^
You can build new Docker images using the ``create`` option to the docker-idaes script.
For example:

    ./docker-idaes create

You need to have the IDAES installation activated. The script will automatically
find the current version and attempt to build a Docker image with the same version.
If it detects an existing image, it will skip the image build. Next, the script will
try to use ``docker save`` to save the image as a compressed archive. This will
also be skipped if an existing image file, with the same version as the "idaes"
Python package, is detected.