Docker Container
----------------
This page documents information needed by developers for working with the IDAES
docker container.

As is expected by Docker, the main file for creating the Docker
image is the "Dockerfile" in the top-level directory.

docker-idaes script
^^^^^^^^^^^^^^^^^^^
You can build new Docker images using the ``create`` option to the
`docker-idaes` script. For example::

    ./docker-idaes create

You need to have the IDAES installation activated. The script will automatically
find the current version and attempt to build a Docker image with the same version.
If it detects an existing image, it will skip the image build. Next, the script will
try to use ``docker save`` to save the image as a compressed archive. This will
also be skipped if an existing image file, with the same version as the "idaes"
Python package, is detected.

Pushing an image to S3
^^^^^^^^^^^^^^^^^^^^^^
The Docker images are stored on Amazon S3. Before you can upload a new image,
you need to be part of the "IDAES-admin" group that is part of Amazon's
IAM (Identity Access Management) system. Please contact one of the core
developers to learn how to join this IAM group.

Once you have the IAM keys, you need to create a file ``~/.aws/credentials``
that has the access key id and key from the IAM account. It will look like this::

    [default]
    aws_access_key_id = IDGOESHERE
    aws_secret_access_key = accesskeygoeshere

The values for the ID and Access key are available from the AWS "IAM"
service console.

Next you need to use the AWS command-line tools to copy the local image
up to Amazon S3. For example, if the image was version "1.0.1", you would
use the following command::

    aws s3 cp idaes-pse-docker-1.0.1.tgz \
           s3://idaes/idaes-pse/idaes-pse-docker-1.0.1.tgz

If the new image should be the latest, you also need to do an S3 -> S3 copy to
create a new latest image::

    aws s3 cp s3://idaes/idaes-pse/idaes-pse-docker-1.0.1.tgz \
              s3://idaes/idaes-pse/idaes-pse-docker-latest.tgz
