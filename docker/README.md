# Docker image for HILA

Create docker image:

    docker build -t hila -f Dockerfile .

Launch image interactively with docker compose

    docker compose run --rm hila-applications

## Developing with docker

The applications folder is automatically mounted from the local host to the docker image when launching the service hila-applications

    ../applications:/HILA/applications

This allows one to develop HILA applications directly from source and launch them in the docker image with ease.

When developing hila libraries and hilapp one can also launch the service hila-source which mounts the HILA/libraries and HILA/hilapp/src folders to the container

    docker compose run --rm hila-source