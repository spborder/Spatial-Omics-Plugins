FROM python:3.11

LABEL maintainer="Sam Border - Computational Microscopy and Imaging Laboratory. <samuel.border@medicine.ufl.edu>"

RUN apt-get update && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    ca-certificates \
    libcurl4-openssl-dev \
    libexpat1-dev \
    graphviz \
    unzip \
    libhdf5-dev \
    software-properties-common \
    libssl-dev \
    # Standard build tools \
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    libmemcached-dev && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update ##[edited]

RUN apt-get install libxml2-dev libxslt1-dev -y

# Required for opencv-python (cv2)
RUN apt-get update && apt-get install ffmpeg libsm6 libxext6 -y

RUN apt-get update && apt-get install libpq-dev gdal-bin libgdal-dev -y
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

WORKDIR /
# Make Python3 the default and install pip.  Whichever is done last determines
# the default python version for pip.

#
RUN which  python && \
    python --version

ENV build_path=$PWD/build
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

ENV plugin_path=.
RUN mkdir -p ${plugin_path}

RUN apt-get update && \
    apt-get install -y --no-install-recommends memcached && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY . $plugin_path
WORKDIR $plugin_path

# Installing packages in setup.py
# Make sure to add links for large-image whls
RUN pip install --upgrade pip

RUN pip install --upgrade --ignore-installed pip setuptools && \
    pip install . --find-links https://girder.github.io/large_image_wheels && \
    rm -rf /root/.cache/pip/*


# Show what was installed
RUN python --version && pip --version && pip freeze

# Defining entrypoint
WORKDIR $plugin_path/cli
LABEL entry_path=$plugin_path/cli

# Testing entrypoint
RUN python -m slicer_cli_web.cli_list_entrypoint --list_cli
RUN python -m slicer_cli_web.cli_list_entrypoint CellSegmentation --help
RUN python -m slicer_cli_web.cli_list_entrypoint RegisterImages --help

ENV PYTHONUNBUFFERED=TRUE

ENTRYPOINT ["/bin/bash","docker-entrypoint.sh"]