FROM satijalab/seurat:5.0.0

LABEL maintainer="Sam Border CMI Lab <samuel.border@medicine.ufl.edu>"

RUN apt-get update && \
    apt-get install --yes --no-install-recommends software-properties-common gpg-agent && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get --yes --no-install-recommends -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confold" dist-upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git \
    wget \
    curl \
    ca-certificates \
    libcurl4-openssl-dev \
    libexpat1-dev \
    unzip \
    libhdf5-dev \
    libpython3-dev \
    python3.10-dev \
    python3.10-distutils \
    libssl-dev \
    libffi-dev \
    # Standard build tools \
    build-essential \
    cmake \
    autoconf \
    automake \
    libtool \
    pkg-config \
    # useful later \
    libmemcached-dev && \
    #apt-get autoremove && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


RUN apt-get update ##[edited]

WORKDIR /
# Make Python3 the default and install pip.  Whichever is done last determines
# the default python version for pip.

# Make a specific version of python the default and install pip
RUN rm -f /usr/bin/python && \
    rm -f /usr/bin/python3 && \
    ln `which python3.10` /usr/bin/python && \
    ln `which python3.10` /usr/bin/python3 && \
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python get-pip.py && \
    rm get-pip.py
    #ln `which pip3` /usr/bin/pip 

RUN which  python && \
    python --version

# Required for R package installations
RUN apt-get update && \
    apt-get install -y libv8-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk-dev \
    libgsl-dev \
    libpcre2-dev \
    libudunits2-dev \
    libgdal-dev \
    libpq-dev \
    unixodbc \
    unixodbc-dev \
    libfontconfig1-dev \
    libcairo2-dev \ 
    libharfbuzz-dev \
    libfribidi-dev \
    libmagick++-dev

## Taken from Azimuth Dockerfile
RUN mkdir lzf
WORKDIR /lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf_filter.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf_filter.h
RUN mkdir lzf
WORKDIR /lzf/lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf_c.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf_d.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzfP.h https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf.h
WORKDIR /lzf
RUN gcc -O2 -fPIC -shared lzf/*.c lzf_filter.c -I /usr/include/hdf5/serial/ -lhdf5_serial -o liblzf_filter.so
WORKDIR /
ENV HDF5_PLUGIN_PATH=/lzf

## Installing R packages
COPY install_R_packages.r .
RUN R -e 'remotes::install_version("Matrix",version="1.6.4",repos="https://cran.r-project.org",dependencies=TRUE)'
RUN R -e 'install.packages("SeuratObject",version=">= 5.0.2",repos="https://cran.r-project.org",dependencies=TRUE)'
RUN Rscript install_R_packages.r BiocManager BSgenome.Hsapiens.UCSC.hg38 glmGamPoi GenomeInfoDb GenomicRanges TFBSTools JASPAR2020 EnsDb.Hsapiens.v86 IRanges Rsamtools S4Vectors
RUN R -e 'remotes::install_github("satijalab/azimuth",ref="master",dependencies=TRUE)'
RUN R -e 'remotes::install_github("ctlab/fgsea",dependencies=TRUE)'
RUN R -e 'remotes::install_github("JEFworks-Lab/STdeconvolve",dependencies=TRUE)'

ENV build_path=/build
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Copying over plugin files
ENV plugin_path=.
RUN mkdir -p $plugin_path

RUN apt-get update && \
    apt-get install -y --no-install-recommends memcached && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY . $plugin_path
WORKDIR $plugin_path

RUN pip install --no-cache-dir --upgrade --ignore-installed pip setuptools && \
    pip install --no-cache-dir .  && \
    rm -rf /root/.cache/pip/*

# Show what was installed
RUN python --version && pip --version && pip freeze

# Defining entrypoint
WORKDIR $plugin_path/cli
LABEL entry_path=$plugin_path/cli

# Testing entrypoint
RUN python -m slicer_cli_web.cli_list_entrypoint --list_cli
RUN python -m slicer_cli_web.cli_list_entrypoint CellDeconvolution --help
RUN python -m slicer_cli_web.cli_list_entrypoint SpotAnnotation --help

ENV PYTHONBUFFERED=TRUE

ENTRYPOINT ["/bin/bash","docker-entrypoint.sh"]