# syntax=docker/dockerfile:1.4
# Wiggle: A comprehensive sandbox environment for document processing, data science, and web automation
# The section without Taiga Special is from original Wiggle Dockerfile.
# Monorepo: api/api/sandbox_gateway/images/wiggle/default/Dockerfile
FROM --platform=linux/amd64 ubuntu:24.04

#############################################################
## Taiga Special
## Reminder:
## The container WILL NOT have internet access at runtime, 
## so anything the environment needs at runtime should be
## downloaded in the Dockerfile
#############################################################
ARG USERNAME="model"

# Remove default ubuntu user (GID/UID 1000)
RUN touch /var/mail/ubuntu && chown ubuntu /var/mail/ubuntu && userdel -r ubuntu


########################################################################################
## Taiga Special
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## !!!! Setting up user for the model to use first !!!!!
## This should happen before any packages are installed, to ensure user ID consistency
## The uid needs to be 1000, and the group needs to be 1000
## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##########################################################
RUN groupadd -g 1000 ${USERNAME} && \
    useradd -m -u 1000 -g 1000 ${USERNAME} && \
    gpasswd -d ${USERNAME} root || true

# Environment configuration
ENV IS_SANDBOX="yes" \
    PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive \
    JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64

# Layer 1: Core dependencies, runtimes, and fast-installing packages
# Pre-configure shared-mime-info to prevent hanging during installation
# IMPORTANT: DO NOT install AGPL-licensed packages via apt in this Docker image.
# AGPL packages require distributing the entire application source code under AGPL,
# which is incompatible with our licensing requirements.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    echo "shared-mime-info shared-mime-info/trigger-cache boolean false" | debconf-set-selections && \
    DEBIAN_FRONTEND=noninteractive apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    # Base utilities
    curl \
    wget \
    gnupg2 \
    software-properties-common \
    ca-certificates \
    git \
    vim \
    nano \
    sudo \
    build-essential \
    zip \
    unzip \
    file \
    bc \
    netcat-openbsd \
    libxml2-utils \
    # NSS and certificate bridge tools
    libnss3-tools \
    p11-kit \
    p11-kit-modules \
    # Python and Node.js
    #############################################################
    ## Taiga Special
    #############################################################
    #python3.11 \
    #python3.11-venv \
    #python3.11-dev \
    #python3-pip \
    #############################################################
    ## Taiga Special
    #############################################################
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    python3-wheel \
    python3-uno \
    pipx \
    nodejs \
    npm \
    # Python library dependencies
    libcairo2-dev \
    pkg-config \
    libgl1 \
    libglib2.0-0 \
    # Java runtime
    default-jre-headless \
    # Fonts
    && rm -rf /var/lib/apt/lists/*

# Create python symlink for compatibility
RUN ln -s /usr/bin/python3 /usr/bin/python



#############################################################
## Taiga Special: Python setup, version management, and uv
#############################################################
# Skip pip upgrade to avoid conflicts with system pip
# Use system Python explicitly to avoid conda Python
RUN /usr/bin/python3 -m pip install --break-system-packages uv
#############################################################
## Taiga Special: Python setup, version management, and uv
#############################################################

# Install core data science dependencies first using system Python
RUN --mount=type=cache,target=/root/.cache/uv \
    /usr/bin/python3 -m uv pip install --system --break-system-packages \
    numpy \
    pandas \
    matplotlib \
    scipy \
    scikit-learn \
    scikit-image


# Switch to root to install packages
USER root
# Install system dependencies first
RUN apt-get update && apt-get install -y wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
# Install mambaforge with Python 3.11+ base environment
# Using Miniforge3 (newer) instead of Mambaforge to get Python 3.11+
RUN wget https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Miniforge3-24.3.0-0-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh && \
    /opt/conda/bin/mamba update -n base python -y
# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"
# Configure conda channels (mambaforge already has conda-forge)
RUN conda config --add channels bioconda && \
    conda config --set channel_priority strict
# Create separate environments for each tool to avoid version conflicts
# Create all environments in a single RUN command to better handle dependencies
RUN mamba create -n blast -y blast && \
    mamba create -n samtools -y samtools && \
    mamba create -n bwa -y bwa && \
    mamba create -n bowtie2 -y bowtie2 && \
    mamba create -n bedtools -y bedtools && \
    mamba create -n fastqc -y fastqc && \
    mamba create -n seqtk -y seqtk && \
    mamba create -n star -y -c bioconda star && \
    mamba create -n subread -y -c bioconda subread && \
    mamba create -n multiqc -y -c bioconda multiqc && \
    mamba create -n bcftools    -y -c bioconda bcftools && \
    mamba create -n gatk        -y -c bioconda gatk4 && \
    mamba create -n eggnogmapper -y -c bioconda eggnog-mapper && \
    mamba create -n orthofinder  -y -c bioconda orthofinder && \
    mamba create -n ivar         -y -c bioconda ivar && \
    mamba create -n freebayes    -y -c bioconda freebayes && \
    mamba clean --all -y
RUN  mamba create -n nextflow -y -c bioconda nextflow 
# Create R environment with base R
# Install R base first, then add packages incrementally to avoid solver issues
RUN mamba create -n R -y -c conda-forge r-base=4.3.* && \
    mamba clean --all -y
# Install R packages in smaller groups to avoid dependency conflicts
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             source /opt/conda/etc/profile.d/mamba.sh && \
             conda activate R && \
             mamba install -y -c conda-forge r-tidyverse r-devtools r-biocmanager zlib" && \
    mamba clean --all -y
# Create Python environment with common bioinformatics and data science packages
RUN mamba create -n python -y -c conda-forge python=3.11 && \
    mamba clean --all -y
# Install canonical Python packages for bioinformatics and data science
RUN bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             source /opt/conda/etc/profile.d/mamba.sh && \
             conda activate python && \
             mamba install -y -c conda-forge \
                 numpy \
                 pandas \
                 scipy \
                 matplotlib \
                 seaborn \
                 scikit-learn \
                 jupyterlab \
                 ipython \
                 biopython \
                 pysam \
                 pyvcf3 \
                 requests \
                 plotly \
                 statsmodels \
                 networkx \
                 h5py" && \
    mamba clean --all -y
# Optional: Install Bioconductor packages (comment out if build time is too long)
RUN echo '#!/usr/bin/env Rscript\n\
if (!require("BiocManager", quietly = TRUE))\n\
     install.packages("BiocManager", repos="https://cloud.r-project.org")\n\
 BiocManager::install(version = "3.20", ask = FALSE, update = TRUE)\n\
# # Install commonly used Bioconductor packages\n\
 BiocManager::install(c(\n\
     "DESeq2",\n\
     "edgeR",\n\
     "limma",\n\
     "GenomicRanges",\n\
     "Biostrings"\n\
 ), ask = FALSE, update = FALSE)' > /tmp/install_bioc.R && \
     bash -c "source /opt/conda/etc/profile.d/conda.sh && \
              source /opt/conda/etc/profile.d/mamba.sh && \
              conda activate R && \
              Rscript /tmp/install_bioc.R" && \
     rm /tmp/install_bioc.R
# Create workspace directory if it doesn't exist
RUN mkdir -p /workspace
# No need to change conda ownership - it's already world-readable and executable
# The model user can use conda commands without owning the files
# Initialize conda for bash shell (mamba uses the same init as conda)
RUN /opt/conda/bin/conda init bash && \
    echo "source /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc && \
    echo "source /opt/conda/etc/profile.d/mamba.sh" >> /etc/bash.bashrc
# Switch back to root for installing packages
USER root
# Ensure PATH includes conda for the user
ENV PATH="/opt/conda/bin:$PATH"
# Initialize conda/mamba for the model user's bashrc
RUN touch /home/${USERNAME}/.bashrc && \
    echo "source /opt/conda/etc/profile.d/conda.sh" >> /home/${USERNAME}/.bashrc && \
    echo "source /opt/conda/etc/profile.d/mamba.sh" >> /home/${USERNAME}/.bashrc && \
    chown ${USERNAME}:${USERNAME} /home/${USERNAME}/.bashrc



##########################################################################
## Taiga Special
## Creating a directory with root only permissions
## for MCP and tools. 
## !!!!! The `model` user should not have access to this directory !!!!
##########################################################################
WORKDIR /mcp_server
RUN chmod -R 0700 /mcp_server

##########################################################################
## Taiga Special
## Copying over tools, pyproject, etc, into this
## restricted access directory
##########################################################################
# Copy taiga core package from named 'taiga' context (../../)
# COPY --from=taiga taiga-core ./taiga-core

# Copy the rubric example from the taiga context
# COPY --from=taiga examples/rubric ./examples/rubric

# Copy the bioinformatics example from default context (current directory when building)
# COPY . ./examples/bioinformatics

# Install the taiga package and rubric example using system Python explicitly
# Use the system Python 3 explicitly, not conda's Python
#RUN /usr/bin/python3 -m uv pip install --system --break-system-packages -e ./taiga-core
#RUN /usr/bin/python3 -m uv pip install --system --break-system-packages -e ./examples/rubric
#RUN /usr/bin/python3 -m uv run rubric hello

##########################################################
## Taiga Special
## Configuring a workdir for the `model` user
##########################################################
ENV WORKDIR=/workdir

RUN mkdir -p ${WORKDIR}

RUN chown -R ${USERNAME}:${USERNAME} ${WORKDIR}

WORKDIR ${WORKDIR}

# Switch to the model user (UID 1000) for runtime
USER ${USERNAME}