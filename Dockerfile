# Use Ubuntu 20.04 as base image
FROM ubuntu:20.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app

# Install system dependencies and Miniconda
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh

# Add conda to path
ENV PATH="/opt/conda/bin:${PATH}"

# Create conda environment with R and basic packages
RUN conda create -n ace_env -c conda-forge -c bioconda -c r \
    r-base=4.2 \
    r-essentials \
    bioconductor-genomicranges \
    bioconductor-rtracklayer \
    bioconductor-qdnaseq \
    bioconductor-biobase \
    r-optparse \
    r-data.table \
    r-ggplot2 \
    r-reshape2 \
    r-gridextra \
    && conda clean -afy

# Install BiocManager and ACE in the conda environment
RUN conda run -n ace_env R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cran.rstudio.com/"); BiocManager::install("ACE", dependencies=TRUE)'

# Verify installation
RUN conda run -n ace_env R -e 'library(ACE); packageVersion("ACE")'

# Create activation script
RUN echo '#!/bin/bash\n\
source /opt/conda/etc/profile.d/conda.sh\n\
conda activate ace_env\n\
exec "$@"' > /usr/local/bin/entrypoint.sh && \
chmod +x /usr/local/bin/entrypoint.sh

# Set entrypoint to activate conda environment
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["R"] 