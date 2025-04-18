FROM bioconductor/bioconductor:RELEASE_3_20

LABEL maintainer="Sora Yonezawa <d246887@hiroshima-u.ac.jp>"
LABEL description="sorayone56/tssr-r-env: R environment for TSSr analysis"

# Install system dependencies
RUN apt-get update && apt-get -y upgrade && apt-get install -y --no-install-recommends \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    build-essential \
    libgmp-dev \
    libmagick++-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libgit2-dev \
    pandoc

# Install R packages
COPY install-r-packages_v2.sh /tmp/
RUN chmod +x /tmp/install-r-packages_v2.sh && \
    /tmp/install-r-packages_v2.sh && \
    # Remove the temporary script after execution
    rm /tmp/install-r-packages_v2.sh

ENV R_LIBS_USER=/home/rstudio/R/library
RUN mkdir -p /home/rstudio/R/library && chmod 777 -R /home/rstudio/R

# Set working directory
WORKDIR /work