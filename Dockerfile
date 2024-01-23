FROM rocker/tidyverse:4.2

LABEL maintainer = "Ryan Corbett (corbettr@chop.edu)"

COPY scripts/install_github.r .

COPY scripts/install_bioc.r .

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2 \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev

# libmagick++-dev is needed for coloblindr to install
RUN apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev

# install R packages from CRAN
RUN install2.r \
	BiocManager \
	circlize \
  data.table \
  ggalluvial \
  ggpubr \
  ggthemes \
	optparse \
	pheatmap \
	RColorBrewer \
	survival \
  survMisc \
  survminer \
  tidytext
  
  # install R packages from Bioconductor 
RUN ./install_bioc.r \
  biomaRt \
  ComplexHeatmap \
  GenomicRanges
  
RUN ./install_github.r \
	clauswilke/colorblindr
	
RUN ./install_github.r  'thomasp85/patchwork' --ref 'c67c6603ba59dd46899f17197f9858bc5672e9f4'
	
WORKDIR /rocker-build/

ADD Dockerfile . 