FROM rocker/tidyverse:4.2

LABEL maintainer = "Ryan Corbett (corbettr@chop.edu)"

COPY scripts/install_github.r .

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
  data.table \
  ggpubr \
  ggthemes \
	optparse \
	pheatmap \
	RColorBrewer \
	survival \
  survMisc \
  survminer \
  tidytext
  
RUN ./install_github.r \
	clauswilke/colorblindr
	
WORKDIR /rocker-build/

ADD Dockerfile . 