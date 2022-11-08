FROM rocker/tidyverse:4.2

LABEL maintainer = "Ryan Corbett (corbettr@chop.edu)"

### Install apt-getable packages to start
#########################################
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Install required R packages from CRAN
RUN install2.r \
		ggpubr	\
		openxlsx \
		patchwork \
		survival \

ADD Dockerfile . 