# Use an official R runtime as a parent image
FROM rocker/shiny:latest
#  "shiny"

# Install system dependencies including those for tidyverse
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git libxml2-dev libmagick++-dev \
    wget libgomp1 \
    libssl-dev \
    make \
    pandoc && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

#    libcurl4-gnutls-dev \
#    libcurl4-openssl-dev \
# Install Blast
# FROM ncbi/blast:2.13.0

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.13.0/ncbi-blast-2.13.0+-x64-linux.tar.gz -O /tmp/blast.tar.gz && \
	tar zxvf /tmp/blast.tar.gz -C /opt/ && rm /tmp/blast.tar.gz

ENV PATH="/opt/ncbi-blast-2.13.0+/bin:${PATH}"

# Install Basic Utility R Packages

RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('KEGGREST')"
#RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "remotes::install_version('tidyverse', version = '2.0.0', dependencies= T)"
RUN R -e "remotes::install_version('arrow', version = '14.0.0', dependencies= T)"
RUN R -e "remotes::install_version('data.table', version = '1.14.8', dependencies= T)"
RUN R -e "remotes::install_version('rBLAST', repos = 'https://mhahsler.r-universe.dev', dependencies= T)"
RUN R -e "remotes::install_version('leaflet', version = '2.1.2', dependencies= T)"
RUN R -e "remotes::install_version('sp', version = '1.6-1', dependencies= T)"
RUN R -e "remotes::install_version('magrittr', version = '2.0.3', dependencies= T)"
RUN R -e "remotes::install_version('shinythemes', version = '1.2.0', dependencies= T)"
RUN R -e "remotes::install_version('DT', version = '0.28', dependencies= T)"
RUN R -e "remotes::install_version('shinyjs', version = '2.1.0', dependencies= T)"
RUN R -e "remotes::install_version('vembedr', version = '0.1.5', dependencies= T)"
RUN R -e "remotes::install_version('multidplyr', version = '0.1.3', dependencies= T)"
RUN R -e "remotes::install_version('shinybusy', version = '0.3.1', dependencies= T)"

RUN rm -rf /srv/shiny-server/*
COPY /app/* /srv/shiny-server/

USER shiny

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
