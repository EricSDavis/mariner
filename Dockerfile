# Use Bioconductor base image
FROM bioconductor/bioconductor_docker:devel

# Install mariner
RUN R -e 'remotes::install_github("EricSDavis/mariner@dev", dependencies=TRUE)'

WORKDIR .
CMD ["/init"]
