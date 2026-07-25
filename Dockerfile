# Use Bioconductor base image
FROM bioconductor/bioconductor_docker:devel

# Install mariner from the build context rather than pulling a branch from
# GitHub, so the image matches the commit that built it. Previously this
# installed EricSDavis/mariner@dev unconditionally, which meant a build
# triggered by main still produced an image of dev.
#
# remotes is installed explicitly because it no longer ships with the base
# image -- that is what broke this build ("there is no package called
# 'remotes'"). BiocManager is still preinstalled, and its repositories are
# needed for the Bioconductor dependencies.
COPY . /tmp/mariner

RUN R -e 'install.packages("remotes"); \
          remotes::install_local("/tmp/mariner", dependencies = TRUE, \
                                 repos = BiocManager::repositories())' \
    && rm -rf /tmp/mariner

WORKDIR /home/rstudio
CMD ["/init"]
