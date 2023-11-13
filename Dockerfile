FROM bioconductor/bioconductor_docker:RELEASE_3_17

# Update apt-get
RUN apt-get update \
        && apt-get install -y nano git  libncurses-dev \
        ## Install the python package magic wormhole to send files
        && pip install magic-wormhole           \
        ## Remove packages in '/var/cache/' and 'var/lib'
        ## to remove side-effects of apt-get update
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*


# Install CRAN packages
RUN Rscript -e 'install.packages(c("beeswarm","data.table","dplyr","eulerr","ggplot2","gplots","kableExtra", "matrixStats", "RColorBrewer", "reshape2", "vioplot"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("DMRcate", "FlowSorted.Blood.450k", "FlowSorted.Blood.EPIC", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "impute", "limma", "minfi", "missMethyl", "mitch"))'

RUN Rscript -e 'install.packages(c("WGCNA"))'

RUN git clone https://github.com/markziemann/asd_meth.git

COPY ASD_EPIC_DATA /asd_meth

# Set the container working directory
ENV DIRPATH /asd_meth
WORKDIR $DIRPATH

