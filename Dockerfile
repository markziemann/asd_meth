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


# problematic package bug
RUN git clone https://github.com/bmbolstad/preprocessCore.git && \
  cd preprocessCore && \
  R CMD INSTALL --configure-args="--disable-threading" . && \
  cd ..

# BiocFileCache is also problematic
RUN git clone https://github.com/Bioconductor/BiocFileCache.git && \
  cd BiocFileCache && \
  R CMD INSTALL --configure-args="--disable-threading" . && \
  cd ..

# Install CRAN packages
RUN Rscript -e 'install.packages(c("beeswarm","data.table","dplyr","eulerr","ggplot2","gplots","kableExtra", "matrixStats", "RColorBrewer", "reshape2", "vioplot","cellranger","readxl","rematch"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("DMRcate","FlowSorted.Blood.450k", "FlowSorted.Blood.EPIC", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "impute", "limma", "minfi", "missMethyl", "mitch","DMRcatedata"))'

# WGCNA separately
RUN Rscript -e 'install.packages(c("WGCNA"))'

# get a clone of the codes
RUN git clone https://github.com/markziemann/asd_meth.git

# copy in the EPIC array data
COPY ASD_EPIC_DATA /asd_meth/ASD_EPIC_DATA
COPY cp_asd_data/204375410103 /asd_meth/cp_asd_data/204375410103
COPY cp_asd_data/204375410107 /asd_meth/cp_asd_data/204375410107

# Set the container working directory
ENV DIRPATH /asd_meth
WORKDIR $DIRPATH

