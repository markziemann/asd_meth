# asd_meth

Analysis of methylation data of twins where one or both of the twins has autism spectrum disorder.
This is a longitudinal study that includes newborm blood spots, blood at ASD assessment, and
buccal swabs in later life.

## Requirements

Ubuntu 20+ with Docker installed.
64GB or more RAM.
8 CPU threads or more.

## Steps to run

1. Pull the docker image which has all the necessary software and data:

```
docker pull mdziemann/asd_meth
```

2. Start up an interactive shell in the docker container:

```
docker run -it asd_meth /bin/bash
```

3. Run the master script:

```
Rscript main.R
```

The master script will subsequently execute 6 R markdown scripts:

* asd_covariates.Rmd

* limma_blood.Rmd, limma_buccal.Rmd, limma_guthrie.Rmd

* downstream1_ados.Rmd

* downstream2_ados.Rmd

4. Copy the analysis directory from the container to the host machine.
It will contain all the results in the HTML and PDF files which includes figures and tables.

```
docker cp $(docker ps -alq):/asd_meth .
```

## Reporting problems

Please use github issues to report bugs.
Include any error messages or other relevant information to facilitate troublieshooting.
