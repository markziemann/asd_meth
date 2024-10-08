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

## Tracing scripts to figures

All the figure panels and tables from the paper.

| Data |  output |  script |
| --- | --- | --- |
| Fig 1A | no PDF, get the png from limma_blood.html and limma_guthrie.html | limma_blood.Rmd limma_guthrie.Rmd |
| Fig 1B | pca_blood.pdf and pca_guthrie.pdf | limma_blood.Rmd limma_guthrie.Rmd |
| Fig 2 | manhat_limma_bl.pdf manhat_limma_gu.pdf | downstream1_ados.Rmd |
| Fig 3AD | hist_bl_ados.pdf hist_gu_ados.pdf | downstream2_ados.Rmd |
| Fig 3EF | genedens_ados.pdf |  downstream2_ados.Rmd |
| Table 2 | downstream2_ados.html | downstream2_ados.Rmd |
| Fig 4A | blgu_mitch_go_ados.pdf | downstream2_ados.Rmd |
| Fig 4B | adosbarplot_hypo_go.pdf | downstream2_ados.Rmd |
| Fig 4C | blgu_mitch_go_ados.pdf | downstream2_ados.Rmd |
| Fig 4D | bitter.pdf | downstream2_ados.Rmd |

## Reporting problems

Please use github issues to report bugs.
Include any error messages or other relevant information to facilitate troublieshooting.
