# this is the main R script which generates all the results.
# Run this inside the docker container
system("mkdir -p ~/.cache/R/ExperimentHub")
system("git pull")
# look at the covariates
if (file.exists("asd_covariates.html") ) {
  message("asd_covariates.Rmd script completed")
} else {
  rmarkdown::render("asd_covariates.Rmd")
}

# limma analysis of blood data
rmarkdown::render("limma_blood.Rmd")

# limma analysis of guthrie card data
rmarkdown::render("limma_guthrie.Rmd")

# limma analysis of buccal swab data
rmarkdown::render("limma_buccal.Rmd")

# basic downstream analysis
# DMR calling and ORA pathway enrichment
rmarkdown::render("downstream1_ados.Rmd")

# deep downstream analysis
# pathway enrichment with mitch
rmarkdown::render("downstream2_ados.Rmd")
