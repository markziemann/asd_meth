
# look at the covariates
rmarkdown::render("asd_covariates.Rmd")

# limma analysis of blood data
rmarkdown::render("limma_blood.Rmd")

# limma analysis of guthrie card data
rmarkdown::render("limma_guthrie.Rmd")

# compartment enrichment
ml <- list.files(".",pattern="downstream1")
ml <- ml[grep("Rmd",ml)]
lapply(ml,rmarkdown::render)

# pathway enrichment 2
ml <- list.files(".",pattern="downstream2")
ml <- ml[grep("Rmd",ml)]
lapply(ml,rmarkdown::render)
