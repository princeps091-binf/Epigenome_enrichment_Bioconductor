library(renv)

renv::init(bioconductor = T)

renv::install("tidyverse")
renv::install("furrr")
renv::install("Matrix")
renv::install("viridis")
renv::install("seriation")
renv::install("svglite")
renv::install("Rcpp")
renv::install("r-lib/later")
renv::install("bioc::GenomicRanges")
renv::install("bioc::liftOver")
renv::install("bioc::AnnotationHub")
renv::install("bioc::rtracklayer")
renv::snapshot()

