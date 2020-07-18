library(devtools)
library(covr)
source('z:/RForge/Nuweb.R')
gs <- as.package("../gsDesignOC")


# initial setup
use_testthat()

nuweb(gs)
document(gs)
run_examples(gs)  # or dev_example("ran.CMData")
load_all(gs)

test(gs)
cov <- package_coverage(gs$path)
shine(gs)


check(gs, check_dir = "c:/RForge", check_version = TRUE, cran = TRUE, manual=TRUE)
install(gs)
