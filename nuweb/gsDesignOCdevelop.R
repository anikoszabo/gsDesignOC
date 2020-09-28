library(devtools)
library(covr)
source('nuweb/Nuweb.R')
gs <- as.package("../gsDesignOC")


# initial setup
#use_testthat()

nuweb(gs)
shell("cd c:/gsDesignOC/nuweb/ && texify --pdf --quiet --run-viewer gsDesignOC.tex")

document(gs)
run_examples(gs)
load_all(gs)

test(gs)
cov <- package_coverage(gs$path)
shine(cov)


check(gs, check_dir = "c:/Temp", cran = TRUE, manual=TRUE)
install(gs)
