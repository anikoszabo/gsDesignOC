library(devtools)
library(covr)
source('nuweb/Nuweb.R')
gs <- as.package("../gsDesignOC")


# initial setup
#use_testthat()
#use_readme_rmd()

nuweb(gs)
shell("cd c:/gsDesignOC/nuweb/ && texify --pdf --quiet --run-viewer gsDesignOC.tex")

document(gs)
run_examples(gs)
load_all(gs)

test(gs)
cov <- package_coverage(gs$path)
report(cov)

build_readme()

check(gs, check_dir = "c:/Temp", cran = TRUE, manual=TRUE)
install(gs)

# make a gsDesign object
library(gsDesign)
g <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "non-binding")
gg <- as.gsDesign(g)

g0 <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "binding")
g1 <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "none")

g2 <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "non-binding")
gg2 <- as.gsDesign(g2)
