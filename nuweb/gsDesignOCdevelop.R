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

g0 <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "none")
g1 <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "non-binding")
g2 <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "binding")

gg0 <- as.gsDesign(g0)
gg1 <- as.gsDesign(g1)
gg2 <- as.gsDesign(g2)

skeleton <- list(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0),
                 power=0.9, power.futility=0.9, power.efficacy=0.9, sig.level=0.025)
b0 <- calc.bounds(c(skeleton,  futility.type = "none"), alpha.seq=c(0.005, 0.01, 0.01))
b1 <- calc.bounds(c(skeleton,  futility.type = "non-binding"), alpha.seq=c(0.005, 0.01, 0.01))
b2 <- calc.bounds(c(skeleton,  futility.type = "binding"), alpha.seq=c(0.005, 0.01, 0.01))

