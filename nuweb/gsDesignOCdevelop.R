library(devtools)
library(covr)
source('nuweb/Nuweb.R')
gs <- as.package("../gsDesignOC")


# initial setup
#use_testthat()
use_readme_rmd()

nuweb(gs)
shell("cd c:/gsDesignOC/nuweb/ && texify --pdf --quiet --run-viewer gsDesignOC.tex")

document(gs)
run_examples(gs)
load_all(gs)

test(gs)
cov <- package_coverage(gs$path)
report(cov)


check(gs, check_dir = "c:/Temp", cran = TRUE, manual=TRUE)
install(gs)

# make a gsDesign object
library(gsDesign)
g <- gsDesignOC(n.stages=3, rE.seq = c(2,1.5,1), rF.seq=c(-1,-0.5,0), futility.type = "binding")

g2 <- gsDesign(k=3,test.type=4)

gg <- list(
  k = g$n.stages,
  test.type = 3, #binding
  alpha = g$sig.level,
  beta = 1-g$power,
  astar = 0,
  delta = (qnorm(g$sig.level, lower=FALSE) + qnorm(g$power) )/g$n.fix,
  n.fix = g$n.fix,
  timing = g$info/max(g$info),
  r=18,
  n.I = g$n,
  maxn.IPlan = 0,
  nFixSurv = 0,
  nSurv = 0,
  endpint = NULL,
  delta1 = 1,
  delta0 = 0,
  upper = list(name="Optimized",
               param=g$rE.seq,
               parname="r(efficacy)",
               spend=g$spending,
               bound=g$upper),
  lower = list(name="Optimized",
               param=g$rF.seq,
               parname="r(futility)",
               bound=g$lower,
               spend=g$spending),
  theta = c(0, (qnorm(g$sig.level, lower=FALSE) + qnorm(g$power) )/g$n.fix),
  en = oc(g, r_EN=c(0,1))$EN.vec
)
class(gg) <- c("gsDesign")
