% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsDesignOC.R
\name{oc}
\alias{oc}
\title{Operating characteristics of a potential design}
\usage{
oc(x, r_EN = x$r_EN, r_EN.w = x$r_EN.w)
}
\arguments{
\item{x}{object of class \code{gsDesignOC}}

\item{r_EN}{numeric vector of parameter values for which expected sample size calculation is performed. Defaults to values used to construct design.}

\item{r_EN.w}{numeric vector of positive weights for each value of \code{r_EN}. Defaults to equal weights. Defaults to values used to construct design.}
}
\value{
a list with the following elements:
\describe{
\item{ave.EN}{numeric; weighted average of expected sample sizes under the requested alternatives}
\item{EN.vec}{numeric vector; expected sample sizes under each of the requested alternatives}
\item{efficacy.cumcross}{numeric vector; cumulative probability of stopping for efficacy at or before stage k
 under rE.seq[k], including 1 at the end.}
\item{futility.cumcross}{numeric vector}{numeric vector; cumulative probability of stopping for
 futility at or before stage k under rF.seq[k], including 0 at the end.}
}
}
\description{
Calculates the average expected information under a weighted set of values for the alternative hypothesis, and
the probability of stopping for futility or effecacy under the design alternatives. If the futility boundary
is non-binding, then the lower boundary is not included in the calculation of the expected sample size and
efficacy stopping probabilities, but it is included in the futility stopping calculations.
}
\examples{
g <- gsDesignOC(n.stages = 2, rE.seq = c(1.5, 1), rF.seq = -1, power.efficacy=0.8,
          power.futility=0.8, power=0.9,
          futility.type = "non-binding", r_EN=c(1.5, 1, 0))
oc(g)

oc(g, r_EN = c(2,1), r_EN.w = c(0.5,0.5))

}
\author{
Aniko Szabo
}
\keyword{design}
