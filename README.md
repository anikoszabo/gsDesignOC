
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gsDesignOC

<!-- badges: start -->

<!-- badges: end -->

The goal of gsDesignOC is to construct optimal group-sequential designs
that maintain pre-specified operating characteristics. You can browse
its [source code](https://github.com/anikoszabo/gsDesignOC).

## Installation

In the future, this package will be released on CRAN, but currently the
development version is available on [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anikoszabo/gsDesignOC")
```

## Basic usage

``` r
library(gsDesignOC)
g <- gsDesignOC(n.stages = 2, rE.seq = c(2,1), rF.seq = c(-1,0), r_EN = c(0, 1),
                sig.level = 0.05, power = 0.8, power.efficacy = 0.8, power.futility = 0.9,
                futility.type = "non-binding")
g
#> $n.stages
#> [1] 2
#> 
#> $rE.seq
#> [1] 2 1
#> 
#> $rF.seq
#> [1] -1  0
#> 
#> $sig.level
#> [1] 0.05
#> 
#> $n.fix
#> [1] 1
#> 
#> $power
#> [1] 0.8
#> 
#> $power.efficacy
#> [1] 0.8
#> 
#> $power.futility
#> [1] 0.9
#> 
#> $futility.type
#> [1] "non-binding"
#> 
#> $upper
#> [1] 2.411574 1.680736
#> 
#> $lower
#> [1] -0.3450637  1.6807356
#> 
#> $info
#> [1] 2.645819 6.295366
#> 
#> $spending
#> [1] 0.007941925 0.042058075
#> 
#> $n
#> [1] 0.427949 1.018246
#> 
#> attr(,"class")
#> [1] "gsDesignOC"
```

## Extended Example

### Group sequential trial without gsDesignOC

We want to design a trial to test the hypotheses
\(H_0:\mu_1 \leq \mu_2\) at a one-sided 5% significance level, so that
80% power is achieved when \(\mu_1 - \mu_2 = 0.5\) (here \(\mu_i\) is
the mean of group \(i\)). We want to add two interim analyses, so that
we can stop the trial early if the difference between the groups is
larger than expected. The `gsDesign` package (among others) provides
tools to construct such group-sequential trials, however one needs to
specify the timing of the interim analyses and an error-spending
function. Different choices give different sample sizes and decision
cutoffs, however there is little guidance on how to make these choices.

First, lets compute the required sample size without interim analyses
(assuming equal standard deviation of 1 in each group). Since
`power.t.test` returns the per-group sample size, we will multiply the
result by 2 to get total sample size.

``` r
n0 <- power.t.test(delta=0.5, sd=1, alternative = "one.sided", sig.level=0.05, power=0.8)$n * 2
n0
#> [1] 100.3016
```

We will construct two group-sequential designs: one with the default
settings, and one with different custom settings.

``` r
library(gsDesign)
#> Loading required package: ggplot2
g1 <- gsDesign(k=3, test.type=1, alpha=0.05, beta=0.2, n.fix=n0)
g1
#> One-sided group sequential design with
#> 80 % power and 5 % Type I Error.
#>               
#>   Analysis  N   Z   Nominal p  Spend
#>          1  34 2.79    0.0026 0.0026
#>          2  68 2.29    0.0110 0.0099
#>          3 102 1.68    0.0465 0.0375
#>      Total                    0.0500 
#> 
#> ++ alpha spending:
#>  Hwang-Shih-DeCani spending function with gamma = -4.
#> 
#> Boundary crossing probabilities and expected sample size
#> assume any cross stops the trial
#> 
#> Upper boundary (power or Type I Error)
#>           Analysis
#>    Theta      1      2      3 Total  E{N}
#>   0.0000 0.0026 0.0099 0.0375  0.05 101.3
#>   0.2483 0.0889 0.3224 0.3886  0.80  84.8
```

``` r
g2 <- gsDesign(k=3, test.type=1, alpha=0.05, beta=0.2, n.fix=n0, sfu="Pocock", timing=c(0.5, 0.7))
g2
#> One-sided group sequential design with
#> 80 % power and 5 % Type I Error.
#>               
#>   Analysis  N   Z   Nominal p  Spend
#>          1  58 1.95    0.0257 0.0257
#>          2  81 1.95    0.0257 0.0128
#>          3 115 1.95    0.0257 0.0115
#>      Total                    0.0500 
#> 
#> ++ alpha spending:
#>  Pocock boundary.
#> 
#> Boundary crossing probabilities and expected sample size
#> assume any cross stops the trial
#> 
#> Upper boundary (power or Type I Error)
#>           Analysis
#>    Theta      1      2      3 Total  E{N}
#>   0.0000 0.0257 0.0128 0.0115  0.05 112.7
#>   0.2483 0.4727 0.1729 0.1544  0.80  81.6
```

Recall that we wanted to the interim analysis to stop early if the
actual difference is larger than hypothesized, but we never had to
define what we actually mean by that. Let’s say we think that the effect
could be up to 2-fold higher (\(\mu_1-\mu_2 = 1\)) or, more
realistically, 1.5-fold higher (\(\mu_1-\mu_2 = 0.75\)). We can look at
the probability of stopping by each analysis (ie at or before) for all
the alternatives of interest:

``` r
cumprob <- function(design, multipliers){
  p <- gsProbability(d=design, theta=design$delta * multipliers)
  cum_p <- apply(p$upper$prob, 2, cumsum)
  dimnames(cum_p) <- list("Analysis"=1:nrow(cum_p), "Multiplier"= multipliers)
  list(p=cum_p, EN=setNames(p$en, multipliers))}
(cp1 <- cumprob(design=g1, multipliers = c(0,1,1.5,2)))
#> $p
#>         Multiplier
#> Analysis           0          1       1.5         2
#>        1 0.002606123 0.08894838 0.2662871 0.5394692
#>        2 0.012492890 0.41135911 0.7863006 0.9651068
#>        3 0.050000001 0.80000000 0.9817161 0.9995825
#> 
#> $EN
#>         0         1       1.5         2 
#> 101.30246  84.83531  66.09186  50.75215
(cp2 <- cumprob(design=g2, multipliers = c(0,1,1.5,2)))
#> $p
#>         Multiplier
#> Analysis          0         1       1.5         2
#>        1 0.02570173 0.4727318 0.8082434 0.9649518
#>        2 0.03849209 0.6456199 0.9285302 0.9948104
#>        3 0.04999998 0.8000000 0.9837051 0.9997151
#> 
#> $EN
#>         0         1       1.5         2 
#> 112.72865  81.59832  64.17558  58.30298
```

We can see that both of our designs stop by the 3rd (final) analysis
with probability 5% and 80% for the null hypothesis (multiplier 0) and
alternative hypothesis (multiplier 1), respectively. However, `g1` has a
fairly low probability of stopping at the first interim analysis even if
the true effect is twice as large as the design effect, while `g2` has
high stopping probabilities for even the 1.5-fold larger effect at the
first interim, and te second interim analysis seems superfluous.
Notably, it also has a larger sample size.

### Optimizing the trial with gsDesignOC

The operating-characteristic guided design in `gsDesignOC` provides an
alternative approach that pre-specifies the probability of stopping by
each interim analysis for a sequence of target alternatives, and then
finds the design satisfying these operating-characteristic constraints
that has the smallest expected sample size.

Here, we will build a design that stops with 80% probability at stage 1
or 2, respectively, when the true effect is 2-fold or 1.5-fold higher
than the target effect. Among all the designs with this property, we
will select the one with the smallest average of expected sample sizes
under the original null and alternative hypotheses.

``` r
library(gsDesignOC)
g <- gsDesignOC(n.stages = 3, rE.seq=c(2,1.5,1), sig.level=0.05, power=0.8, n.fix=n0,
                power.efficacy = 0.8, r_EN = c(0,1))
gg <- as.gsDesign(g)
gg
#> One-sided group sequential design with
#> 80 % power and 5 % Type I Error.
#>               
#>   Analysis  N   Z   Nominal p  Spend
#>          1  45 2.47    0.0068 0.0068
#>          2  67 2.24    0.0125 0.0090
#>          3 103 1.70    0.0443 0.0342
#>      Total                    0.0500 
#> 
#> ++ alpha spending:
#>  Piecewise linear spending function with line points = 0.431607983655741 0.649699912404647 1 0.136037087582927 0.316445695166114 1.
#> 
#> Boundary crossing probabilities and expected sample size
#> assume any cross stops the trial
#> 
#> Upper boundary (power or Type I Error)
#>           Analysis
#>    Theta      1      2      3 Total  E{N}
#>   0.0000 0.0068 0.0090 0.0342  0.05 102.2
#>   0.2483 0.2081 0.2271 0.3648  0.80  82.5
```

We can see that the timing of the analyses is selected automatically,
and it is different from both previous designs.

The following summary shows that the target stopping probabilities are
reached, and the average of the first two expected sample sizes is
smaller than for `g1` and `g2`, even though `g1` does not reach the
target interim analysis probabilities at the first stage (and `g2`
overshoots them greatly).

``` r
(cpopt <- cumprob(gg, multipliers = c(0, 1, 1.5, 2)))
#> $p
#>         Multiplier
#> Analysis           0         1       1.5         2
#>        1 0.006801854 0.2081210 0.5057218 0.8000065
#>        2 0.015822303 0.4352238 0.8000038 0.9678276
#>        3 0.050000023 0.8000000 0.9819091 0.9995966
#> 
#> $EN
#>         0         1       1.5         2 
#> 102.18227  82.54563  62.71714  50.06290
c(g1 = mean(cp1$EN[1:2]), g2 = mean(cp2$EN[1:2]), gg = mean(cpopt$EN[1:2]))
#>       g1       g2       gg 
#> 93.06888 97.16348 92.36395
```
