
#'One-sample Z-based sample size for one-sided test
#'
#'The \code{ztest.n} function finds the sample size of a one-sided Z-test with given power
#'and significance level
#'
#'@param delta numeric; targeted effect size
#'@param sd numeric; standard deviation, assumed known
#'@param sig.level numeric; one-sided significance level
#'@param power numeric; target power
#'@return numeric value with the (fractional) sample size achieving the target power
#'@author Aniko Szabo
#'@keywords internal
#'@examples
#'
#'ztest.n(delta=1, sd=1, sig.level=0.05, power=0.9)

ztest.n <- function(delta, sd, sig.level, power){
  za <- qnorm(sig.level, lower=FALSE)
  zb <- qnorm(power)
  n <- (za+zb)^2 * sd^2/delta^2
  n
  }

#'Conversion between boundary and nominal significance level
#'
#'The \code{nom.to.bnd} and \code{bnd.to.nom} functions perform conversion between the boundary and
#' matching significance level.
#'
#'@param nominal.level numeric vector or matrix with two rows of the nominal significance level for each stage.
#'  When a matrix, row 1 corresponds to the efficacy boundary and row 2 to the futility. Defaults to NULL.
#'@param cutoffs numeric vector or numeric matrix with two rows of the z-test cutoffs for each stage.
#'  When a matrix, row 1 corresponds to the efficacy boundary and row 2 to the futility. Defaults to NULL.
#'@param n integer vector of sample sizes of each stage. Its sum is the total study sample size.
#'@param theta.null numeric, the null hypothesis being tested
#'@param theta.alt numeric, the alternative hypothesis being tested for the futility bound calculation
#'@keywords internal
#'@examples
#'
#'## only efficacy
#'(b1 <- nom.to.bnd(nominal.level = c(0.01, 0.02), n = c(30, 50)))
#'bnd.to.nom(b1, n=c(30,50))

nom.to.bnd <- function(nominal.level, n, theta.null, theta.alt=NULL){
  nominal.level <- rbind(nominal.level) # raises to matrix if one row
  m <- nrow(nominal.level)
  k <- ncol(nominal.level)
  cutoffs <- matrix(NA, nrow=m, ncol=k)
  cutoffs[1, ] <- qnorm(nominal.level[1,], lower=FALSE)
  if (k > 1) {
    # flip direction + shift hypothesis
    cutoffs[2, ] <- qnorm(nominal.level[2,], lower=TRUE) + sqrt(n)*(theta.alt -theta.null)
    }
  cutoffs[1:k,]  # drops to vector if one row
}

#'@describeIn nom.to.bnd Conversion from nominal significance level to boundary
#'@export
#'@keywords internal
#'@examples
#'## both efficacy and futility
#'(b2 <- nom.to.bnd(nominal.level = rbind(c(0.01, 0.02), c(0.05,0.1)),
#'                 n = c(30, 50)))
#'bnd.to.nom(b2, n=c(30,50))

bnd.to.nom <- function(cutoffs, n, theta.null, theta.alt=NULL){
  cutoffs <- rbind(cutoffs) # raises to matrix if one row
  m <- nrow(cutoffs)
  k <- ncol(cutoffs)
  noms <- matrix(NA, nrow=m, ncol=k)
  noms[1,] <- pnorm(cutoffs[1,], lower=FALSE)
  if (k > 1) {
    # flip direction + shift hypothesis
    noms[2,] <- pnorm(cutoffs[2,] - sqrt(n)*(theta.alt-theta.null), lower=TRUE)
    }
  noms[1:k,] # drops to vector if one row
}

