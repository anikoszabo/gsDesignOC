
#'One-sample Z-based sample size for one-sided test
#'
#'The \code{ztest.n} function finds the sample size of a one-sided Z-test with given power
#'and significance level
#'
#'@export
#'@param delta numeric; targeted effect size
#'@param sd numeric; standard deviation, assumed known
#'@param sig.level numeric; one-sided significance level
#'@param power numeric; target power
#'@return numeric value with the (fractional) sample size achieving the target power
#'@author Aniko Szabo
#'@references Szabo, A, Tarima, S (?) Operating-characteristic guided design of group-sequential trials.
#'@keywords
#'@examples
#'
#'ztest.n(delta=1, sd=1, sig.level=0.05, power=0.9)

ztest.n <- function(delta, sd, sig.level, power){
  za <- qnorm(sig.level, lower=FALSE)
  zb <- qnorm(power)
  n <- (za+zb)^2 * sd^2/delta^2
  n
  }
