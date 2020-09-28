
#'Fixed sample-size information for one-sided test
#'
#'The \code{ztest.I} function finds the information required of a one-sided Z-test with given power
#'and significance level
#'#ztest.I(delta=1, sig.level=0.05, power=0.9)
#'
#'@param delta numeric; targeted standardized effect size
#'@param sig.level numeric; one-sided significance level
#'@param power numeric; target power
#'@return numeric value with the (fractional) information achieving the target power
#'@author Aniko Szabo
#'@keywords internal
#'@importFrom stats qnorm
#'

ztest.I <- function(delta, sig.level, power){
  za <- qnorm(sig.level, lower.tail=FALSE)
  zb <- qnorm(power)
  I <- (za+zb)^2 /delta^2
  I
}
