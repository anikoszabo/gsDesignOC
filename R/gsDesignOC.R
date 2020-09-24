
#'Find optimal group-sequential design
#'
#'The \code{gsDesignOC} function finds the analysis times, boundaries, and sample size
#'for a group-sequential trial
#'
#'@export
#'@param thA numeric; effect size under the alternative hypothesis
#'@param thA.seq monotone numeric vector of values getting closer to thA from outside the
#'th0-thA range (ie, if thA > th0, they should form a decreasing sequence with all values > thA). #'For the k-th value the study will stop for efficacy at or before the k-th stage with probability \code{power.efficacy}.
#'@param th0 numeric; effect size under the null hypothesis
#'@param th0.seq monotone numeric vector of values getting closer to th0 from outside the
#'th0-thA range (ie, if thA > th0, they should form an increasing sequence with all values < th0).
#'The study will stop for futility at or before the k-th stage with probability \code{power.futility}.
#'Its length should be equal to that of \code{thA.seq}. The default value of \code{NULL} implies no
#'futility monitoring.
#'@param sig.level numeric; the one-sided significance level. The default is 0.025.
#'@param power numeric; the desired study power. The default is 0.9.
#'@param power.efficacy numeric; the probability of stopping for efficacy at stage k under \code{thA.seq}
#'@param power.futility numeric; the probability of stopping for futility at stage k under \code{th0.seq}
#'@param futility.type character string, one of \code{c("none", "non-binding","binding")} or their
#'unique abbreviations. Specifies whether a futility boundary should not be used at all ("none"), or if it
#'is used, then whether the effect of the futility boundary should be included
#'in the efficacy power/type I error calculations ("binding"), or not ("non-binding").
#'@param mix.theta numeric vector; specifies the set of alternatives under which the design is optimized.
#' Defaults to \code{thA}.
#'@param mix.w numeric vector of length equal to that of \code{mix.theta}. The weights of the
#'elements of \code{mix.theta} in the optimality criterion. It will be normalized to a sum of 1.
#'Defaults to equal weights.
#'@param control list of parameters controlling the estimation altorithm to replace the default
#'values returned by the \code{glsControl} function.
#'@return an object of class \code{gsDesignOC}
#'@author Aniko Szabo
#'@references Szabo, A, Tarima, S (?) Operating-characteristic guided design of group-sequential trials.
#'@keywords nonparametric design
#'@examples
#'
#'gsDesignOC(0.3, thA.seq = c(1, 0.5))
#'
#'@name gsDesignOC

gsDesignOC <- function(thA, thA.seq, th0=0, th0.seq=NULL,
                       sig.level = 0.025,
                       power=0.9, power.efficacy=power, power.futility = power,
                       futility.type=c("none","non-binding","binding"),
                       mix.theta = thA,
                       mix.w = rep(1, length(mix.theta)),
                       control=list()){
  
    if (any(diff(c(thA.seq, thA)) * sign(th0-thA) <= 0))
      stop("'thA.seq' should approach thA in a monotone sequence outside the th0-thA range")

    if (!is.null(th0.seq)){
      if (length(th0.seq) != length(thA.seq))
        stop("If specified, 'th0.seq' should have the same length as 'thA.seq'")
      if (any(diff(c(th0.seq, th0)) * sign(thA-th0) <= 0))
        stop("'th0.seq' should approach th0 in a monotone sequence outside the th0-thA range")
    }

    if (length(mix.w) != length(mix.theta))
        stop("'mix.w' should have the same length as 'mix.theta'")

    if (power.efficacy > power)
      stop("The value of 'power.efficacy' should not exceed the value of 'power'.")

    if (power.futility > 1-sig.level)
      stop("The value of 'power.futility' should not exceed the value of 1-'sig.level'.")

    futility.type <- match.arg(futility.type)
    if (futility.type != "none"){
      if (is.na(th0.seq)) stop("`th0.seq` should be specified if a futility bound is requested")
    }

    controlvals <- ocControl()
    if (!missing(control))
      controlvals[names(control)] <- control

  

  # create skeleton gsDesignOC object
  res <- list(thA=thA, thA.seq=thA.seq,
              th0=th0, th0.seq=th0.seq,
              sig.level = sig.level,
              power=power, power.efficacy = power.efficacy, power.futility = power.futility,
              futility.type=futility.type)
  class(res) <- "gsDesignOC"

  
    # n.vec is stage-specific sample sizes
    .cp <- function(n.vec){

      n <- sum(n.vec)
      timing <- cumsum(n.vec/n)
      if (!is.null(sig.level_final)) timing <- head(timing, -1)

      res$n <- n
      res$timing <- timing
      res$bounds <- calc.bounds(res)
      dp <- oc(res, EN_theta=mix.theta, mix.w=mix.w)
      Q <-
        controlvals$optim.penalty * abs(dp$typeI - sig.level)/sqrt(sig.level*(1-sig.level)) +
        controlvals$optim.penalty * abs(dp$power - power)/sqrt(power*(1-power)) +
        dp$EN
      Q
    }
  

  k <- length(thA.seq) + 1
  n.guess <- ztest.n(delta=c(thA.seq,thA)-th0, sd=1,
                     sig.level=sig.level, power = power)
  thA.seq <- c(thA.seq, thA)

  oo <- optim(diff(c(0,n.guess)), fn=.cp)

  n <- sum(oo$par)
  timing <- cumsum(oo$par)/n
  if (!is.null(sig.level_final)) timing <- head(timing, -1)

  res$n <- n
  res$timing <- timing
  return(res)
}


#'Operating characteristics of a potential design
#'
#'The values supplied in the function call replace the defaults and a list with all possible #'arguments is returned. The returned list is used as the \code{control} argument to the
#'\code{gsDesignOC} function.
#'
#'@export
#'@param x object of class \code{gsDesignOC}
#'@param EN_theta numeric vector of parameter values for which expected sample size calculation is performed
#'@param mix.w numeric vector of positive weights for each value of \code{EN_theta}. Defaults to equal weights.
#'@return a list with the following elements:
#'\describe{
#'\item{typeI}{type I error, ie the probability of not crossing the efficacy boundary under theta=th0}
#'\item{power}{power, ie the probability of crossing the efficacy boundary under theta=thA}
#'\item{EN}{expected sample size under the requested alternatives}
#'}
#'@author Aniko Szabo
#'@keywords design
#'@examples
#'
#'

oc <- function(x, EN_theta=x$thetaA,  mix.w = rep(1, length(EN_theta))){

  if (length(mix.w) != length(EN_theta))
    stop("`EN_theta` and `mix.w` should have the same length")
  if (!all(mix.w > 0))
    stop("`mix.w` should have only positive elements")

  # under H0
  res0 <- overall.crossprob(x, theta=x$theta0)
  # under Ha
  resA <- overall.crossprob(x, theta=x$thetaA)

  # under EN-alternatives
  en_vec <- numeric(length(EN_theta))
  for (i in seq_along(EN_theta)){
    th <- EN_theta[i]
    # check if we already calculated the EN at this value
    if (abs(th - x$theta0) < .Machine$double.eps){
      en_vec[i] <- res0$EN
    } else if (abs(th - x$thetaA) < .Machine$double.eps) {
      en_vec[i] <- resA$EN
    } else {
      en_vec[i] <- overall.crossprob(x, theta=th)
    }
  }
  en <- en_vec %*% mix.w / sum(mix.w)
  list(typeI = res0$cross.up,
       power = resA$cross.up,
       EN = en)
}

#'Control values for gsDesignOC
#'
#'The values supplied in the function call replace the defaults and a list with all possible
#'arguments is returned. The returned list is used as the \code{control} argument to the
#'\code{gsDesignOC} function.
#'
#'@export
#'@param optim.penalty numeric; relative weight of terms ensuring target type I and type II
#'errors versus sample size in optimization routine
#'@return a list with components for each of the possible arguments
#'@author Aniko Szabo
#'@keywords design
#'@examples
#'
#'ocControl(optim.penalty = 100)

ocControl <- function(optim.penalty = 1000){
  list(optim.penalty = optim.penalty)
}

#'Calculate the information times and efficacy/futility boundary values given the alpha-spending sequence and desired operating characteristics

#'@export
#'@param x an object of class \code{gsDesignOC}

calc.bounds <- function(x, alpha.seq){

  n.stages <- length(x$thA.seq) + 1
  alpha.cum <- cumsum(alpha.seq)
  ub <- rep(20, n.stages)
  lb <- rep(-20, n.stages)
  ivec <- rep(NA, n.stages)


  
  exc <- function(u, I, stage,  theta, target){
    gg <- gsDesign::gsProbability(k=stage,
                                  theta=theta,
                                  n.I =c(ivec[1:(stage-1)],I),
                                  a = c(lb[1:(stage-1)], -20),
                                  b = c(ub[1:(stage-1)], u))
     sum(gg$upper$prob[1:stage]) - target
  }
  
  
  uI <- function(I, stage){
    res <- uniroot(exc, interval=c(qnorm(x$sig.level, lower=FALSE), 20),
                   I = I, stage=stage, theta=x$th0, target = alpha.cum[stage],
                   extendInt = "downX")
   res$root
  }
  
  
  exc_low <- function(l, stage, theta, target){
    gg <- gsDesign::gsProbability(k=stage,
                                  theta=theta,
                                  n.I =ivec[1:stage],
                                  a = c(head(lb, stage-1), l),
                                  b = ub[1:stage])
     sum(gg$lower$prob[1:stage]) - target
  }
  
  
  lI <- function(stage, theta){
    res <- uniroot(exc_low, interval=c(-20, qnorm(x$sig.level, lower=FALSE)),
                   stage=stage, theta=theta, target = x$power.futility,
                   extendInt = "upX")
   res$root
  }
  


  ivec[1] <- ztest.I(delta=x$thA.seq[1]-x$th0, sig.level=alpha.seq[1], power=x$power.efficacy)
  ub[1] <- qnorm(alpha.seq[1], lower=FALSE)


  if (n.stages > 1){
    for (k in 2:n.stages){
    
      if (x$futility.type == "binding"){
      
        lb[k-1] <- lI(stage=k-1, theta=x$th0.seq[k-1])
      
      }

      
      if (k == n.stages){
        power.target <- x$power
        .th <- x$thA
      } else {
        power.target <- x$power.efficacy
        .th <- x$thA.seq[k]
      }

      minI <- max(ztest.I(delta = .th - x$th0, power=power.target, sig.level=alpha.cum[k]),
                  ivec[k-1])
      resI <- uniroot(function(ii)exc(uI(ii, k), ii, stage=k, theta=.th,
                                      target = power.target),
                      interval=c(minI, 2*minI),
                      extendInt = "upX")
      ivec[k] <- resI$root
      ub[k] <- lb[k] <- uI(resI$root, k)
      

    
    }
  }

  if (x$futility.type == "non-binding"){
    
      for (k in 1:(n.stages-1)){
        lb[k] <- lI(stage=k, theta=x$th0.seq[k])
      }
    
  }

  x$upper <- ub
  x$lower <- lb
  x$info <- ivec
  return(x)
}

