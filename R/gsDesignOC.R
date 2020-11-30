
#'Find optimal group-sequential design
#'
#'The \code{gsDesignOC} function finds the analysis times, boundaries, and sample size
#'for a group-sequential trial
#'
#'@export
#'@param n.stages integer; number of stages planned. One stage implies no interim analysis.
#'@param rE.seq monotone decreasing numeric vector of length \code{n.stages} ending with 1.
#' If it has length \code{n.stages-1}, a 1 will be appended to the end.
#' For the k-th interim value the study will stop for efficacy at or before the k-th stage with probability
#' \code{power.efficacy} if the true effect is \code{rE.seq[k]} times the effect under the alternative hypothesis.
#'@param rF.seq monotone increasing numeric vector of length \code{n.stages} ending with 0.
#' If it has length \code{n.stages-1}, a 0 will be appended to the end.
#' The study will stop for futility at or before the k-th stage with probability \code{power.futility}.
#' if the true effect is \code{rE.seq[k]} times the effect under the alternative hypothesis, ie
#' the true effect is actually in the opposite direction of the hypothesized effect.
#' The default value of \code{NULL} implies no futility monitoring.
#'@param n.fix numeric; the sample size for a one-stage design without interim tests. Defaults to 1,
#' in which case the resulting sample sizes can be interpreted as being relative to the single-stage study sample size.
#'@param sig.level numeric; the one-sided significance level. The default is 0.025.
#'@param power numeric; the desired study power. The default is 0.9.
#'@param power.efficacy numeric; the probability of stopping for efficacy at stage k under \code{rE.seq}
#'@param power.futility numeric; the probability of stopping for futility at stage k under \code{rF.seq}
#'@param futility.type character string, one of \code{c("none", "non-binding","binding")} or their
#'unique abbreviations. Specifies whether a futility boundary should not be used at all ("none"), or if it
#'is used, then whether the effect of the futility boundary should be included
#'in the efficacy power/type I error calculations ("binding"), or not ("non-binding").
#'@param r_EN numeric vector; specifies the set of alternatives under which the design is optimized.
#' It is interpreted multiplicatively compared to the design alternative hypothesis.
#' Defaults to 1, implying minimization under the alternative hypothesis.
#'@param r_EN.w numeric vector of length equal to that of \code{r_EN}. The weights of the
#'elements of \code{r_EN} in the optimality criterion. It will be normalized to a sum of 1.
#'Defaults to equal weights.
#'@param control list of parameters controlling the estimation altorithm to replace the default
#'values returned by the \code{glsControl} function.
#'@return an object of class \code{gsDesignOC}
#'@author Aniko Szabo
#'@references Szabo, A, Tarima, S (?) Operating-characteristic guided design of group-sequential trials.
#'@keywords nonparametric design
#'@examples
#'
#'gsDesignOC(2, rE.seq = c(1.5,1), rF.seq = c(-0.5,0),
#'           power.efficacy=0.8, power.futility=0.8, power=0.9,
#'           futility.type = "non-binding")
#'
#'@name gsDesignOC

gsDesignOC <- function(n.stages, rE.seq, rF.seq=NULL, n.fix=1,
                       sig.level = 0.025, power=0.9,
                       power.efficacy=power, power.futility = power,
                       futility.type=c("none","non-binding","binding"),
                       r_EN = 1,
                       r_EN.w = rep(1, length(r_EN)),
                       control=ocControl()){
  
    if (length(rE.seq) == n.stages - 1){
      rE.seq <- c(rE.seq, 1)
    } else if (length(rE.seq) == n.stages){
      if (!isTRUE(all.equal(rE.seq[n.stages], 1)))
        stop("Invalid input: When 'rE.seq' has length `n.stages`, its last element should equal 1.")
    } else {
      stop("Invalid input: `rE.seq` should have length equal to `n.stages` or `n.stages-1`.")
    }
    if (any(diff(rE.seq) >= 0))
      stop("Invalid input: `rE.seq` should form a decreasing sequence")

    if (length(rF.seq) == n.stages - 1){
      rF.seq <- c(rF.seq, 0)
    } else if (length(rF.seq) == n.stages){
      if (!isTRUE(all.equal(rF.seq[n.stages], 0)))
        stop("Invalid input: When `rF.seq` has length `n.stages`, its last element should equal 0.")
    } else if (!is.null(rF.seq)) {
      stop("Invalid input: When not NULL, `rF.seq` should have length equal to `n.stages` or `n.stages-1`.")
    }
    if (any(diff(rF.seq) <= 0))
      stop("Invalid input: `rF.seq` should form an increasing sequence")

    if (length(r_EN.w) != length(r_EN))
        stop("Invalid input: `r_EN.w` should have the same length as `r_EN`")
    if (any(r_EN.w < 0))
        stop("Invalid input: `r_EN.w` should have only positive values")

    if (power >= 1 || power <= sig.level)
      stop("Invalid input: `power` should be between `sig.level` and 1 (exclusive).")
    if (power.efficacy > power || power.efficacy <= sig.level)
      stop("Invalid input: The value of `power.efficacy` should be between `sig.level` and `power`.")

    if (power.futility > 1-sig.level || power.futility <= 0)
      stop("Invalid input: The value of `power.futility` should be between 0 and 1-`sig.level`.")

    futility.type <- match.arg(futility.type)
    if (futility.type != "none"){
      if (is.null(rF.seq)) stop("Invalid input: `rF.seq` should be specified if a futility bound is requested")
    }

    controlvals <- ocControl()
    if (!missing(control))
      controlvals[names(control)] <- control

  

  # create skeleton gsDesignOC object
  res <- list(n.stages = n.stages, rE.seq=rE.seq, rF.seq = rF.seq,
              sig.level = sig.level, n.fix = n.fix,
              power=power, power.efficacy = power.efficacy, power.futility = power.futility,
              futility.type=futility.type)
  class(res) <- "gsDesignOC"

  if (n.stages == 1){
    alpha.seq <- sig.level
  } else if (control$optim.method == "direct"){
  
    
      .cp <- function(y.vec){

        y <- c(y.vec, 0)  # add y_K
        alpha.seq <- sig.level * exp(y) / sum(exp(y))

        rescp <- calc.bounds(x=res, alpha.seq = alpha.seq)
        dp <- oc(rescp, r_EN=r_EN, r_EN.w=r_EN.w)
        Q <- dp$ave.EN
        Q
      }
    

    if (n.stages == 2){
      oo <- optimize(.cp, interval=c(-5,5))

      y.res <- oo$minimum
      alpha.seq <- sig.level * exp(c(y.res,0)) / sum(exp(c(y.res,0)))
    } else {
      y.start <- -log(seq(n.stages, 2, by=-1))

      oo <- optim(y.start, fn=.cp, method="Nelder-Mead")

      y.res <- oo$par
      alpha.seq <- sig.level * exp(c(y.res,0)) / sum(exp(c(y.res,0)))
    }
  
  } else if (control$optim.method == "dynamic"){
  

    .opt.spend <- function(x){
      if (x$n.stages == 1){
        return(x$sig.level)  # we spend it all
      }

      
        #create skeleton for a design with one fewer stages
        xx <- list(n.stages = x$n.stages - 1,
                   rE.seq = head(rE.seq, -1),
                   rF.seq = head(rF.seq, -1),
                   power=x$power.efficacy, power.efficacy = x$power.efficacy, power.futility = x$power.futility,
                   futility.type=x$futility.type)
        class(xx) <- "gsDesignOC"

      en <- function(y){
        delta.alpha <- exp(y) * x$sig.level
        cum.alpha <- x$sig.level - delta.alpha
        # figure out the optimal spending within the first K-1 stages
        xx$sig.level <- cum.alpha
        prev.spend <-.opt.spend (xx)
        # calculate design using these spending values
        xx2 <- calc.bounds(x, alpha.seq = c(prev.spend, delta.alpha))
        dp <- oc(xx2, r_EN = r_EN, r_EN.w = r_EN.w)
        en.res <- dp$ave.EN
        attr(en.res, "spending") <- c(prev.spend, delta.alpha)
        en.res
      }
      
        res.spend <- optimize(en, interval=c(-5, 0))
        opt.spending <- attr(res.spend$objective, "spending")
        return(opt.spending)
      
    }
  
    alpha.seq <- .opt.spend(x=res)
  
  } else {
    stop(sprintf("Optimization method %s is not available. Choose 'dynamic' or 'direct'.", control$optim.method))
  }

  res <- calc.bounds(x=res, alpha.seq)
  res$n <- n.fix * res$info /(qnorm(sig.level, lower.tail=FALSE) + qnorm(power))^2

  return(res)
}


#'Operating characteristics of a potential design
#'
#'Calculates the average expected information under a weighted set of values for the alternative hypothesis, and
#' the probability of stopping for futility or effecacy under the design alternatives. If the futility boundary
#' is non-binding, then the lower boundary is not included in the calculation of the expected sample size and
#' efficacy stopping probabilities, but it is included in the futility stopping calculations.
#'
#'@export
#'@param x object of class \code{gsDesignOC}
#'@param r_EN numeric vector of parameter values for which expected sample size calculation is performed
#'@param r_EN.w numeric vector of positive weights for each value of \code{r_EN}. Defaults to equal weights.
#'@return a list with the following elements:
#'\describe{
#'\item{ave.EN}{numeric; weighted average of expected sample sizes under the requested alternatives}
#'\item{EN.vec}{numeric vector; expected sample sizes under each of the requested alternatives}
#'\item{efficacy.cumcross}{numeric vector; cumulative probability of stopping for efficacy at or before stage k
#'  under rE.seq[k], including 1 at the end.}
#'\item{futility.cumcross}{numeric vector}{numeric vector; cumulative probability of stopping for
#'  futility at or before stage k under rF.seq[k], including 0 at the end.}
#'}
#'@author Aniko Szabo
#'@keywords design
#'@importFrom stats optimize optim
#'@importFrom utils head tail
#'@examples
#'g <- gsDesignOC(n.stages = 2, rE.seq = c(1.5, 1), rF.seq = -1, power.efficacy=0.8,
#'           power.futility=0.8, power=0.9,
#'           futility.type = "non-binding", r_EN=c(1.5, 1, 0))
#'oc(g, r_EN = c(1.5, 1, 0))
#'
#'

oc <- function(x, r_EN=1,  r_EN.w = rep(1, length(r_EN))){

  if (length(r_EN.w) != length(r_EN))
    stop("'r_EN' and 'r_EN.w' should have the same length")
  if (!all(r_EN.w > 0))
    stop("'r_EN.w' should have only positive elements")

  n.EN <- length(r_EN)
  n.A <- length(x$rE.seq)
  n.0 <- length(x$rF.seq)

  # crossing probabilities for futility
  ggF <- gsDesign::gsProbability(k = x$n.stages,
                                  theta = c(r_EN, x$rE.seq, x$rF.seq),
                                  n.I = x$info,
                                  a = x$lower,
                                  b = x$upper)

  # crossing probabilities for efficacy and EN
  if (x$futility.type == "non-binding") {
    #ignore lower boundary (except last stage) for EN & efficacy stopping
    ggE <- gsDesign::gsProbability(k = x$n.stages,
                                  theta = c(r_EN, x$rE.seq, x$rF.seq),
                                  n.I = x$info,
                                  a = c(rep(-20, x$n.stages-1), x$lower[x$n.stages]),
                                  b = x$upper)
   } else {
    ggE <- ggF
   }

  # expected sample size calculations
  en_vec <- ggE$en[1:n.EN]
  # rescale to n.fix
  n.fix <- if (is.null(x$n.fix)) 1 else x$n.fix
  en_vec <- n.fix * en_vec /(qnorm(x$sig.level, lower.tail=FALSE) + qnorm(x$power))^2

  # weighted average
  en <- c(en_vec %*% r_EN.w / sum(r_EN.w))

  # cumulative stopping probability calculations
  if (n.A > 0){
    p.up <- ggE$upper$prob[, n.EN + (1:n.A), drop=FALSE]
    cump.up <- apply(p.up, 2, cumsum)
    thA.cumcross <- setNames(diag(cump.up), x$rE.seq)
  } else {
    thA.cumcross <- NULL
  }

  if (n.0 > 0){
    p.low <- ggF$lower$prob[, n.EN + n.A + (1:n.0), drop=FALSE]
    cump.low <- apply(p.low, 2, cumsum)
    th0.cumcross <- setNames(diag(cump.low), x$rF.seq)
  } else {
    th0.cumcross <- NULL
  }

  list(ave.EN = en,
       EN.vec = setNames(en_vec, r_EN),
       efficacy.cumcross = thA.cumcross,
       futility.cumcross = th0.cumcross)
}

#'Control values for gsDesignOC
#'
#'The values supplied in the function call replace the defaults and a list with all possible
#'arguments is returned. The returned list is used as the \code{control} argument to the
#'\code{gsDesignOC} function.
#'
#'@export
#'@param optim.method character; defines the optimization method used: \code{"dynamic"} uses a recursive dynamic algorithm #' which selects the alpha-spending one stage at a time, \code{"direct"} uses the "Nelder-Mead" algorithm from
#' \code{\link{optim}} to simultaneously search over all possible alpha-spending sequences.
#'@return a list with components for each of the possible arguments
#'@author Aniko Szabo
#'@keywords design
#'@examples
#'
#'ocControl(optim.method = "direct")

ocControl <- function(optim.method = c("direct", "dynamic")){
  optim.method <- match.arg(optim.method)
  list(optim.method = optim.method)
}

#'Calculate the information times and efficacy/futility boundary values given the alpha-spending sequence and desired operating characteristics


#'@param x a list with desired operating characteristics. Should have \code{rE.seq}, \code{sig.level},
#' \code{power}, \code{power.efficacy}, and \code{futility.type} components; and \code{rF.seq} and \code{power.futility}
#' components if lower boundary is requested. Objects of class \code{gsDesignOC} (potentially incomplete)  have
#' these components.
#'@param alpha.seq numeric vector of stage-specific alpha-spending; values should add up to \code{x$sig.level}.
#'@return The value of \code{x} augmented with the following components:
#'\describe{
#'\item{info}{numeric vector; information times for the analyses}
#'\item{lower}{numeric vector; lower Z-score boundary for the futility decisions}
#'\item{upper}{numeric vector; upper Z-score boundary for the efficacy decisions}
#'\item{spending}{numeric vector}{the value of alpha-spending sequence \code{alpha.seq}}
#'}
#'@keywords internal
#'@importFrom stats uniroot qnorm
#'@importFrom utils head tail

calc.bounds <- function(x, alpha.seq){

  alpha.cum <- cumsum(alpha.seq)
  ub <- rep(20, x$n.stages)
  lb <- rep(-20, x$n.stages)
  ivec <- rep(NA, x$n.stages)


  
  exc <- function(u, I, stage,  theta, target){
    gg <- gsDesign::gsProbability(k=stage,
                                  theta=theta,
                                  n.I =c(head(ivec, stage-1),I),
                                  a = c(head(lb, stage-1), -20),
                                  b = c(head(ub, stage-1), u))
     sum(gg$upper$prob[1:stage]) - target
  }
  
  
  uI <- function(I, stage){
    res <- uniroot(exc, interval=c(qnorm(x$sig.level, lower.tail=FALSE), 20),
                   I = I, stage=stage, theta=0, target = alpha.cum[stage],
                   extendInt = "downX")
   res$root
  }
  
  
  exc_low <- function(l, stage, theta, target, I=ivec[stage]){
    gg <- gsDesign::gsProbability(k=stage,
                                  theta=theta,
                                  n.I =c(head(ivec, stage-1), I),
                                  a = c(head(lb, stage-1), l),
                                  b = ub[1:stage])
     sum(gg$lower$prob[1:stage]) - target
  }
  
  
  lI <- function(stage, theta, I=ivec[stage]){
    res <- uniroot(exc_low, interval=c(-20, qnorm(x$sig.level, lower.tail=FALSE)),
                   stage=stage, theta=theta, target = x$power.futility, I=I,
                   extendInt = "upX")
   res$root
  }
  


  ivec[1] <- ztest.I(delta=x$rE.seq[1], sig.level=alpha.seq[1], power=x$power.efficacy)
  ub[1] <- qnorm(alpha.seq[1], lower.tail=FALSE)


  if (x$n.stages > 1){
    for (k in 2:x$n.stages){
    
      
      if (k == x$n.stages){
        power.target <- x$power
      } else {
        power.target <- x$power.efficacy
      }
      .th <- x$rE.seq[k]

      
      
        if (x$futility.type == "binding"){
          lb[k-1] <- lI(stage=k-1, theta=x$rF.seq[k-1])
          
            typeII.overspent <- exc_low(lb[k-1], stage=k-1, theta=.th, target = 1 - power.target)
            if (typeII.overspent > 0){
              exc_low_i <- function(ii)exc_low(lI(I=ii, theta=x$rF.seq[k-1], stage=k-1), ii, stage=k-1, theta=.th,
                                               target = 1-power.target)
              resI_low <- uniroot(exc_low_i, interval=c(ivec[k-1], 2*ivec[k-1]), extendInt = "downX")
              ivec[k-1] <- resI_low$root
              ub[k-1] <- uI(I = resI_low$root, stage = k-1)
              lb[k-1] <- lI(I = resI_low$root, theta = x$rF.seq[k-1], stage = k-1)
            }
          
        } else {
          lb[k-1] <- -20
        }
      
      
      exc_i <- function(ii)exc(uI(ii, k), ii, stage=k, theta=.th,
                                      target = power.target)

      minI <- max(ztest.I(delta = .th, power=power.target, sig.level=alpha.cum[k]),
                  ivec[k-1])
      curr.exc <- exc_i(minI)

      if (curr.exc > 0){
        # already past target power
        ivec[k] <- minI
        ub[k] <- lb[k] <- uI(minI, k)
      } else {
        resI <- uniroot(exc_i, interval=c(minI, 2*minI), extendInt = "upX")
        ivec[k] <- resI$root
        ub[k] <- lb[k] <- uI(resI$root, k)
      }
      
    
    }
  }

  if (x$futility.type == "non-binding"){
    
      for (k in 1:(x$n.stages-1)){
        lb[k] <- lI(stage=k, theta=x$rF.seq[k])
      }
    
  }

  x$upper <- ub
  x$lower <- lb
  x$info <- ivec
  x$spending <- alpha.seq
  return(x)
}

