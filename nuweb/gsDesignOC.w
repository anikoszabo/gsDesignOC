\documentclass[reqno]{amsart}
\usepackage[margin=1in]{geometry}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\renewcommand{\NWtarget}[2]{\hypertarget{#1}{#2}}
\renewcommand{\NWlink}[2]{\hyperlink{#1}{#2}}
\providecommand{\ds}{\displaystyle}
\newcommand{\I}{\mathcal{I}}
\newcommand{\R}{\mathcal{R}}
\newcommand{\A}{\mathcal{A}}
\newcommand{\C}{\mathcal{C}}
\newtheorem{theorem}{Theorem}
\newtheorem{corollary}{Corollary}


\title{Operating-characteristic guided design of group-sequential trials}
    \author{Aniko Szabo}
    \date{\today}

\begin{document}
\begin{abstract}
Group-sequential designs are commonly used for clinical trials to allow early stopping for efficacy or futility. While the design of a single-stage randomized trial is guided by a target power for an alternative hypothesis of interest, the addition of interim analyses is driven by technical choices that are less understandable for clinicians. For example, the commonly used Lan-DeMets methodology requires specification of the timing of analyses and error spending functions. Since the rationale and effect of these technical choices is often unclear, the operating characteristics of the final design are explored under various values of the parameter of interest, and the design is then adjusted until desired properties are obtained.

In this work we develop methods for constructing designs that achieve the desired operating characteristics without the need to specify error spending functions or the timing of analyses. Specifically, we consider designing a study for the mean difference $\theta$ of a normally distributed outcome with known variance. The null hypothesis $H_0: \theta=\theta_0$ is tested versus $H_a: \theta=\theta_A$, with power $\pi$ at a significance level $\alpha$. The interim analyses are designed so that for a pre-specified sequence $\theta_{Ak}$ the study stops for efficacy at stage $k$ with probability $\pi$ if $\theta=\theta_{Ak}$. If stopping for futility is also considered, then the requirement to stop for futility at stage $k$ with probability $\pi_F$ if $\theta=\theta_{0k}$ for pre-specified sequence $\theta_{0k}$  can also be added.
We show that under some monotonicity restrictions, such designs exist for any choice of the timing of interim analyses. Specific designs can be selected by imposing additional optimality requirements, such as minimizing the expected sample size under the target alternative $\theta_A$, or the average sample size under a weighted selection of the alternatives.
\end{abstract}
\maketitle

\section{Setup and notations}

All the calculations are set up assuming a setting with a normally distributed outcome with a known variance, but many common settings can be reduced to / approximated by this special case.

Let $X_i \sim N(\theta, \sigma^2)$ be a sequence of independent identically distributed measurements. We are designing a study to test the following hypotheses:
\begin{equation}
H_0: \theta = \theta_0 \quad\text{versus}\quad H_A: \theta\geq\theta_0
\end{equation}
with type I error $\alpha$ and power $\pi$ at a pre-specified value $\theta=\theta_A$.

The test statistic for $H_0$ based on $n$ values is $Z(n) = (\bar{X}_n - \theta_0)\sqrt{n}\big/\sigma = \sqrt{\I_n}(\bar{X}_n - \theta_0)$, where $\I_n = [\sigma^2/n]^{-1}$ is the Fisher information after $n$ observations. Then for any value of $\theta$
\begin{equation}
Z(n) \sim N\big(\sqrt{\I_n}(\theta-\theta_0), 1\big),
\end{equation}
and for $n_1 \leq n_2$
\begin{equation}
Cov(Z_{n_1}, Z_{n_2}) = \sqrt{\I_{n_1} / \I_{n_2}}.
\end{equation}
This information-based formulation holds (approximately) for a variety of potential designs and outcomes, including two-arm randomized and cross-over studies with normal or binary outcomes.


The study will use a group-sequential design with $K$ analyses, conducted at relative times $t_1,\ldots,t_K=1$. The $k^\text{th}$ analysis will be conducted after $n_k = N \times t_k$ observations have been obtained, where $N$ is the maximum sample size of the study. Based on the test statistic $Z_k = Z(n_k)$ at  analysis $k=1,\ldots,K-1$, the trial can be stopped early for efficacy (reject $H_0$) if $Z_k \geq u_k$, stopped early for futility (fail to reject $H_0$) if $Z_k \leq l_k$, or continued to the next stage. At the final, $K^\text{th}$ analysis,  $H_0$ would either be rejected if $Z_K \geq u_K=l_K$, or we would fail to reject $H_0$ otherwise. As a special case, setting $l_k=-\infty$ will result in a study when interim stopping for futility is not considered.

The decision boundaries $u_k, l_k, k=1,\ldots,K$, and the maximum information target ($\sim$ sample  size) $\I_\text{max}$ have to be selected so that the overall type I error and power of the study are maintained at their design values, but there are many valid choices. In this \emph{operating characteristic driven} approach we propose to select the boundaries based on the cumulative probabilities of stopping for efficacy or futility at each stage.

\begin{table}\caption{Notations}\label{T:notation}
\begin{center}
\begin{tabular}{p{4cm}|c|c}
Characteristic & Notation & Event \\ \hline
\multicolumn{3}{l}{\textbf{No or non-binding futility boundary}}\\
Rejection \emph{at} stage $k$ & $\R_k$ &
  $\{Z_1<u_1,\ldots,Z_{k-1}<u_{k-1}, Z_k \geq u_k \}$   \\
Continuation \emph{at} stage $k$ & $\C_k$ &
  $\{Z_1<u_1,\ldots,Z_{i-1}<u_{k-1}, Z_k < u_k \}$  \\ \hline
\multicolumn{3}{l}{\textbf{Non-binding futility boundary}}\\
Acceptance \emph{at} stage $k$ & $\A_k$ &
  $\{l_1 < Z_1<u_1,\ldots,l_{k-1} < Z_{k-1}<u_{k-1}, Z_i \leq l_k \}$  \\ \hline
\multicolumn{3}{l}{\textbf{Binding futility boundary}}\\
Rejection \emph{at} stage $k$ & $\R^*_k$ &
  $\{l_1 < Z_1<u_1,\ldots,l_{k-1} <Z_{k-1}<u_{k-1}, Z_k \geq u_k \}$   \\
Continuation \emph{at} stage $k$ & $\C^*_k$ &
  $\{l_1 < Z_1<u_1,\ldots,l_{k-1} <Z_{k-1}<u_{k-1}, l_k < Z_k < u_k \}$  \\
Acceptance \emph{at} stage $k$ & $\A_k$ &
  $\{l_1 < Z_1<u_1,\ldots,l_{k-1} < Z_{k-1}<u_{k-1}, Z_i \leq l_k \}$  \\
\end{tabular}
\end{center}
\end{table}

\begin{table}\caption{Target operating characteristics}\label{T:OCdef}
\begin{center}
\begin{tabular}{p{4cm}|c|c|c}
Characteristic & Event & $\theta$ & Probability target\\ \hline
\multicolumn{4}{l}{\textbf{No or non-binding futility boundary}}\\
Overall type I error &
  $\ds\bigcup_{i=1}^K \R_i$ & $\theta_0$ &
  $\leq\alpha$\\
Overall power &
  $\ds\bigcup_{i=1}^K\ \R_i$ & $\theta_A$ &
  $\geq\pi$\\
Stop by stage $k$ for efficacy  &
  $\ds\bigcup_{i=1}^k \R_i$ & $\theta_{Ak}$ &
  $\geq\pi_E$\\ \hline
\multicolumn{4}{l}{\textbf{Non-binding futility boundary}}\\
Stop by stage $k$ for futility  &
  $\ds\bigcup_{i=1}^k \A_i$ & $\theta_{0k}$ &
  $\geq\pi_F$\\ \hline
\multicolumn{4}{l}{\textbf{Binding futility boundary}}\\
Overall type I error &
  $\ds\bigcup_{i=1}^K \R^*_i$ & $\theta_0$ &
  $\leq\alpha$\\
Overall power &
  $\ds\bigcup_{i=1}^K \R^*_i$ & $\theta_A$ &
  $\geq\pi$\\
Stop by stage $k$ for efficacy  &
  $\ds\bigcup_{i=1}^k \R^*_i$ & $\theta_{Ak}$ &
  $\geq\pi_E$\\
Stop by stage $k$ for futility  &
  $\ds\bigcup_{i=1}^k \A_i$ & $\theta_{0k}$ &
  $\geq\pi_F$\\
\end{tabular}
\end{center}
\end{table}

\begin{theorem}\label{Th:exists} A design satisfying the operating characteristic requirements of Table~\ref{T:OCdef} exists, if the following conditions are satisfied:
\begin{itemize}
\item[C1.] $\pi_E \leq \pi$
\item[C2.] $\theta_{A1} \geq \theta_{A2} \geq \cdots \geq \theta_{A,K-1} \geq \theta_A$
\end{itemize}
Additionally, if a futility boundary is needed:
\begin{itemize}
\item[C3.]  $\pi_F \leq 1-\alpha$
\item[C4.] $\theta_{01} \leq \theta_{02} \leq \cdots \leq \theta_{0,K-1} \leq \theta_0$
\end{itemize}
\end{theorem}
\begin{proof} We will consier the cases of no, non-binding, and binding futility boundaries separately.

\textbf{I. no futility boundary} We will use proof by induction over the number of stages $K$. When $K=1$, only the overall type I error and power need to be considered, and the well-known single-stage design achieving the overall information $$\I_\text{fix} = \frac{(z_\alpha + z_{1-\beta})^2}{(\theta_A - \theta_0)^2}$$ will satisfy the requirements with $u_K = z_\alpha$.

Now suppose we can construct the desired design for $K-1$ stages, and we try to add an additional stage. For the first $K-1$ stages select the boundary $u_1,\ldots,u_{K-1}$ and information times $\I_{n_1},\ldots,\I_{n_{K-1}}$ to achieve over these $K-1$ stages stage-specific stopping probabilities $\pi_E$, overall power $\pi_E$, and overall type I error $\alpha_{K-1} = \alpha-\Delta\alpha_K$, where $0 < \Delta\alpha_K <\alpha$ is an arbitrary value. With these choices, the requirements for all the stage-specific probabilities for the $K$-stage design are satisfied. We need to select $I_{n_K}$ and $u_K$ to satistfy the overall power and type I error restrictions.

Consider the function $a(u \mid \I_N) = Pr\big(\{\bigcup_{i=1}^{K-1}\R_i \}\bigcup \{Z_1 < u_1,\ldots, Z_{K-1} < u_{K-1}, Z(N) \geq u \} \mid \theta_{0}\big)$, that gives the type I error if the boundary is set at $u$ for the last stage for any given $I_N > I_{n_{K-1}}$. Under $H_0$, $Z(N) ~ \sim N(0,1)$, so $a(z_\alpha\mid\I_N) \geq \alpha$.  On the other hand, $a(\infty | \I_N) = \alpha - \Delta\alpha_K < \alpha$ by the design of the first $K-1$ stages. Since $a$ is monotone in $u$, for any $I_N$ there exists a cutoff $u_K(\I_N)$ for which $a(u_K(\I_N)) = \alpha$, ie for which the overall type I error is maintained at the desired level.

Next consider the function $b(\I_N) = Pr\big(\{\bigcup_{i=1}^{K-1}\R_i\} \bigcup \{Z_1 < u_1, \ldots, Z_{K-1}< u_{K-1}, Z(N) \geq u_K(\I_N)  \mid \theta_{A}\big)$ that gives the power to reject $H_0$ if the final analysis is set at information $\I_N > I_{n_{K-1}}$. If $b(\I_{K-1}) \geq \pi$, then we can set $\I_N = \I_{K-1}$.
%Since $\theta_A \leq \theta_{A,K-1}$, $b(\I_{K-1}) \leq \pi_E \leq \pi$.
Otherwise, note that $b(\infty) = 1$, and $b$ is monotone in $\I_N$. Thus we can select a value $n_K$ such that $b(\I_{n_K}) = \pi$, which will provide the desired design.

\textbf{II. Non-binding futility boundary}
The overall power, type I error, and stage-specific efficacy power requirements do not depend on a non-binding futility boundary, so we can start with a design constructed without a futility boundary. We then need to calculate the lower bounds $l_k$ to satisfy the futility stopping probabilities. For stage 1, we need $Pr(Z_1 < l_1 | \theta_{01}) = \pi_F$, where $Z_1 \sim N(\sqrt{I_1}(\theta_{01}-\theta_0), 1)$. Thus $l_1 = z_{1-\pi_F} + \sqrt{I}(\theta_{01}-\theta_0)$ satisfies the target power. Given conditions C3-C4, $l_1 \leq z_\alpha$, while $u_1 \geq z_\alpha$, so $l_1 \leq u_1$.

At stage $k$, consider the function $c(l) = Pr(l_1 < Z_1< u_1, \ldots, l_{k-1} < Z_{k-1}< u_{k-1}, Z_k < l \mid \theta_{0k})$. Since $\theta_{0k} \leq \theta_0$, $c(z_\alpha) \geq 1-\alpha \geq \pi_F$, while $c(-\infty) = 0$, so a value $l_k$ can be selected for which $c(l_k) = \pi_F$.

\textbf{III. Binding futility boundary}
The construction follows the induction-based approach of part I. During the inductive step in addition to $I_K$ and $u_K$, we need to find a value for $l_{K-1}$ that will satisfy the futility boundary restriction. This can be done as shown in part II before proceeding with the selection of the parameters of the last stage. A potential complication compared to the efficacy-only situation, is that by adding $l_{K-1}$ we might overspend the type II error, $Pr(\bigcup_{i=1}^{K-1}\A_i | \theta_A) > 1-\pi$. In this case we would not be able to select a value for $\I_N$ such that $b^*(\I_N)\geq \pi$, because even $b^*(\infty) < \pi$. One of the possible solutions is to adjust the timing of the $(K-1)st$ stage to reduce its beta-spending. As $\I_{K-1}$ is increased beyond the value that provided power $\pi_E$ at significance level $\alpha-\Delta\alpha_K$, if we maintain the significance level, $P(\R^*_{K-1}|\theta_{A,K-1})$ will increase, so the power will be maintained. However $Pr(\A_{K-1}|\theta_A)$ will decrease, so we can find a value of $I_{K-1}$ for which $Pr(\bigcup_{i=1}^{K-1}\A_i | \theta_A) = 1-\pi$. Starting from this modified design with $K-1$ stages, we can now construct the required $K$-stage design.

\end{proof}

In the proof of the theorem we have obtained the following uniqueness result:

\begin{corollary}\label{Th:construct} Under the assumptions of Theorem~\ref{Th:exists},  an alpha-spending sequence $\Delta\alpha_k > 0$, $k=1,\ldots,K$ such that $\sum_{k=1}^K\Delta\alpha_k = \alpha$ uniquely determines a design with the required operating characteristics.
\end{corollary}

\clearpage

\section{Main function: calculate design}

@O ../R/gsDesignOC.R
@{
#'Find optimal group-sequential design
#'
#'The \code{gsDesignOC} function finds the analysis times, boundaries, and sample size
#'for a group-sequential trial
#'
#'@@export
#'@@param thA numeric; effect size under the alternative hypothesis
#'@@param thA.seq monotone numeric vector of values getting closer to thA from outside the
#'th0-thA range (ie, if thA > th0, they should form a decreasing sequence with all values > thA). #'For the k-th value the study will stop for efficacy at or before the k-th stage with probability \code{power.efficacy}.
#'@@param th0 numeric; effect size under the null hypothesis
#'@@param th0.seq monotone numeric vector of values getting closer to th0 from outside the
#'th0-thA range (ie, if thA > th0, they should form an increasing sequence with all values < th0).
#'The study will stop for futility at or before the k-th stage with probability \code{power.futility}.
#'Its length should be equal to that of \code{thA.seq}. The default value of \code{NULL} implies no
#'futility monitoring.
#'@@param sig.level numeric; the one-sided significance level. The default is 0.025.
#'@@param power numeric; the desired study power. The default is 0.9.
#'@@param power.efficacy numeric; the probability of stopping for efficacy at stage k under \code{thA.seq}
#'@@param power.futility numeric; the probability of stopping for futility at stage k under \code{th0.seq}
#'@@param futility.type character string, one of \code{c("none", "non-binding","binding")} or their
#'unique abbreviations. Specifies whether a futility boundary should not be used at all ("none"), or if it
#'is used, then whether the effect of the futility boundary should be included
#'in the efficacy power/type I error calculations ("binding"), or not ("non-binding").
#'@@param mix.theta numeric vector; specifies the set of alternatives under which the design is optimized.
#' Defaults to \code{thA}.
#'@@param mix.w numeric vector of length equal to that of \code{mix.theta}. The weights of the
#'elements of \code{mix.theta} in the optimality criterion. It will be normalized to a sum of 1.
#'Defaults to equal weights.
#'@@param control list of parameters controlling the estimation altorithm to replace the default
#'values returned by the \code{glsControl} function.
#'@@return an object of class \code{gsDesignOC}
#'@@author Aniko Szabo
#'@@references Szabo, A, Tarima, S (?) Operating-characteristic guided design of group-sequential trials.
#'@@keywords nonparametric design
#'@@examples
#'
#'gsDesignOC(0.3, thA.seq = c(1, 0.5))
#'
#'@@name gsDesignOC

gsDesignOC <- function(thA, thA.seq, th0=0, th0.seq=NULL,
                       sig.level = 0.025,
                       power=0.9, power.efficacy=power, power.futility = power,
                       futility.type=c("none","non-binding","binding"),
                       mix.theta = thA,
                       mix.w = rep(1, length(mix.theta)),
                       control=list()){
  @< Check inputs @>

  # create skeleton gsDesignOC object
  res <- list(thA=thA, thA.seq=thA.seq,
              th0=th0, th0.seq=th0.seq,
              sig.level = sig.level,
              power=power, power.efficacy = power.efficacy, power.futility = power.futility,
              futility.type=futility.type)
  class(res) <- "gsDesignOC"

  @< Define optimization function @>

  n.stages <- length(thA.seq) + 1
  if (n.stages == 1){
    alpha.seq <- sig.level
  } else if (n.stages == 2){
    oo <- optimize(.cp, interval=c(-5,5))

    y.res <- oo$minimum
    alpha.seq <- sig.level * exp(c(y.res,0)) / sum(exp(c(y.res,0)))
  } else {
    y.start <- -log(seq(n.stages, 2, by=-1))

    oo <- optim(y.start, fn=.cp, method="Nelder-Mead")

    y.res <- oo$par
    alpha.seq <- sig.level * exp(c(y.res,0)) / sum(exp(c(y.res,0)))
  }
  res <- calc.bounds(x=res, alpha.seq)

  return(res)
}

@}

@d Check inputs @{
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

@}


Using the \texttt{calc.bounds} function, we can construct a design matching all the type I error/power requirements given a sequence of stage-specific positive alpha-spending values $\Delta\alpha_k$ that add up to $\alpha$. To allow unconstrained optimization, we can reparametrize it using the multinomial logit:
$$ y_k = \log\frac{\Delta\alpha_k}{\Delta\alpha_K}, \quad k=1,\ldots,K-1,$$
with reverse transition
$$ \Delta\alpha_k = \frac{\alpha e^{y_k}}{\sum_{i=1}^K e^{y_i}}, \quad k=1,\ldots, K,$$
where $y_K := 0$.

@D Define optimization function @{
  .cp <- function(y.vec){

    y <- c(y.vec, 0)  # add y_K
    alpha.seq <- sig.level * exp(y) / sum(exp(y))

    res <- calc.bounds(x=res, alpha.seq = alpha.seq)
    dp <- oc(res, EN_theta=mix.theta, mix.w=mix.w)
    Q <- dp$ave.EN
    Q
  }
@}

\section{Key support functions}



@O ../R/gsDesignOC.R
@{
#'Operating characteristics of a potential design
#'
#'Calculates the average expected information under a weighted set of values for the alternative hypothesis, and
#' the probability of stopping for futility or effecacy under the design alternatives. If the futility boundary
#' is non-binding, then the lower boundary is not included in the calculation of the expected sample size and
#' efficacy stopping probabilities, but it is included in the futility stopping calculations.
#'
#'@@export
#'@@param x object of class \code{gsDesignOC}
#'@@param EN_theta numeric vector of parameter values for which expected sample size calculation is performed
#'@@param mix.w numeric vector of positive weights for each value of \code{EN_theta}. Defaults to equal weights.
#'@@return a list with the following elements:
#'\describe{
#'\item{ave.EN}{numeric; weighted average of expected sample sizes under the requested alternatives}
#'\item{EN.vec}{numeric vector; expected sample sizes under each of the requested alternatives}
#'\item{thA.cumcross}{numeric vector; cumulative probability of stopping for efficacy at or before stage k
#'  under thA.seq[k], including thA at the end.}
#'\item{th0.cumcross}{numeric vector}{numeric vector; cumulative probability of stopping for futility at or before stage k
#'  under th0.seq[k], including th0 at the end.}
#'}
#'@@author Aniko Szabo
#'@@keywords design
#'@@examples
#'
#'

oc <- function(x, EN_theta=x$thA,  mix.w = rep(1, length(EN_theta))){

  if (length(mix.w) != length(EN_theta))
    stop("`EN_theta` and `mix.w` should have the same length")
  if (!all(mix.w > 0))
    stop("`mix.w` should have only positive elements")

  n.stages <- length(x$thA.seq) + 1
  n.EN <- length(EN_theta)
  n.A <- length(c(x$thA.seq, x$thA))
  n.0 <- length(c(x$th0.seq, x$th0))

  # crossing probabilities for futility
  ggF <- gsDesign::gsProbability(k = n.stages,
                                  theta = c(EN_theta, x$thA.seq, x$thA, x$th0.seq, x$th0),
                                  n.I = x$info,
                                  a = x$lower,
                                  b = x$upper)

  # crossing probabilities for efficacy and EN
  if (x$futility.type == "non-binding") {
    #ignore lower boundary (except last stage) for EN & efficacy stopping
    ggE <- gsDesign::gsProbability(k = n.stages,
                                  theta = c(EN_theta, x$thA.seq, x$thA, x$th0.seq, x$th0),
                                  n.I = x$info,
                                  a = c(rep(-20, n.stages-1), x$lower[n.stages]),
                                  b = x$upper)
   } else {
    ggE <- ggF
   }

  # expected sample size calculations
  en_vec <- ggE$en[1:n.EN]

  en <- c(en_vec %*% mix.w / sum(mix.w))

  # cumulative stopping probability calculations
  p.up <- ggE$upper$prob[, n.EN + (1:n.A)]
  cump.up <- apply(p.up, 2, cumsum)
  thA.cumcross <- diag(cump.up)

  p.low <- ggF$lower$prob[, n.EN + n.A + (1:n.0)]
  cump.low <- apply(p.low, 2, cumsum)
  th0.cumcross <- diag(cump.low)

  list(ave.EN = en,
       EN.vec = en_vec,
       thA.cumcross = thA.cumcross,
       th0.cumcross = th0.cumcross)
}
@| oc@}

@o ../R/gsDesignOC.R
@{
#'Control values for gsDesignOC
#'
#'The values supplied in the function call replace the defaults and a list with all possible
#'arguments is returned. The returned list is used as the \code{control} argument to the
#'\code{gsDesignOC} function.
#'
#'@@export
#'@@param optim.penalty numeric; relative weight of terms ensuring target type I and type II
#'errors versus sample size in optimization routine
#'@@return a list with components for each of the possible arguments
#'@@author Aniko Szabo
#'@@keywords design
#'@@examples
#'
#'ocControl(optim.penalty = 100)

ocControl <- function(optim.penalty = 1000){
  list(optim.penalty = optim.penalty)
}
@}

\section{Boundary achieving the desired operating characteristics}

Given the alpha-spending sequence and using the probability targets defined in Table~\ref{T:OCdef}, we can derive the information times $\I_k$, and boundaries $u_k$ and $l_k$ using the recursive construction process described in the proof of Theorem~\ref{Th:exists}.



@O ../R/gsDesignOC.R
@{
#'Calculate the information times and efficacy/futility boundary values given the alpha-spending sequence and desired operating characteristics

#'@@export
#'@@param x an object of class \code{gsDesignOC}

calc.bounds <- function(x, alpha.seq){

  n.stages <- length(x$thA.seq) + 1
  alpha.cum <- cumsum(alpha.seq)
  ub <- rep(20, n.stages)
  lb <- rep(-20, n.stages)
  ivec <- rep(NA, n.stages)

@< Define exceedance probability functions  @>
@< Construct first stage @>

  if (n.stages > 1){
    for (k in 2:n.stages){
    @< Construct next stage @>
    }
  }

  if (x$futility.type == "non-binding"){
    @< Add non-binding lower bound @>
  }

  x$upper <- ub
  x$lower <- lb
  x$info <- ivec
  return(x)
}

@|calc.bounds @}

The first stage is based on the fixed sample-size formula for alternative $\theta_{A1}$, power $\pi_E$, and significance level $\Delta\alpha_1$.
@D Construct first stage @{
  ivec[1] <- ztest.I(delta=x$thA.seq[1]-x$th0, sig.level=alpha.seq[1], power=x$power.efficacy)
  ub[1] <- qnorm(alpha.seq[1], lower=FALSE)
@}

The cumulative rejection / acceptance probabilities depend on the type of the futility boundary.
@D Define exceedance probability functions  @{
  @< Define upper bound exceedance @>
  @< Define function to find u(I) @>
  @< Define lower bound exceedance @>
  @< Define function to find l @>
@}

\subsection{Upper bound exceedance}

The functions $a$ and $b$ in the proof of Theorem~\ref{Th:exists} have a similar form, so we will define the following function:
$$e(u \mid \I) = Pr\big(
  \bigcup_{i=1}^{k-1}\{l_1< Z_1<u_1,\ldots, l_{i-1}< Z_{i-1} < u_{i-1}, Z_i \geq u_i \}
 \bigcup \{l_1< Z_1<u_1,\ldots, l_{k-1}< Z_{k-1} < u_{k-1}, Z(\I) \geq u \} \mid \theta\big).$$

If $l_i$ is set to $-\infty$, then that is equivalent to not having a lower boundary.
We will need to solve equations of the form $e(u)=c$, so we will actually define a function for the difference of the LHS and RHS.

@d Define upper bound exceedance  @{
exc <- function(u, I, stage,  theta, target){
  gg <- gsDesign::gsProbability(k=stage,
                                theta=theta,
                                n.I =c(head(ivec, stage-1),I),
                                a = c(head(lb, stage-1), -20),
                                b = c(head(ub, stage-1), u))
   sum(gg$upper$prob[1:stage]) - target
}
@}

To find $u(I)$ we need to solve $a(u|\I) = \sum_{i=1}^k \Delta\alpha_i=\tilde{\alpha}_k$. Since $a(z_\alpha\mid\I_N) \geq \alpha$ and $a(\infty | \I_N) = \tilde{\alpha}_{k-1}$, the solution for fixed $\I_N$ is between $z_\alpha$ and $20$ (which is essentially $\infty$).

@d Define function to find u(I) @{
uI <- function(I, stage){
  res <- uniroot(exc, interval=c(qnorm(x$sig.level, lower=FALSE), 20),
                 I = I, stage=stage, theta=x$th0, target = alpha.cum[stage],
                 extendInt = "downX")
 res$root
}
@}



\subsection{Construct next stage}

@d Construct next stage @{
  @< Determine target power and alternatives @>
  @< Find lower bound for previous stage @>
  @< Find I for next stage @>
@}


To find $I_k$, we need to solve $b(\I) = \pi_E$, if $k<K$ and $=\pi$ for $k=K$.
@d Determine target power and alternatives @{
if (k == n.stages){
  power.target <- x$power
  .th <- x$thA
} else {
  power.target <- x$power.efficacy
  .th <- x$thA.seq[k]
}
@}


Note that $b(I_{k-1}) \leq \pi_E$ and $b(\I_{fix}(\theta_{Ak}-\theta_0, \tilde{\alpha}_k, \pi_E)) \leq \pi_E$, we can start the search at the maximum of these two values.

@D Find I for next stage @{
exc_i <- function(ii)exc(uI(ii, k), ii, stage=k, theta=.th,
                                target = power.target)

minI <- max(ztest.I(delta = .th - x$th0, power=power.target, sig.level=alpha.cum[k]),
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
@}

\subsection{Add lower bound}

The lower bound will usually be computed only when the information and upper bound for a stage has been determined.

@d Define lower bound exceedance  @{
exc_low <- function(l, stage, theta, target, I=ivec[stage]){
  gg <- gsDesign::gsProbability(k=stage,
                                theta=theta,
                                n.I =c(head(ivec, stage-1), I),
                                a = c(head(lb, stage-1), l),
                                b = ub[1:stage])
   sum(gg$lower$prob[1:stage]) - target
}
@}


@d Define function to find l @{
lI <- function(stage, theta, I=ivec[stage]){
  res <- uniroot(exc_low, interval=c(-20, qnorm(x$sig.level, lower=FALSE)),
                 stage=stage, theta=theta, target = x$power.futility, I=I,
                 extendInt = "upX")
 res$root
}
@}

We need to find the lower bound for stage $k-1$, before computing the size of the next stage. If the boundary is not binding, we need to overwrite $l_k:= u_k$ that was indicating the "last" stage so far to $l_k=-\infty$.
@D Find lower bound for previous stage @{
  if (x$futility.type == "binding"){
    lb[k-1] <- lI(stage=k-1, theta=x$th0.seq[k-1])
    @< Check beta-spending and adjust previous stage if needed @>
  } else {
    lb[k-1] <- -20
  }
@}

With a binding boundary it is possible that the addition of $l_{K-1}$ overspends the type II error under $\theta_{Ak}$, and the target power for the next stage is not achievable anymore. In this case, the information for the $(K-1)$st stage needs to be increased.
@D Check beta-spending and adjust previous stage if needed @{
  typeII.overspent <- exc_low(lb[k-1], stage=k-1, theta=.th, target = 1 - power.target)
  if (typeII.overspent > 0){
    exc_low_i <- function(ii)exc_low(lI(I=ii, theta=x$th0.seq[k-1], stage=k-1), ii, stage=k-1, theta=.th,
                                     target = 1-power.target)
    resI_low <- uniroot(exc_low_i, interval=c(ivec[k-1], 2*ivec[k-1]), extendInt = "downX")
    ivec[k-1] <- resI_low$root
    ub[k-1] <- uI(I = resI_low$root, stage = k-1)
    lb[k-1] <- lI(I = resI_low$root, theta = x$th0.seq[k-1], stage = k-1)
  }
@}


For a non-binding boundary, the lower bound is added only after the entire upper bound has been calculated.
@D Add non-binding lower bound @{
  for (k in 1:(n.stages-1)){
    lb[k] <- lI(stage=k, theta=x$th0.seq[k])
  }
@}


\section{Utility help-functions}



@o ../R/Utility.R
@{
#'Fixed sample-size information for one-sided test
#'
#'The \code{ztest.I} function finds the information required of a one-sided Z-test with given power
#'and significance level
#'
#'@@param delta numeric; targeted standardized effect size
#'@@param sig.level numeric; one-sided significance level
#'@@param power numeric; target power
#'@@return numeric value with the (fractional) information achieving the target power
#'@@author Aniko Szabo
#'@@keywords internal
#'@@examples
#'
#'ztest.I(delta=1, sig.level=0.05, power=0.9)

ztest.I <- function(delta, sig.level, power){
  za <- qnorm(sig.level, lower=FALSE)
  zb <- qnorm(power)
  I <- (za+zb)^2 /delta^2
  I
}
@| ztest.I @}

@O ../R/Utility.R
@{
#'Conversion between boundary and nominal significance level
#'
#'The \code{nom.to.bnd} and \code{bnd.to.nom} functions perform conversion between the boundary and
#' matching significance level.
#'
#'@@param nominal.level numeric vector or matrix with two rows of the nominal significance level for each stage.
#'  When a matrix, row 1 corresponds to the efficacy boundary and row 2 to the futility. Defaults to NULL.
#'@@param cutoffs numeric vector or numeric matrix with two rows of the z-test cutoffs for each stage.
#'  When a matrix, row 1 corresponds to the efficacy boundary and row 2 to the futility. Defaults to NULL.
#'@@param n integer vector of sample sizes of each stage. Its sum is the total study sample size.
#'@@param theta.null numeric, the null hypothesis being tested
#'@@param theta.alt numeric, the alternative hypothesis being tested for the futility bound calculation
#'@@keywords internal
#'@@examples
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

#'@@describeIn nom.to.bnd Conversion from nominal significance level to boundary
#'@@export
#'@@keywords internal
#'@@examples
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

@}




\end{document}
