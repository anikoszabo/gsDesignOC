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
  $\{Z_1<u_1,\ldots,Z_{k-1}<u_{k-1}, Z_k < u_k \}$  \\ \hline
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
  $\geq\pi_{Ek}$\\ \hline
\multicolumn{4}{l}{\textbf{Non-binding futility boundary}}\\
Stop by stage $k$ for futility  &
  $\ds\bigcup_{i=1}^k \A_i$ & $\theta_{0k}$ &
  $\geq\pi_{Fk}$\\ \hline
\multicolumn{4}{l}{\textbf{Binding futility boundary}}\\
Overall type I error &
  $\ds\bigcup_{i=1}^K \R^*_i$ & $\theta_0$ &
  $\leq\alpha$\\
Overall power &
  $\ds\bigcup_{i=1}^K \R^*_i$ & $\theta_A$ &
  $\geq\pi$\\
Stop by stage $k$ for efficacy  &
  $\ds\bigcup_{i=1}^k \R^*_i$ & $\theta_{Ak}$ &
  $\geq\pi_{Ek}$\\
Stop by stage $k$ for futility  &
  $\ds\bigcup_{i=1}^k \A_i$ & $\theta_{0k}$ &
  $\geq\pi_{Fk}$\\
\end{tabular}
\end{center}
\end{table}

%% TODO: switch C1-C2 and C3-C4
%% Some restriction $(\pi_{Ek} - \pi_{E,k-1})(\pi_{Ak} - \pi_{A,k-1}) < 0$ is probably needed to ensure
%% that at least one of them changes

\begin{theorem}\label{Th:exists} A design satisfying the operating characteristic requirements of Table~\ref{T:OCdef} exists, if the following conditions are satisfied:
\begin{itemize}
\item[C1.] $\pi_{E1}\leq \pi_{E2} \leq\cdots\leq \pi_{E,K-1} \leq \pi$
\item[C2.] $\theta_{A1} \geq \theta_{A2} \geq \cdots \geq \theta_{A,K-1} \geq \theta_A$
\end{itemize}
Additionally, if a futility boundary is needed:
\begin{itemize}
\item[C3.]  $\pi_{F1} \leq\pi_{F2} \leq\cdots\leq \pi_{F,K-1} \leq 1-\alpha$
\item[C4.] $\theta_{01} \leq \theta_{02} \leq \cdots \leq \theta_{0,K-1} \leq \theta_0$
\end{itemize}
\end{theorem}
\begin{proof} We will consider the cases of no, non-binding, and binding futility boundaries separately.

\textbf{I. no futility boundary} We will use proof by induction over the number of stages $K$. When $K=1$, only the overall type I error and power need to be considered, and the well-known single-stage design achieving the overall information \begin{equation}\label{E:Ifix}
\I_\text{fix} = \frac{(z_\alpha + z_{1-\beta})^2}{(\theta_A - \theta_0)^2}
\end{equation}
will satisfy the requirements with $u_K = z_\alpha$.

Now suppose we can construct the desired design for $K-1$ stages, and we try to add an additional stage. For the first $K-1$ stages select the boundary $u_1,\ldots,u_{K-1}$ and information times $\I_{n_1},\ldots,\I_{n_{K-1}}$ to achieve over these $K-1$ stages stage-specific stopping probabilities $\pi_{E,K-2}$ at $\theta_{A1},\ldots,\theta_{A,K-2}$, overall power $\pi_{E,K-1}$ at $\theta_{A,K-1}$, and overall type I error $\alpha_{K-1} = \alpha-\Delta\alpha_K$, where $0 < \Delta\alpha_K <\alpha$ is an arbitrary value. With these choices, the requirements for all the stage-specific probabilities for the $K$-stage design are satisfied. We need to select $I_{n_K}$ and $u_K$ to satisfy the overall power and type I error restrictions.

Consider the function $a(u \mid \I_N) = Pr\big(\{\bigcup_{i=1}^{K-1}\R_i \}\bigcup \{Z_1 < u_1,\ldots, Z_{K-1} < u_{K-1}, Z(N) \geq u \} \mid \theta_{0}\big)$, that gives the type I error if the boundary is set at $u$ for the last stage for any given $I_N > I_{n_{K-1}}$. Under $H_0$, $Z(N) ~ \sim N(0,1)$, so $a(z_\alpha\mid\I_N) \geq \alpha$.  On the other hand, $a(\infty | \I_N) = \alpha - \Delta\alpha_K < \alpha$ by the design of the first $K-1$ stages. Since $a$ is monotone in $u$, for any $I_N$ there exists a cutoff $u_K(\I_N)$ for which $a(u_K(\I_N)) = \alpha$, ie for which the overall type I error is maintained at the desired level.

Next consider the function $b(\I_N) = Pr\big(\{\bigcup_{i=1}^{K-1}\R_i\} \bigcup \{Z_1 < u_1, \ldots, Z_{K-1}< u_{K-1}, Z(N) \geq u_K(\I_N)  \mid \theta_{A}\big)$ that gives the power to reject $H_0$ if the final analysis is set at information $\I_N > I_{n_{K-1}}$ and the cutoff $u_K(\I_N)$ defined above is used. If $b(\I_{K-1}) \geq \pi$, then we can set $\I_N = \I_{K-1}$.
%Since $\theta_A \leq \theta_{A,K-1}$, $b(\I_{K-1}) \leq \pi_E \leq \pi$.
Otherwise, note that $b(\infty) = 1$, and $b$ is monotone in $\I_N$. Thus we can select a value $n_K$ such that $b(\I_{n_K}) = \pi$, which will provide the desired design.

\textbf{II. Non-binding futility boundary}
The overall power, type I error, and stage-specific efficacy power requirements do not depend on a non-binding futility boundary, so we can start with a design constructed without a futility boundary. We then need to calculate the lower bounds $l_k$ to satisfy the futility stopping probabilities. For stage 1, we need $Pr(Z_1 < l_1 | \theta_{01}) = \pi_{F1}$, where $Z_1 \sim N(\sqrt{I_1}(\theta_{01}-\theta_0), 1)$. Thus $l_1 = z_{1-\pi_{F1}} + \sqrt{I_1}(\theta_{01}-\theta_0)$ satisfies the target power. Given conditions C3--C4, $l_1 \leq z_\alpha$, while $u_1 \geq z_\alpha$, so $l_1 \leq u_1$.

At stage $k$, consider the function $c(l) = Pr(l_1 < Z_1< u_1, \ldots, l_{k-1} < Z_{k-1}< u_{k-1}, Z_k < l \mid \theta_{0k})$. Since $\theta_{0k} \leq \theta_0$, $c(z_\alpha) \geq 1-\alpha \geq \pi_{Fk}$, while $c(-\infty) = 0$, so a value $l_k$ can be selected for which $c(l_k) = \pi_{Fk}$.

\textbf{III. Binding futility boundary}
The construction follows the induction-based approach of part I. During the inductive step in addition to $I_K$ and $u_K$, we need to find a value for $l_{K-1}$ that will satisfy the futility boundary restriction. This can be done as shown in part II before proceeding with the selection of the parameters of the last stage. A potential complication compared to the efficacy-only situation is that by adding $l_{K-1}$ we might overspend the type II error, $Pr(\bigcup_{i=1}^{K-1}\A_i | \theta_A) > 1-\pi$. In this case we would not be able to select a value for $\I_N$ such that $b^*(\I_N)\geq \pi$, because even $b^*(\infty) < \pi$. One of the possible solutions is to adjust the timing of the $(K-1)st$ stage to reduce its beta-spending. As $\I_{K-1}$ is increased beyond the value that provided power $\pi_E$ at significance level $\alpha-\Delta\alpha_K$, if we maintain the significance level, $P(\R^*_{K-1}|\theta_{A,K-1})$ will increase, so the power will be maintained. However $Pr(\A_{K-1}|\theta_A)$ will decrease, so we can find a value of $I_{K-1}$ for which $Pr(\bigcup_{i=1}^{K-1}\A_i | \theta_A) = 1-\pi$. Starting from this modified design with $K-1$ stages, we can now construct the required $K$-stage design.

\end{proof}

In the proof of the theorem we have obtained the following uniqueness result:

\begin{corollary}\label{Th:construct} Under the assumptions of Theorem~\ref{Th:exists},  an alpha-spending sequence $\Delta\alpha_k > 0$, $k=1,\ldots,K$ such that $\sum_{k=1}^K\Delta\alpha_k = \alpha$ uniquely determines a design with the required operating characteristics.
\end{corollary}

\clearpage

\section{Main function: calculate design}

We want to set up the design as a function of \texttt{n.fix}, the sample size required for the single-stage study, and report the results in terms of sample size. Since $\I \propto n$, we can do the calculations in terms of information, but then rescale everything as
\begin{equation}\label{E:n_I_conversion}
 n = \frac{\I}{\I_{fix}} n_{fix}.
\end{equation}

The second modification from the theoretical development is that instead of standardized effect sizes we will use the relative size of the interim effects, scaling everything to $\theta_A - \theta_0$:
\begin{gather*}
r^E_k = \frac{\theta_{Ak} - \theta_{0}}{\theta_A - \theta_0},\\
r^F_k = \frac{\theta_{0k} - \theta_{0}}{\theta_A - \theta_0}.
\end{gather*}

Keeping the fold-changes $r_k$ constant, a $C$-fold decrease in the effect size $\theta_A - \theta_0$ is equivalent to a $C$-fold increase in the required sample size/information. So we can build the design assuming $\theta_0 = 0$ and $\theta_A=1$, and then adjust the final sample size as
\begin{equation}\label{E:n_theta_conversion}
n(\theta_0, \theta_1) = \frac{n(0,1)}{(\theta_A - \theta_0)^2}.
\end{equation}

Putting together equations \ref{E:Ifix}, \ref{E:n_I_conversion}, and \ref{E:n_theta_conversion}, we have
\begin{equation}\label{E:n_comb}
n(\theta_0, \theta_1) = \frac{\I(0,1)}{(z_\alpha+z_{1-\beta})^2}n_{fix},
\end{equation}
where $\I(0,1)$ is computed assuming $\theta_0 = 0$ and $\theta_A=1$.
@O ../R/gsDesignOC.R
@{
#'Find optimal group-sequential design
#'
#'The \code{gsDesignOC} function finds the analysis times, boundaries, and sample size
#'for a group-sequential trial
#'
#'@@export
#'@@param n.stages integer; number of stages planned. One stage implies no interim analysis.
#'@@param rE.seq monotone non-increasing numeric vector of length \code{n.stages} ending with 1.
#' If it has length \code{n.stages-1}, a 1 will be appended to the end.
#' For the k-th interim value the study will stop for efficacy at or before the k-th stage with probability
#' \code{power.efficacy[k]} if the true effect is \code{rE.seq[k]} times the effect under the alternative hypothesis.
#'@@param rF.seq monotone non-decreasing numeric vector of length \code{n.stages} ending with 0.
#' If it has length \code{n.stages-1}, a 0 will be appended to the end.
#' The study will stop for futility at or before the k-th stage with probability \code{power.futility[k]}.
#' if the true effect is \code{rE.seq[k]} times the effect under the alternative hypothesis, ie
#' the true effect is actually in the opposite direction of the hypothesized effect.
#' The default value of \code{NULL} implies no futility monitoring.
#'@@param n.fix numeric; the sample size for a one-stage design without interim tests. Defaults to 1,
#' in which case the resulting sample sizes can be interpreted as being relative to the single-stage study sample size.
#'@@param sig.level numeric; the one-sided significance level. The default is 0.025.
#'@@param power numeric; the desired study power. The default is 0.9.
#'@@param power.efficacy numeric value or non-decreasing vector of length \code{n.stages} with #'values between `sig.level` and `power`, ending in `power`. If it has length \code{n.stages-1}, #'the value of `power` will be appended to the end. The k-th element is the desired probability #'of stopping for efficacy at stage k under \code{rE.seq[k]}
#'@@param power.futility numeric value or non-decreasing vector  of length \code{n.stages} with values  between `0` and `1-sig.level`, ending with `1-sig.level`.
#'If it has length \code{n.stages-1}, the value of `1-sig.level` will be appended to the end. #'The k-th element is the desired probability of stopping for futility at stage k under \code{rF.seq[k]}
#'@@param futility.type character string, one of \code{c("none", "non-binding","binding")} or their
#'unique abbreviations. Specifies whether a futility boundary should not be used at all ("none"), or if it
#'is used, then whether the effect of the futility boundary should be included
#'in the efficacy power/type I error calculations ("binding"), or not ("non-binding").
#'@@param r_EN numeric vector; specifies the set of alternatives under which the design is optimized.
#' It is interpreted multiplicatively compared to the design alternative hypothesis.
#' Defaults to 1, implying minimization under the alternative hypothesis.
#'@@param r_EN.w numeric vector of length equal to that of \code{r_EN}. The weights of the
#'elements of \code{r_EN} in the optimality criterion. It will be normalized to a sum of 1.
#'Defaults to equal weights.
#'@@param spending (optional) numeric vector of stage-wise alpha-spending with positive values
#' summing to \code{sig.level} used as starting values for direct optimization. If used in conjunction with
#' \code{optim.method="none"} within \code{ocControl}, a non-optimized design using this alpha-spending
#' will be returned. The starting value will be selected automatically for the default value of NULL.
#'@@param control list of parameters controlling the estimation altorithm to replace the default
#'values returned by the \code{glsControl} function.
#'@@return an object of class \code{gsDesignOC}
#'@@author Aniko Szabo
#'@@references Szabo, A, Tarima, S (?) Operating-characteristic guided design of group-sequential trials.
#'@@keywords nonparametric design
#'@@examples
#'
#'gsDesignOC(2, rE.seq = c(1.5,1), rF.seq = c(-0.5,0),
#'           power.efficacy=0.8, power.futility=0.8, power=0.9,
#'           futility.type = "non-binding")
#'
#' ## non-constant target power
#'gsDesignOC(2, rE.seq = c(1,1), rF.seq = c(0,0),
#'           power.efficacy=0.5, power.futility=0.5, power=0.9,
#'           futility.type = "binding")

#' ## pre-specifiy alpha-spending
#'gsDesignOC(2, rE.seq = c(1.5,1), rF.seq = c(-0.5,0),
#'           power.efficacy=0.8, power.futility=0.8, power=0.9,
#'           futility.type = "non-binding", spending=c(0.005, 0.02),
#'           control = ocControl(optim.method="none"))
#'

#'@@name gsDesignOC

gsDesignOC <- function(n.stages, rE.seq, rF.seq=NULL, n.fix=1,
                       sig.level = 0.025, power=0.9,
                       power.efficacy=power, power.futility = power,
                       futility.type=c("none","non-binding","binding"),
                       r_EN = 1, r_EN.w = rep(1, length(r_EN)),
                       spending = NULL,
                       control=ocControl()){
  @< Check inputs @>

  # create skeleton gsDesignOC object
  res <- list(n.stages = n.stages, rE.seq=rE.seq, rF.seq = rF.seq,
              sig.level = sig.level, n.fix = n.fix,
              power=power, power.efficacy = power.efficacy, power.futility = power.futility,
              futility.type=futility.type, r_EN = r_EN, r_EN.w = r_EN.w)
  class(res) <- "gsDesignOC"

  if (n.stages == 1){
    alpha.seq <- sig.level
  } else if (control$optim.method == "none"){
    alpha.seq <- spending
  } else if (control$optim.method == "direct"){
  @< Find optimal design using direct search @>
  } else if (control$optim.method == "dynamic"){
  @< Find optimal design using recursive search @>
  } else {
    stop(sprintf("Optimization method %s is not available. Choose 'dynamic' or 'direct'.", control$optim.method))
  }

  res <- calc.bounds(x=res, alpha.seq)
  res$n <- n.fix * res$info /(qnorm(sig.level, lower.tail=FALSE) + qnorm(power))^2

  return(res)
}

@}


@D Check inputs @{
  if (length(rE.seq) == n.stages - 1){
    rE.seq <- c(rE.seq, 1)
  } else if (length(rE.seq) == n.stages){
    if (!isTRUE(all.equal(rE.seq[n.stages], 1)))
      stop("Invalid input: When 'rE.seq' has length `n.stages`, its last element should equal 1.")
  } else {
    stop("Invalid input: `rE.seq` should have length equal to `n.stages` or `n.stages-1`.")
  }
  if (any(diff(rE.seq) > 0))
    stop("Invalid input: `rE.seq` should form a non-increasing sequence")

  if (length(rF.seq) == n.stages - 1){
    rF.seq <- c(rF.seq, 0)
  } else if (length(rF.seq) == n.stages){
    if (!isTRUE(all.equal(rF.seq[n.stages], 0)))
      stop("Invalid input: When `rF.seq` has length `n.stages`, its last element should equal 0.")
  } else if (!is.null(rF.seq)) {
    stop("Invalid input: When not NULL, `rF.seq` should have length equal to `n.stages` or `n.stages-1`.")
  }
  if (any(diff(rF.seq) < 0))
    stop("Invalid input: `rF.seq` should form an non-decreasing sequence")

  if (length(r_EN.w) != length(r_EN))
      stop("Invalid input: `r_EN.w` should have the same length as `r_EN`")
  if (any(r_EN.w < 0))
      stop("Invalid input: `r_EN.w` should have only positive values")

  if (power >= 1 || power <= sig.level)
    stop("Invalid input: `power` should be between `sig.level` and 1 (exclusive).")

  if (length(power.efficacy) == 1) power.efficacy <- rep(power.efficacy, n.stages - 1)
  if (length(power.efficacy) == n.stages - 1){
    # add overall power at the end of interim power values
    power.efficacy <- c(power.efficacy, power)
  }
  if (length(power.efficacy) == n.stages){
     if (!isTRUE(all.equal(power.efficacy[n.stages], power)))
      stop("Invalid input: When 'power.futility' has length `n.stages`, its last element should be equal to `power`.")} else {
    stop("Invalid input: `power.efficacy` should have length 1, `n.stages-1`, or `n.stages`")
  }
  if (any(diff(c(sig.level, power.efficacy)) < 0))
    stop("Invalid input: `power.efficacy` should form a non-decreasing sequence with values  between `sig.level` and `power`.")

  if (length(power.futility) == 1) power.futility <- rep(power.futility, n.stages - 1)
  if (length(power.futility) == n.stages - 1){
    # add 1-alpha at the end of interim power values
    power.futility <- c(power.futility, 1-sig.level)
  }
  if (length(power.futility) == n.stages){
     if (!isTRUE(all.equal(power.futility[n.stages], 1-sig.level)))
      stop("Invalid input: When 'power.futility' has length `n.stages`, its last element should be equal to `1-sig.level`.")} else {
    stop("Invalid input: `power.futility` should have length 1, `n.stages-1`, or `n.stages`")
  }
  if (any(diff(c(0, power.futility)) < 0))
    stop("Invalid input: `power.futility` should form a non-decreasing sequence with values  between 0 and 1-`sig.level`.")

  futility.type <- match.arg(futility.type)
  if (futility.type != "none" && is.null(rF.seq))
    stop("Invalid input: `rF.seq` should be specified if a futility bound is requested")
  if (futility.type == "none" && !is.null(rF.seq))
    message("Note: `rF.seq` is ignored, as no futility bound is requested.")


  controlvals <- ocControl()
  if (!missing(control))
    controlvals[names(control)] <- control

  if (!is.null(spending)){
    if (any(spending <= 0))
      stop("Invalid input: `spending` should contain only positive values")
    if (!isTRUE(all.equal(sum(spending), sig.level)))
      stop("Invalid input: the values of `spending` should sum to `sig.level`")
  }
  if (control$optim.method == "none" & is.null(spending)){
    stop("Invalid input: with `optim.method` set to `none`, the alpha-spending sequence has to be specified in `spending`")
  }

@}

@O ../tests/testthat/test_main.R @{
test_that("Invalid rE-sequence inputs give an input error",{
  expect_error(gsDesignOC(n.stages = 2, rE.seq=c(3, 2, 1)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=c(3, 2)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=c(0.5, 1)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(2, 0.5)), "^Invalid input:")
})
test_that("Invalid rF-sequence inputs give an input error",{
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, rF.seq=c(-2, -1, 0)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, rF.seq=c(-2, -1)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, rF.seq=c(1, 0)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(3,2), rF.seq=c(-1, -2)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(3,2), rF.seq=NULL, futility.type="binding"), "^Invalid input:")
})
test_that("Integers and rounding in r-sequence do NOT give an input error",{
  res1 <- tryCatch(gsDesignOC(n.stages = 1, rE.seq = 1L), error=function(e)e)
  expect_true(!is(res1, "error") || !grepl(res1$message, "^Invalid input:"))

  res2 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = c(2L,1L)), error=function(e)e)
  expect_true(!is(res2, "error") || !grepl(res2$message, "^Invalid input:"))

  res3 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = c(2, 0.99999999999)), error=function(e)e)
  expect_true(!is(res3, "error") || !grepl(res3$message, "^Invalid input:"))

  res4 <- tryCatch(gsDesignOC(n.stages = 1, rE.seq = 1, rF.seq = 0L), error=function(e)e)
  expect_true(!is(res4, "error") || !grepl(res4$message, "^Invalid input:"))

  res5 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = 2, rF.seq = c(-1L,0L)), error=function(e)e)
  expect_true(!is(res5, "error") || !grepl(res5$message, "^Invalid input:"))

  res6 <- tryCatch(gsDesignOC(n.stages = 2, rE.seq = 2, rF.seq = c(-1, 0.000000001)), error=function(e)e)
  expect_true(!is(res6, "error") || !grepl(res6$message, "^Invalid input:"))
})
test_that("Invalid r_EN weights give an input error",{
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, r_EN=1, r_EN.w = c(1,2)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, r_EN=c(0, 1), r_EN.w = 1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 2, rE.seq=2, r_EN=c(0, 1), r_EN.w = c(-1, 1)), "^Invalid input:")
})
test_that("Invalid power inputs give an input error",{
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= 1.2), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, sig.level = 0.1, power= 0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= -0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= 0.8, power.efficacy = 0.9), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, power= 0.8, power.efficacy = -0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(2,1.5), power= 0.8, power.efficacy = c(0.8,0.5)), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, sig.level=0.1, power.futility = 0.95), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 1, rE.seq=2, sig.level=0.1, power.futility = -0.1), "^Invalid input:")
  expect_error(gsDesignOC(n.stages = 3, rE.seq=c(2,1.5), rF.seq=c(-1,-0.5), power= 0.8, power.futility = c(0.8,0.5)), "^Invalid input:")})

@}

\subsection{Print result}

@O ../R/gsDesignOC.R @{
#'@@export
#'@@rdname gsDesignOC

print.gsDesignOC <- function(x,...){
  cat("Group sequential design with ")
  cat(100 * x$power, "% power and", 100 * x$sig.level,
        "% Type I Error.\n\n")
  if (x$futility.type == "non-binding") {
    cat("Type I error computations assume trial continues if lower bound is crossed.\n")
    cat("Power computations assume trial stops if a bound is crossed.\n\n")
  }
  else if (x$futility.type =="binding")
    cat("Power and type I error computations assume trial stops if a bound is crossed.\n\n")

  # modified from print.gsDesign
  if (x$n.fix != 1) {
      ntxt <- "N "
      nval <- ceiling(x$n)
      nspace <- log10(x$n[x$n.stages])
      for (i in 1:nspace) {
          cat(" ")
      }
      cat("            ")
  }
  else {
      ntxt <- "Ratio*"
      nval <- round(x$n, 3)
      cat("        Sample\n")
      cat("         Size ")
  }
  if (x$futility.type != "none") {
    if (min(x$lower) < 0) {
      cat(" ")
    }
    cat(" ---Lower bounds--  ---Upper bounds--")
    y <- cbind(1:x$n.stages,
               nval,
               round(x$lower, 2),
               round(stats::pnorm(x$lower), 4),
               round(x$upper, 2),
               round(stats::pnorm(-x$upper), 4))
    colnames(y) <- c("Stage", ntxt, "Z   ",
            "   Nominal p", "  Z   ", "  Nominal p")
    }
  else {
    y <- cbind(1:x$n.stages,
               nval,
               round(x$upper, 2),
               round(stats::pnorm(-x$upper), 4))
    colnames(y) <- c("Stage", ntxt, "Z  ",   "Nominal p")
  }
  rownames(y) <- rep(" ", x$n.stages)
  cat("\n")
  print(y)
  if (x$n.fix == 1) {
    cat("* Sample size ratio compared to fixed design with no interim\n")
  }

  ox <- oc(x)
  cat("\nThe average expected sample size is ", format(round(ox$ave.EN,3), nsmall=4), ":\n")
  nmat <- cbind(round(x$r_EN, 2),
                round(x$r_EN.w,3),
                round(ox$EN.vec,3))
  colnames(nmat) <- c("Alternative", "Weight", "E{N}")
  rownames(nmat) <- rep("", length(x$r_EN))
  print(nmat)
  cat("\n\n")

  cat("Optimized under the restrictions on the cumulative probabilities of stopping at or before each stage:\n")
  cat("Stopping for efficacy\n")
  umat <- cbind(1:x$n.stages,
                round(x$rE.seq, 2),
                round(x$power.efficacy, 4),
                round(ox$efficacy.cumcross, 4))
  colnames(umat) <- c("Stage", "Alternative","Target","Actual")
  rownames(umat) <- rep("", x$n.stages)
  print(umat)

  if (x$futility.type != "none"){
    cat("Stopping for futility\n")
    lmat <- cbind(1:x$n.stages,
                  round(x$rF.seq, 2),
                  round(x$power.futility, 4),
                  round(ox$futility.cumcross, 4))
    colnames(lmat) <- c("Stage", "Alternative","Target","Actual")
    rownames(lmat) <- rep("", x$n.stages)
    print(lmat)
  }
  cat("\n")

  invisible(x)
}
@}

\subsection{gsDesign-compliant output object}

We will try to make an output object of class \texttt{gsDesign}, so the functions in that package will work with the output.

@O ../R/gsDesignOC.R @{
#'Convert fit to a \code{gsDesign} object
#'
#'The \code{as.gsDesign} converts a \code{gsDesignOC} object to a form that works with the
#' \code{gsDesign} package
#'
#'@@export
#'@@param x an object of class \code{gsDesignOC}
#'@@return an object of class \code{gsDesign}
#'@@author Aniko Szabo
#'@@keywords nonparametric design
#'@@examples
#'
#' require(gsDesign)
#' g <- gsDesignOC(n.stages=2, rE.seq=c(2,1))
#' gg <- as.gsDesign(g)
#' gg

as.gsDesign <- function(x){

  timing <- x$info / max(x$info)
  @< Create custom spending functions @>

  gres <- gsDesign::gsDesign(
    k = x$n.stages,
    test.type = switch(x$futility.type, none = 1, binding = 3, "non-binding"=4),
    alpha = x$sig.level,
    beta = 1 - x$power,
    timing = timing,
    n.fix = x$n.fix,
    sfu = upper_sf,
    sfupar = upper_par,
    sfl = lower_sf,
    sflpar = lower_par
  )

  return(gres)
}
@| as.gsDesign @}

We define a spending function with values equal to the optimal $\alpha_k$ at the optimal time-points, and is linearly interpolated between them. For the lower boundary we need to compute the actual $\beta$-spending first.
The \texttt{sfLinear} spending function expects a vector of parameters in which the first half defines the time-points and the second half the corresponding values scaled as the proportion of total error to be spent.

@D Create custom spending functions @{
  upper_sf <- gsDesign::sfLinear
  upper_par <- c(timing, cumsum(x$spending)/x$sig.level)

  if (x$futility.type != "none"){
    probs <- gsDesign::gsProbability(k=x$n.stages, n.I=x$info, a=x$lower, b=x$upper, theta=c(0,1))
    beta.spend <- cumsum(probs$lower$prob[,2])
    lower_par <- c(timing,  beta.spend/tail(beta.spend,1))
    lower_sf <- gsDesign::sfLinear
  } else {
    lower_sf <- NULL
    lower_par <- NULL
  }
@}


\subsection{Direct optimization}

@D Find optimal design using direct search @{
  @< Define optimization function @>

  if (n.stages == 2){
    oo <- optimize(.cp, interval=c(-5,5))

    y.res <- oo$minimum
    alpha.seq <- sig.level * exp(c(y.res,0)) / sum(exp(c(y.res,0)))
  } else {

    if (!is.null(spending)){
      y.start <- log(head(spending,-1)/tail(spending,1))
    } else {
      y.start <- -log(seq(n.stages, 2, by=-1))
    }

    oo <- optim(y.start, fn=.cp, method="Nelder-Mead")

    y.res <- oo$par
    alpha.seq <- sig.level * exp(c(y.res,0)) / sum(exp(c(y.res,0)))
  }
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

    rescp <- calc.bounds(x=res, alpha.seq = alpha.seq)
    dp <- oc(rescp, r_EN=r_EN, r_EN.w=r_EN.w)
    Q <- dp$ave.EN
    Q
  }
@}


\subsection{Recursive optimization}
@D Find optimal design using recursive search @{
@< Define recursive optimization function @>
  alpha.seq <- .opt.spend(x=res)
@}

The recursive optimization also implements the construction process in the proof of Theorem~\ref{Th:exists}, but will select the optimal $\Delta\alpha_K$ when adding the last stage.

For the first $K-1$ stages we will select the boundary $u_1,\ldots,u_{K-1}$ and information times $\I_{n_1},\ldots,\I_{n_{K-1}}$ to achieve over these $K-1$ stages stage-specific stopping probabilities $\pi_E$ at $\theta_{A1},\ldots,\theta_{A,K-2}$, overall power $\pi_E$ at $\theta_{A,K-1}$, and overall type I error $\alpha_{K-1} = \alpha-\Delta\alpha_K$.

The \texttt{.opt.spend} function will return the optimal alpha-spending sequence.

@D Define recursive optimization function @{
  .opt.spend <- function(x){
    if (x$n.stages == 1){
      return(x$sig.level)  # we spend it all
    }

    @< Find optimal DeltaAlpha to spend @>
  }
@}


The search for $0 < \Delta\alpha_K \leq \alpha$ will be parametrized via $y = \log(\Delta\alpha_K / \alpha)$ with $y\leq 0$ to restrict it to the correct range, and better spread out small values.

@D Find optimal DeltaAlpha to spend @{
  #create skeleton for a design with one fewer stages
  xx <- list(n.stages = x$n.stages - 1,
             rE.seq = head(rE.seq, -1),
             rF.seq = head(rF.seq, -1),
             power = tail(x$power.efficacy, 1),
             power.efficacy = head(x$power.efficacy, -1),
             power.futility = head(x$power.futility, -1),
             futility.type=x$futility.type)
  class(xx) <- "gsDesignOC"
@< Define function to calculate EN given DeltaAlpha @>
  res.spend <- optimize(en, interval=c(-5, 0))
  opt.spending <- attr(res.spend$objective, "spending")
  return(opt.spending)
@}

@D Define function to calculate EN given DeltaAlpha @{
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
#'@@param r_EN numeric vector of parameter values for which expected sample size calculation is performed. Defaults to values used to construct design.
#'@@param r_EN.w numeric vector of positive weights for each value of \code{r_EN}. Defaults to equal weights. Defaults to values used to construct design.
#'@@return a list with the following elements:
#'\describe{
#'\item{ave.EN}{numeric; weighted average of expected sample sizes under the requested alternatives}
#'\item{EN.vec}{numeric vector; expected sample sizes under each of the requested alternatives}
#'\item{efficacy.cumcross}{numeric vector; cumulative probability of stopping for efficacy at or before stage k
#'  under rE.seq[k], including 1 at the end.}
#'\item{futility.cumcross}{numeric vector}{numeric vector; cumulative probability of stopping for
#'  futility at or before stage k under rF.seq[k], including 0 at the end.}
#'}
#'@@author Aniko Szabo
#'@@keywords design
#'@@importFrom stats optimize optim setNames
#'@@importFrom utils head tail
#'@@examples
#'g <- gsDesignOC(n.stages = 2, rE.seq = c(1.5, 1), rF.seq = -1, power.efficacy=0.8,
#'           power.futility=0.8, power=0.9,
#'           futility.type = "non-binding", r_EN=c(1.5, 1, 0))
#'oc(g)
#'
#'oc(g, r_EN = c(2,1), r_EN.w = c(0.5,0.5))
#'

oc <- function(x, r_EN=x$r_EN,  r_EN.w = x$r_EN.w){

  if (length(r_EN.w) != length(r_EN))
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
#'@@param optim.method character; defines the optimization method used: \code{"dynamic"} uses a recursive dynamic algorithm #' which selects the alpha-spending one stage at a time, \code{"direct"} uses the "Nelder-Mead" algorithm from
#' \code{\link{optim}} to simultaneously search over all possible alpha-spending sequences. If set to \code{"none"}, then
#' no optimization is performed, and the designe corresponding to the starting value of the alpha-spending sequence,
#' which can be specified in the \code{spending} option of the call, will be returned.
#'@@return a list with components for each of the possible arguments
#'@@author Aniko Szabo
#'@@keywords design
#'@@examples
#'
#'ocControl(optim.method = "direct")

ocControl <- function(optim.method = c("direct", "dynamic","none")){
  optim.method <- match.arg(optim.method)
  list(optim.method = optim.method)
}
@}

\section{Boundary achieving the desired operating characteristics}

Given the alpha-spending sequence and using the probability targets defined in Table~\ref{T:OCdef}, we can derive the information times $\I_k$, and boundaries $u_k$ and $l_k$ using the recursive construction process described in the proof of Theorem~\ref{Th:exists}.



@O ../R/gsDesignOC.R
@{
#'Calculate the information times and efficacy/futility boundary values given the alpha-spending sequence and desired operating characteristics


#'@@param x a list with desired operating characteristics. Should have \code{rE.seq}, \code{sig.level},
#' \code{power}, \code{power.efficacy}, and \code{futility.type} components; and \code{rF.seq} and \code{power.futility}
#' components if lower boundary is requested. Objects of class \code{gsDesignOC} (potentially incomplete)  have
#' these components.
#'@@param alpha.seq numeric vector of stage-specific alpha-spending; values should add up to \code{x$sig.level}.
#'@@return The value of \code{x} augmented with the following components:
#'\describe{
#'\item{info}{numeric vector; information times for the analyses}
#'\item{lower}{numeric vector; lower Z-score boundary for the futility decisions}
#'\item{upper}{numeric vector; upper Z-score boundary for the efficacy decisions}
#'\item{spending}{numeric vector}{the value of alpha-spending sequence \code{alpha.seq}}
#'}
#'@@keywords internal
#'@@importFrom stats uniroot qnorm
#'@@importFrom utils head tail

calc.bounds <- function(x, alpha.seq){

  alpha.cum <- cumsum(alpha.seq)
  ub <- rep(20, x$n.stages)
  lb <- rep(-20, x$n.stages)
  ivec <- rep(NA, x$n.stages)

@< Define exceedance probability functions  @>
@< Construct first stage @>

  if (x$n.stages > 1){
    for (k in 2:x$n.stages){
    @< Construct next stage @>
    }
  }

  x$upper <- ub
  x$lower <- lb
  x$info <- ivec
  x$spending <- alpha.seq
  return(x)
}

@|calc.bounds @}

@O ../tests/testthat/test_components.R @{
test_that("calc.bounds for one stage returns fixed-sample test",{
  x <- list(n.stages=1, rE.seq=0.5, sig.level=0.05, power=0.8, power.efficacy=0.8, futility.type="none")
  xx <- calc.bounds(x, alpha.seq=0.05)
  expect_equal(xx$lower, -20)
  expect_equal(xx$upper, qnorm(0.95))
  expect_equal(xx$info, 24.7302289280791)
  expect_equal(xx$spending, 0.05)
 })
@}

The first stage is based on the fixed sample-size formula for alternative $\theta_{A1}$, power $\pi_{E1}$, and significance level $\Delta\alpha_1$.
@D Construct first stage @{
  ivec[1] <- ztest.I(delta=x$rE.seq[1], sig.level=alpha.seq[1], power=x$power.efficacy[1])
  ub[1] <- qnorm(alpha.seq[1], lower.tail=FALSE)
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
We will need to solve equations of the form $e(u)=c$, so we will actually define a function for the difference of the LHS and RHS. For non-binding futility, the lower bound is not needed under the null and is needed under the alternative hypotheses, so there is an option to control it.

@d Define upper bound exceedance  @{
exc <- function(u, I, stage,  theta, target, use_lb){
  aa <- if (use_lb) {c(head(lb, stage-1), -20)} else {rep(-20, stage)}
  gg <- gsDesign::gsProbability(k=stage,
                                theta=theta,
                                n.I =c(head(ivec, stage-1),I),
                                a = aa,
                                b = c(head(ub, stage-1), u))
   sum(gg$upper$prob[1:stage]) - target
}
@}

To find $u(I)$ we need to solve $a(u|\I) = \sum_{i=1}^k \Delta\alpha_i=\tilde{\alpha}_k$. Since $a(z_\alpha\mid\I_N) \geq \alpha$ [NOTE: I AM NOT SURE ABOUT THAT] and $a(\infty | \I_N) = \tilde{\alpha}_{k-1}$, the solution for fixed $\I_N$ is between $z_\alpha$ and $20$ (which is essentially $\infty$).

@d Define function to find u(I) @{
uI <- function(I, stage, use_lb){
  res <- uniroot(exc, interval=c(-20, 20),
                 I = I, stage=stage, theta=0, target = alpha.cum[stage],
                 use_lb = use_lb,
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


To find $I_k$, we need to solve $b(\I) = \pi_{Ek}$, where we set $\pi_{EK} = \pi$ for the last stage.
@d Determine target power and alternatives @{
power.target <- x$power.efficacy[k]
.th <- x$rE.seq[k]

@}


Note that $b(I_{k-1}) \leq \pi_E$ and $b(\I_{fix}(\theta_{Ak}-\theta_0, \tilde{\alpha}_k, \pi_E)) \leq \pi_E$, we can start the search at the maximum of these two values.

@D Find I for next stage @{

exc_i <- function(ii)exc(uI(ii, k, use_lb = (x$futility.type=="binding")), ii, stage=k, theta=.th,
                                target = power.target, use_lb = TRUE)

minI <- max(ztest.I(delta = .th, power=power.target, sig.level=alpha.cum[k]),
            ivec[k-1])
curr.exc <- exc_i(minI)

if (curr.exc > 0){
  # already past target power
  ivec[k] <- minI + .Machine$double.eps^0.25
  ub[k] <- lb[k] <- uI(ivec[k], k, use_lb=(x$futility.type=="binding"))
} else {
  resI <- uniroot(exc_i, interval=c(minI, 2*minI), extendInt = "upX")
  ivec[k] <- resI$root
  ub[k] <- lb[k] <- uI(resI$root, k, use_lb=(x$futility.type=="binding"))
}
@}

\subsection{Add lower bound}

The lower bound will usually be computed only when the information and upper bound for a stage has been determined. We will allow for vectorization over \texttt{theta} and \texttt{target}.

@d Define lower bound exceedance  @{
exc_low <- function(l, stage, theta, target, I=ivec[stage]){
  gg <- gsDesign::gsProbability(k=stage,
                                theta=theta,
                                n.I =c(head(ivec, stage-1), I),
                                a = c(head(lb, stage-1), l),
                                b = ub[1:stage])
  #lower$prob has rows by stage and columns by theta
   colSums(gg$lower$prob[1:stage, ,drop=FALSE]) - target
}
@}


@d Define function to find l @{
lI <- function(stage, theta, I=ivec[stage]){
  res <- uniroot(exc_low, interval=c(-20, qnorm(x$sig.level, lower.tail=FALSE)),
                 stage=stage, theta=theta, target = x$power.futility[stage], I=I,
                 extendInt = "upX")
 res$root
}
@}

We need to find the lower bound for stage $k-1$, before computing the size of the next stage. If the boundary is not binding, we need to overwrite $l_k:= u_k$ that was indicating the ``last" stage so far to $l_k=-\infty$.
@D Find lower bound for previous stage @{
  if (x$futility.type != "none"){
    lb[k-1] <- lI(stage=k-1, theta=x$rF.seq[k-1])
    if (x$futility.type == "binding"){
    @< Check beta-spending and adjust previous stage if needed @>
    }
  } else {
    lb[k-1] <- -20
  }
@}

With a binding boundary it is possible that the addition of $l_{K-1}$ overspends the type II error under $\theta_{Ak}$, and the target power for the next stage is not achievable anymore. More generally, the type II error could be overspent for any $\theta_{Aj}, j\geq k$. In this case, the information for the $(K-1)$st stage needs to be increased. Due to the monotonicity of $\theta_{Aj}$, we only need to monitor type II error spending under $\theta_A$. Specifically, we will ensure that each stage $k$ has at least $\beta_{min}(K-k+1)$ type II error available for a small value of $\beta_{min}$. Since \texttt{uniroot} solve with a tolerance of \texttt{.Machine\$double.eps\^0.25}, we need to set $\beta_{min}$ above that.


@D Check beta-spending and adjust previous stage if needed @{
  beta.min <- .Machine$double.eps^0.25
  .beta.limit <- beta.min * (x$n.stages - k + 1)
  typeII.overspent <- exc_low(lb[k-1], stage=k-1, theta=1, target = 1 - x$power - .beta.limit)
  if (typeII.overspent > 0){
    exc_low_i <- function(ii)
      exc_low(lI(I=ii, theta=x$rF.seq[k-1], stage=k-1), ii, stage=k-1, theta=1,
                     target = 1 - x$power - .beta.limit)
    resI_low <- uniroot(exc_low_i, interval=c(ivec[k-1], 2*ivec[k-1]), extendInt = "downX")
    ivec[k-1] <- resI_low$root
    ub[k-1] <- uI(I = resI_low$root, stage = k-1, use_lb = TRUE)
    lb[k-1] <- lI(I = resI_low$root, theta = x$rF.seq[k-1], stage = k-1, use_lb = TRUE)

  }
@}



\section{Utility help-functions}



@o ../R/Utility.R
@{
#'Fixed sample-size information for one-sided test
#'
#'The \code{ztest.I} function finds the information required of a one-sided Z-test with given power
#'and significance level
#'#ztest.I(delta=1, sig.level=0.05, power=0.9)
#'
#'@@param delta numeric; targeted standardized effect size
#'@@param sig.level numeric; one-sided significance level
#'@@param power numeric; target power
#'@@return numeric value with the (fractional) information achieving the target power
#'@@author Aniko Szabo
#'@@keywords internal
#'@@importFrom stats qnorm
#'

ztest.I <- function(delta, sig.level, power){
  za <- qnorm(sig.level, lower.tail=FALSE)
  zb <- qnorm(power)
  I <- (za+zb)^2 /delta^2
  I
}
@| ztest.I @}

@O ../tests/testthat/test_components.R @{
test_that("ztest.I is consistent with asbio::power.z.test",{
  i1 <- ztest.I(delta = 0.5, sig.level = 0.05, power=0.8)
  expect_equal(i1, 24.7302289280791)
  })

test_that("No data is needed if sig.level = power",{
  i2 <- ztest.I(delta = 0.5, sig.level = 0.1, power=0.1)
  expect_equal(i2, 0)
  })
@}

@{
#'Fixed sample-size effect size for one-sided test
#'
#'The \code{ztest.delta} function finds the standardized effect size for a one-sided Z-test with given power
#'and significance level
#'#ztest.delta(n=10, sig.level=0.05, power=0.9)
#'
#'@@param n numeric; sample size
#'@@param sig.level numeric; one-sided significance level
#'@@param power numeric; target power
#'@@return numeric value with the standardized effect size achieving the target power
#'@@author Aniko Szabo
#'@@keywords internal
#'@@importFrom stats qnorm
#'

ztest.delta <- function(n, sig.level, power){
  za <- qnorm(sig.level, lower.tail=FALSE)
  zb <- qnorm(power)
  delta <- abs(za+zb) / sqrt(n)
  delta
}
@| ztest.delta @}

@O ../tests/testthat/test_components.R @{
@%test_that("ztest.delta is consistent with ztest.I",{
@%  i1 <- ztest.I(delta = 0.5, sig.level = 0.05, power=0.8)
@%  d1 <- gsDesignOC:::ztest.delta(n = i1, sig.level = 0.05, power=0.8)
@%  expect_equal(d1, 0.5)
@%  })

@%test_that("No data is needed if sig.level = power",{
@%  d2 <- gsDesignOC:::ztest.delta(n=1, sig.level = 0.1, power=0.1)
@%  expect_equal(d2, 0)
@%})
@}



\end{document}
