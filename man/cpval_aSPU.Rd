\name{cpval_aSPU}
\alias{cpval_aSPU}
\title{
Permutation-And-Asymptotics-Based p-values of the SPU and aSPU Tests
}
\description{
Calculates p-values of the sum-of-powers (SPU) and adaptive SPU (aSPU) tests based on the combination of permutation method and asymptotic distributions of the test statistics (Xu et al, 2016).
}
\usage{
cpval_aSPU(sam1, sam2, pow = c(1:6, Inf), n.iter = 1000, seeds)
}
\arguments{
  \item{sam1}{
  an n1 by p matrix from sample population 1. Each row represents a \eqn{p}-dimensional sample.
}
  \item{sam2}{
  an n2 by p matrix from sample population 2. Each row represents a \eqn{p}-dimensional sample.
}
  \item{pow}{
  a numeric vector indicating the candidate powers \eqn{\gamma} in the SPU tests. It should contain \code{Inf} and both odd and even integers. The default is \code{c(1:6, Inf)}.
}
  \item{n.iter}{
  a numeric integer indicating the number of permutation iterations for calculating the means, variances, covariances of SPU test statistics' asymptotic distributions. The default is 1,000.
}
  \item{seeds}{
  a vector of seeds for each permutation iteration; this is optional.
}
}

\details{
Suppose that the two groups of \eqn{p}-dimensional independent and identically distributed samples \eqn{\{X_{1i}\}_{i=1}^{n_1}} and \eqn{\{X_{2j}\}_{j=1}^{n_2}} are observed; we consider high-dimensional data with \eqn{p \gg n := n_1 + n_2 - 2}. Assume that the covariances of the two sample populations are \eqn{\Sigma_1 = (\sigma_{1, ij})} and \eqn{\Sigma_2 = (\sigma_{2, ij})}. The primary object is to test \eqn{H_{0}: \mu_1 = \mu_2} versus \eqn{H_{A}: \mu_1 \neq \mu_2}. Let \eqn{\bar{X}_{k}} be the sample mean for group \eqn{k = 1, 2}. For a vector \eqn{v}, we denote \eqn{v^{(i)}} as its \eqn{i}th element.

For any \eqn{1 \le \gamma < \infty}, the sum-of-powers (SPU) test statistic is defined as:
\deqn{L(\gamma) = \sum_{i = 1}^{p} (\bar{X}_1^{(i)} - \bar{X}_2^{(i)})^\gamma.}
For \eqn{\gamma = \infty},
\deqn{L (\infty) = \max_{i = 1, \ldots, p} (\bar{X}_1^{(i)} - \bar{X}_2^{(i)})^2/(\sigma_{1,ii}/n_1 + \sigma_{2,ii}/n_2).}
The adaptive SPU (aSPU) test combines the SPU tests and improve the test power:
\deqn{T_{aSPU} = \min_{\gamma \in \Gamma} P_{SPU(\gamma)},}
where \eqn{P_{SPU(\gamma)}} is the p-value of SPU(\eqn{\gamma}) test, and \eqn{\Gamma} is a candidate set of \eqn{\gamma}'s. Note that \eqn{T_{aSPU}} is no longer a genuine p-value.

The asymptotic properties of the SPU and aSPU tests are studied in Xu et al (2016). When using the theoretical means, variances, and covarainces of \eqn{L (\gamma)} to calculate the p-values of SPU and aSPU tests (\eqn{1 \le \gamma < \infty}), the high-dimensional covariance matrix of the samples needs to be consistently estimated; such estimation is usually time-consuming.

Alternatively, assuming that the two sample groups have same covariance, the permutation method can be applied to efficiently estimate the means, variances, and covarainces of \eqn{L (\gamma)}'s asymptotic distributions, which then yield the p-values of SPU and aSPU tests based on the combination of permutation method and asymptotic distributions.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pow}{
  the powers \eqn{\gamma} used for the SPU tests.
}
  \item{spu.stat}{
  the observed SPU test statistics.
}
  \item{spu.e}{
  the asymptotic means of SPU test statistics with finite \eqn{\gamma} under the null hypothesis.
}
  \item{spu.var}{
  the asymptotic variances of SPU test statistics with finite \eqn{\gamma} under the null hypothesis.
}
  \item{spu.corr.odd}{
  the asymptotic correlations between SPU test statistics with odd \eqn{\gamma}.
}
  \item{spu.corr.even}{
  the asymptotic correlations between SPU test statistics with even \eqn{\gamma}.
}
  \item{cov.assumption}{
  the equality assumption on the covariances of the two sample populations; this reminders users that \code{cpval_aSPU()} assumes that the two sample groups have same covariance.
}
  \item{method}{
  this output reminds users that the p-values are obtained using the asymptotic distributions of test statistics.
}
  \item{pval}{
  the p-values of the SPU tests and the aSPU test.
}
}
\references{
Bickel PJ and Levina E (2008). "Regularized estimation of large covariance matrices." \emph{The Annals of Statistics}, \bold{36}(1), 199--227.

Pan W, Kim J, Zhang Y, Shen X, and Wei P (2014). "A powerful and adaptive association test for rare variants." \emph{Genetics}, \bold{197}(4), 1081--1095.

Pourahmadi M (2013). \emph{High-Dimensional Covariance Estimation}. John Wiley & Sons, Hoboken, NJ.

Xu G, Lin L, Wei P, and Pan W (2016). "An adaptive two-sample test for high-dimensional means." \emph{Biometrika}, \bold{103}(3), 609--624.
}
%\author{
%Lifeng Lin <linl@umn.edu>
%}
\note{
The permutation technique assumes that the distributions of the two sample populations are the same under the null hypothesis.
}
\seealso{
\code{\link{apval_aSPU}}, \code{\link{epval_aSPU}}
}
\examples{
library(MASS)
set.seed(1234)
n1 <- n2 <- 50
p <- 200
mu1 <- rep(0, p)
mu2 <- mu1
mu2[1:10] <- 0.2
true.cov <- 0.4^(abs(outer(1:p, 1:p, "-"))) # AR1 covariance
sam1 <- mvrnorm(n = n1, mu = mu1, Sigma = true.cov)
sam2 <- mvrnorm(n = n2, mu = mu2, Sigma = true.cov)
cpval_aSPU(sam1, sam2, n.iter = 100)
}
