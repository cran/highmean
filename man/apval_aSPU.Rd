\name{apval_aSPU}
\alias{apval_aSPU}
\title{
Asymptotics-Based p-values of the SPU and aSPU Tests
}
\description{
Calculates p-values of the sum-of-powers (SPU) and adaptive SPU (aSPU) tests based on the asymptotic distributions of the test statistics (Xu et al 2016+).
}
\usage{
apval_aSPU(sam1, sam2, pow = c(1:6, Inf), eq.cov = TRUE, cov.est,
           cov1.est, cov2.est, bandwidth, bandwidth1, bandwidth2,
           cv.fold = 5, norm = "F")
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
  \item{eq.cov}{
  a logical value. The default is \code{TRUE}, indicating that the two sample populations have same covariance; otherwise, the covariances are assumed to be different. 
}
  \item{cov.est}{
  a consistent estimate of the common covariance matrix when \code{eq.cov} is \code{TRUE}. This can be obtained from various apporoaches (e.g., banding, tapering, and thresholding; see Pourahmadi 2013). If not specified, this function uses a banding approach proposed by Bickel and Levina (2008) to estimate the covariance matrix.
}
  \item{cov1.est}{
  a consistent estimate of the covariance matrix of sample population 1 when \code{eq.cov} is \code{FALSE}. It is similar with the argument \code{cov.est}.
}
  \item{cov2.est}{
  a consistent estimate of the covariance matrix of sample population 2 when \code{eq.cov} is \code{FALSE}. It is similar with the argument \code{cov.est}.
}
  \item{bandwidth}{
  a vector of nonnegative integers indicating the candidate bandwidths to be used in the banding approach (Bickel and Levina, 2008) for estimating the common covariance when \code{eq.cov} is \code{TRUE}. This argument is effective only if \code{cov.est} is not provided. The default is a vector containing 50 candidate bandwidths chosen from \{0, 1, 2, ..., p\}.
}
  \item{bandwidth1}{
  similar with the argument \code{bandwidth}; it is used to specify candidate bandwidths for estimating the covariance of sample population 1 when \code{eq.cov} is \code{FALSE}.
}
  \item{bandwidth2}{
  similar with the argument \code{bandwidth}; it is used to specify candidate bandwidths for estimating the covariance of sample population 2 when \code{eq.cov} is \code{FALSE}.
}
  \item{cv.fold}{
  an integer greater than or equal to 2 indicating the fold of cross-validation. The default is 5. See page 211 in Bickel and Levina (2008).
}
  \item{norm}{
  a character string indicating the type of matrix norm for the calculation of risk function in cross-validation. This argument will be passed to the \code{\link[base]{norm}} function. The default is the Frobenius norm (\code{"F"}).
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
where \eqn{P_{SPU(\gamma)}} is the p-value of SPU(\eqn{\gamma}) test, and \eqn{\Gamma} is a candidate set of \eqn{\gamma}'s. Note that \eqn{T_{aSPU}} is no longer a genuine p-value. The asymptotic properties of the SPU and aSPU tests are studied in Xu et al (2016+).
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pow}{
  the powers \eqn{\gamma} used for the SPU tests.
}
  \item{opt.bw}{
  the optimal bandwidth determined by the cross-validation when \code{eq.cov} was \code{TRUE} and \code{cov.est} was not specified.
}
  \item{opt.bw1}{
  the optimal bandwidth determined by the cross-validation when \code{eq.cov} was \code{FALSE} and \code{cov1.est} was not specified.
}
  \item{opt.bw2}{
  the optimal bandwidth determined by the cross-validation when \code{eq.cov} was \code{FALSE} and \code{cov2.est} was not specified.
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
  the equality assumption on the covariances of the two sample populations; this was specified by the argument \code{eq.cov}.
}
  \item{method}{
  this output reminds users that the p-values are obtained using the asymptotic distributions of test statistics.
}
  \item{pval}{
  the p-values of the SPU tests and the aSPU test.
}
}
\references{
Bickel PJ and Levina E (2008). "Regularized estimation of large covariance matrices." \emph{The Annals of Statistics}, 36(1), 199--227.

Pan W, Kim J, Zhang Y, Shen X, and Wei P (2014). "A powerful and adaptive association test for rare variants." \emph{Genetics}, 197(4), 1081--1095.

Pourahmadi M (2013). \emph{High-Dimensional Covariance Estimation}. John Wiley & Sons, Hoboken, NJ.

Xu G, Lin L, Wei P, and Pan W (2016+). "An adaptive two-sample test for high-dimensional means."
}
%\author{
%Lifeng Lin <linl@umn.edu>
%}
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
# use true covariance matrix
apval_aSPU(sam1, sam2, cov.est = true.cov)
# fix bandwidth as 10
apval_aSPU(sam1, sam2, bandwidth = 10)
# use the optimal bandwidth from a candidate set
#apval_aSPU(sam1, sam2, bandwidth = 0:20)

# the two sample populations have different covariances
#true.cov1 <- 0.2^(abs(outer(1:p, 1:p, "-")))
#true.cov2 <- 0.6^(abs(outer(1:p, 1:p, "-")))
#sam1 <- mvrnorm(n = n1, mu = mu1, Sigma = true.cov1)
#sam2 <- mvrnorm(n = n2, mu = mu2, Sigma = true.cov2)
#apval_aSPU(sam1, sam2, eq.cov = FALSE,
#	bandwidth1 = 10, bandwidth2 = 10)
}
