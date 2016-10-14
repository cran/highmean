\name{apval_Cai2014}
\alias{apval_Cai2014}
\title{
Asymptotics-Based p-value of the Test Proposed by Cai et al (2014)
}
\description{
Calculates p-value of the test for testing equality of two-sample high-dimensional mean vectors proposed by Cai et al (2014) based on the asymptotic distribution of the test statistic.
}
\usage{
apval_Cai2014(sam1, sam2, eq.cov = TRUE)
}
\arguments{
  \item{sam1}{
  an n1 by p matrix from sample population 1. Each row represents a \eqn{p}-dimensional sample.
}
  \item{sam2}{
  an n2 by p matrix from sample population 2. Each row represents a \eqn{p}-dimensional sample.
}
  \item{eq.cov}{
  a logical value. The default is \code{TRUE}, indicating that the two sample populations have same covariance; otherwise, the covariances are assumed to be different. 
}
}

\details{
Suppose that the two groups of \eqn{p}-dimensional independent and identically distributed samples \eqn{\{X_{1i}\}_{i=1}^{n_1}} and \eqn{\{X_{2j}\}_{j=1}^{n_2}} are observed; we consider high-dimensional data with \eqn{p \gg n := n_1 + n_2 - 2}. Assume that the covariances of the two sample populations are \eqn{\Sigma_1 = (\sigma_{1, ij})} and \eqn{\Sigma_2 = (\sigma_{2, ij})}. The primary object is to test \eqn{H_{0}: \mu_1 = \mu_2} versus \eqn{H_{A}: \mu_1 \neq \mu_2}. Let \eqn{\bar{X}_{k}} be the sample mean for group \eqn{k = 1, 2}. For a vector \eqn{v}, we denote \eqn{v^{(i)}} as its \eqn{i}th element.

Cai et al (2014) proposed the following test statistic:
\deqn{T_{CLX} = \max_{i = 1, \ldots, p} (\bar{X}_1^{(i)} - \bar{X}_2^{(i)})^2/(\sigma_{1,ii}/n_1 + \sigma_{2, ii}/n_2),}
This test statistic follows an extreme value distribution under the null hypothesis.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{cov.assumption}{
  the equality assumption on the covariances of the two sample populations; this was specified by the argument \code{eq.cov}.
}
  \item{method}{
  this output reminds users that the p-values are obtained using the asymptotic distributions of test statistics.
}
  \item{pval}{
  the p-value of the test proposed by Cai et al (2014).
}
}
\references{
Cai TT, Liu W, and Xia Y (2014). "Two-sample test of high dimensional means under dependence." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{76}(2), 349--372.
}
\note{
This function does not transform the data with their precision matrix (see Cai et al, 2014). To calculate the p-value of the test statisic with transformation, users can use transformed samples for \code{sam1} and \code{sam2}.
}
\seealso{
\code{\link{epval_Cai2014}}
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
apval_Cai2014(sam1, sam2)

# the two sample populations have different covariances
true.cov1 <- 0.2^(abs(outer(1:p, 1:p, "-")))
true.cov2 <- 0.6^(abs(outer(1:p, 1:p, "-")))
sam1 <- mvrnorm(n = n1, mu = mu1, Sigma = true.cov1)
sam2 <- mvrnorm(n = n2, mu = mu2, Sigma = true.cov2)
apval_Cai2014(sam1, sam2, eq.cov = FALSE)
}
