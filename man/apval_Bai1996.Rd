\name{apval_Bai1996}
\alias{apval_Bai1996}
\title{
Asymptotic p-value of the Test Proposed by Bai and Saranadasa (1996)
}
\description{
Calculates p-value of the test for testing equality of two-sample high-dimensional mean vectors proposed by Bai and Saranadasa (1996) based on the asymptotic distribution of the test statistic.
}
\usage{
apval_Bai1996(sam1, sam2)
}
\arguments{
  \item{sam1}{
  an n1 by p matrix from sample population 1. Each row represents a \eqn{p}-dimensional sample.
}
  \item{sam2}{
  an n2 by p matrix from sample population 2. Each row represents a \eqn{p}-dimensional sample.
}
}

\details{
Suppose that the two groups of \eqn{p}-dimensional independent and identically distributed samples \eqn{\{X_{1i}\}_{i=1}^{n_1}} and \eqn{\{X_{2j}\}_{j=1}^{n_2}} are observed; we consider high-dimensional data with \eqn{p \gg n := n_1 + n_2 - 2}. Assume that the two groups share a common covariance matrix. The primary object is to test \eqn{H_{0}: \mu_1 = \mu_2} versus \eqn{H_{A}: \mu_1 \neq \mu_2}. Let \eqn{\bar{X}_{k}} be the sample mean for group \eqn{k = 1, 2}. Also, let \eqn{S = n^{-1} \sum_{k = 1}^{2} \sum_{i = 1}^{n_{k}} (X_{ki} - \bar{X}_k) (X_{ki} - \bar{X}_k)^T} be the pooled sample covariance matrix from the two groups.

Bai and Saranadasa (1996) proposed the following test statistic:
\deqn{T_{BS} = \frac{(n_1^{-1} + n_2^{-1})^{-1} (\bar{X}_1 - \bar{X}_2)^T (\bar{X}_1 - \bar{X}_2) - tr S}{\sqrt{2 n (n + 1) (n - 1)^{-1} (n + 2)^{-1} [tr S^2 - n^{-1} (tr S)^2]}},}
and its asymptotic distribution is normal under the null hypothesis.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pval}{
  the asymptotic p-value of the test proposed by Bai and Saranadasa (1996).
}
}
\references{
Bai ZD and Saranadasa H (1996). "Effect of high dimension: by an example of a two sample problem." \emph{Statistica Sinica}, 6(2), 311--329.
}
\note{
The asymptotic distribution of the test statistic was derived under normality assumption in Bai and Saranadasa (1996). Also, this function assumes that the two groups of samples have a common covariance matrix.
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
sam2 <- mvrnorm(n = n1, mu = mu2, Sigma = true.cov)
apval_Bai1996(sam1, sam2)
}
