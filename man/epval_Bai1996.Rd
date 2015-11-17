\name{epval_Bai1996}
\alias{epval_Bai1996}
\title{
Empirical permutation-based p-value of the Test Proposed by Bai and Saranadasa (1996)
}
\description{
Calculates p-value of the test for testing equality of two-sample high-dimensional mean vectors proposed by Bai and Saranadasa (1996) based on permutation.
}
\usage{
epval_Bai1996(sam1, sam2, perm.iter = 1000)
}
\arguments{
  \item{sam1}{
  an n1 by p matrix from sample population 1. Each row represents a \eqn{p}-dimensional sample.
}
  \item{sam2}{
  an n2 by p matrix from sample population 2. Each row represents a \eqn{p}-dimensional sample.
}
  \item{perm.iter}{
  a numeric integer indicating the number of permutation iterations. The default is 1,000.
}
}

\details{
See the details in \code{\link{apval_Bai1996}}.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pval}{
  the permutation-based p-value of the test proposed by Bai and Saranadasa (1996).
}
}
\references{
Bai ZD and Saranadasa H (1996). "Effect of high dimension: by an example of a two sample problem." \emph{Statistica Sinica}, 6(2), 311--329.
}
\note{
The permutation technique assumes that the distributions of the two groups of samples are the same under the null hypothesis.
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
# increase perm.iter to reduce Monte Carlo error.
epval_Bai1996(sam1, sam2, perm.iter = 10)
}