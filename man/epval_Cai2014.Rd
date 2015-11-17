\name{epval_Cai2014}
\alias{epval_Cai2014}
\title{
Empirical permutation-based p-value of the Test Proposed by Cai, Liu, and Xia (2014)
}
\description{
Calculates p-value of the test for testing equality of two-sample high-dimensional mean vectors proposed by Cai, Liu, and Xia (2014) based on permutation.
}
\usage{
epval_Cai2014(sam1, sam2, perm.iter = 1000)
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
See the details in \code{\link{apval_Cai2014}}.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pval}{
  the permutation-based p-value of the test proposed by Cai, Liu, and Xia (2014).
}
}
\references{
Cai TT, Liu W, and Xia Y (2014). "Two-sample test of high dimensional means under dependence." \emph{JRSSB}, 76(2), 349--372.
}
\note{
The permutation technique assumes that the distributions of the two groups of samples are the same under the null hypothesis. Also, this function does not transform the data with their precision matrix (see Cai, Liu, and Xia, 2014). To calculate the p-value of the test statisic with transformation, users can input transformed samples to \code{sam1} and \code{sam2}.
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
epval_Cai2014(sam1, sam2, perm.iter = 10)
}
