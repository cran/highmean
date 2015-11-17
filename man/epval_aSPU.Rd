\name{epval_aSPU}
\alias{epval_aSPU}
\title{
Empirical permutation-based p-values of the SPU and aSPU Tests
}
\description{
Calculates p-values of the SPU (sum of powered score) and aSPU (adaptive SPU) tests based on permutation (Pan et al 2014).
}
\usage{
epval_aSPU(sam1, sam2, pow = c(1:6, Inf), perm.iter = 1000)
}
\arguments{
  \item{sam1}{
  an n1 by p matrix from sample population 1. Each row represents a \eqn{p}-dimensional sample.
}
  \item{sam2}{
  an n2 by p matrix from sample population 2. Each row represents a \eqn{p}-dimensional sample.
}
  \item{pow}{
  a numeric vector indicating the candidate values for the power \eqn{\gamma} in SPU tests. It should contain \code{Inf} and both odd and even integers. The default is \code{c(1:6, Inf)}.
}
  \item{perm.iter}{
  a numeric integer indicating the number of permutation iterations. The default is 1,000.
}
}

\details{
See the details in \code{\link{apval_aSPU}}.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pow}{
  the powers \eqn{\gamma} used for the SPU tests.
}
  \item{pval}{
  the permutation-based p-values of the SPU tests and the aSPU test.
}
}
\references{
Bickel PJ and Levina E (2008). "Regularized estimation of large covariance matrices." \emph{The Annals of Statistics}, 36(1), 199--227.

Pan W, Kim J, Zhang Y, Shen X, and Wei P (2014). "A powerful and adaptive association test for rare variants." \emph{Genetics}, 197(4), 1081--1095.

Xu G, Lin L, Wei P, and Pan W (2015+). "An adaptive two-sample test for high-dimensional means."
}
%\author{
%Lifeng Lin <linl@umn.edu>
%}
\note{
This function assumes that the two groups of samples have a common covariance matrix.
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
epval_aSPU(sam1, sam2, perm.iter = 10)
}
