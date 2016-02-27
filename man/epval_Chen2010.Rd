\name{epval_Chen2010}
\alias{epval_Chen2010}
\title{
Empirical Permutation- or Resampling-Based p-value of the Test Proposed by Chen and Qin (2010)
}
\description{
Calculates p-value of the test for testing equality of two-sample high-dimensional mean vectors proposed by Chen and Qin (2010) based on permutation or parametric bootstrap resampling.
}
\usage{
epval_Chen2010(sam1, sam2, eq.cov = TRUE, n.iter = 1000, cov1.est, cov2.est,
               bandwidth1, bandwidth2, cv.fold = 5, norm = "F", seeds)
}
\arguments{
  \item{sam1}{
  an n1 by p matrix from sample population 1. Each row represents a \eqn{p}-dimensional sample.
}
  \item{sam2}{
  an n2 by p matrix from sample population 2. Each row represents a \eqn{p}-dimensional sample.
}
  \item{eq.cov}{
  a logical value. The default is \code{TRUE}, indicating that the two sample populations have same covariance; otherwise, the covariances are assumed to be different. If \code{eq.cov} is \code{TRUE}, the permutation method is used to calculate p-values; otherwise, the parametric bootstrap resampling is used.
}
  \item{n.iter}{
  a numeric integer indicating the number of permutation/resampling iterations. The default is 1,000.
}
  \item{cov1.est}{
  This and the following arguments are only effective when \code{eq.cov = FALSE} and the parametric bootstrap resampling is used to calculate p-values. This argument specifies a consistent estimate of the covariance matrix of sample population 1 when \code{eq.cov} is \code{FALSE}. This can be obtained from various apporoaches (e.g., banding, tapering, and thresholding; see Pourahmadi 2013). If not specified, this function uses a banding approach proposed by Bickel and Levina (2008) to estimate the covariance matrix.
}
  \item{cov2.est}{
  a consistent estimate of the covariance matrix of sample population 2 when \code{eq.cov} is \code{FALSE}. It is similar with the argument \code{cov1.est}.
}
  \item{bandwidth1}{
  a vector of nonnegative integers indicating the candidate bandwidths to be used in the banding approach (Bickel and Levina, 2008) for estimating the covariance of sample population 1 when \code{eq.cov} is \code{FALSE}. This argument is effective when \code{cov1.est} is not provided. The default is a vector containing 50 candidate bandwidths chosen from \{0, 1, 2, ..., p\}.
}
  \item{bandwidth2}{
  similar with the argument \code{bandwidth1}; it is used to specify candidate bandwidths for estimating the covariance of sample population 2 when \code{eq.cov} is \code{FALSE}.
}
  \item{cv.fold}{
  an integer greater than or equal to 2 indicating the fold of cross-validation. The default is 5. See page 211 in Bickel and Levina (2008).
}
  \item{norm}{
  a character string indicating the type of matrix norm for the calculation of risk function in cross-validation. This argument will be passed to the \code{\link[base]{norm}} function. The default is the Frobenius norm (\code{"F"}).
}
  \item{seeds}{
  a vector of seeds for each permutation or parametric bootstrap resampling iteration; this is optional.
}
}

\details{
See the details in \code{\link{apval_Chen2010}}.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{opt.bw1}{
  the optimal bandwidth determined by the cross-validation when \code{eq.cov} was \code{FALSE} and \code{cov1.est} was not specified.
}
  \item{opt.bw2}{
  the optimal bandwidth determined by the cross-validation when \code{eq.cov} was \code{FALSE} and \code{cov2.est} was not specified.
}
  \item{cov.assumption}{
  the equality assumption on the covariances of the two sample populations; this was specified by the argument \code{eq.cov}.
}
  \item{method}{
  this output reminds users that the p-values are obtained using permutation or parametric bootstrap resampling.
}
  \item{pval}{
  the p-value of the test proposed by Chen and Qin (2010).
}
}
\references{
Bickel PJ and Levina E (2008). "Regularized estimation of large covariance matrices." \emph{The Annals of Statistics}, 36(1), 199--227.

Chen SX and Qin YL (2010). "A two-sample test for high-dimensional data with applications to gene-set testing." \emph{The Annals of Statistics}, 38(2), 808--835.

Pourahmadi M (2013). \emph{High-Dimensional Covariance Estimation}. John Wiley & Sons, Hoboken, NJ.
}
\examples{
#library(MASS)
#set.seed(1234)
#n1 <- n2 <- 50
#p <- 200
#mu1 <- rep(0, p)
#mu2 <- mu1
#mu2[1:10] <- 0.2
#true.cov <- 0.4^(abs(outer(1:p, 1:p, "-"))) # AR1 covariance
#sam1 <- mvrnorm(n = n1, mu = mu1, Sigma = true.cov)
#sam2 <- mvrnorm(n = n2, mu = mu2, Sigma = true.cov)
# increase n.iter to reduce Monte Carlo error.
#epval_Chen2010(sam1, sam2, n.iter = 10)

# the two sample populations have different covariances
#true.cov1 <- 0.2^(abs(outer(1:p, 1:p, "-")))
#true.cov2 <- 0.6^(abs(outer(1:p, 1:p, "-")))
#sam1 <- mvrnorm(n = n1, mu = mu1, Sigma = true.cov1)
#sam2 <- mvrnorm(n = n2, mu = mu2, Sigma = true.cov2)
# increase n.iter to reduce Monte Carlo error
#epval_Chen2010(sam1, sam2, eq.cov = FALSE, n.iter = 10,
#	bandwidth1 = 10, bandwidth2 = 10)
}