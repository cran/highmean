\name{apval_Chen2014}
\alias{apval_Chen2014}
\title{
Asymptotic p-value of the Test Proposed by Chen, Li, and Zhong (2014)
}
\description{
Calculates p-value of the test for testing equality of two-sample high-dimensional mean vectors proposed by Chen, Li, and Zhong (2014) based on the asymptotic distribution of the test statistic.
}
\usage{
apval_Chen2014(sam1, sam2)
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
Suppose that the two groups of \eqn{p}-dimensional independent and identically distributed samples \eqn{\{X_{1i}\}_{i=1}^{n_1}} and \eqn{\{X_{2j}\}_{j=1}^{n_2}} are observed; we consider high-dimensional data with \eqn{p \gg n := n_1 + n_2 - 2}. Assume that the two groups share a common covariance matrix. The primary object is to test \eqn{H_{0}: \mu_1 = \mu_2} versus \eqn{H_{A}: \mu_1 \neq \mu_2}. Let \eqn{\bar{X}_{k}} be the sample mean for group \eqn{k = 1, 2}. Also, let \eqn{S = n^{-1} \sum_{k = 1}^{2} \sum_{i = 1}^{n_{k}} (X_{ki} - \bar{X}_k) (X_{ki} - \bar{X}_k)^T} be the pooled sample covariance matrix from the two groups. Moreover, for a vector \eqn{v}, we denote \eqn{v^{(i)}} as its \eqn{i}th element.

Chen, Li, and Zhong (2014) proposed removing estimated zero components in the mean difference through thresholding; they considered
\deqn{T_{CLZ}(s) = \sum_{i = 1}^{p} \left\{ \frac{n_1 n_2 (\bar{X}_1^{(i)} - \bar{X}_2^{(i)})^2}{(n_1 + n_2) \sigma_{ii}} - 1 \right\} I \left\{ \frac{n_1 n_2 (\bar{X}_1^{(i)} - \bar{X}_2^{(i)})^2}{(n_1 + n_2) \sigma_{ii}} > \lambda_{p} (s) \right\},}
where the threshold level is \eqn{\lambda_p(s) := 2 s \log p} and \eqn{I(\cdot)} is the indicator function. Since an optimal choice of the threshold is unknown, they proposed trying all possible threshold values, then choosing the most significant one as their final test statistic:
\deqn{T_{CLZ} = \max_{s \in (0, 1 - \eta)} \{ T_{CLZ}(s) - \hat{\mu}_{T_{CLZ}(s), 0}\}/\hat{\sigma}_{T_{CLZ}(s), 0},}
where \eqn{\hat{\mu}_{T_{CLZ}(s), 0}} and \eqn{\hat{\sigma}_{T_{CLZ}(s), 0}} are estimates of the mean and standard deviation of \eqn{T_{CLZ}(s)} under the null hypothesis. They derived its asymptotic null distribution as an extreme value distribution.
}

\value{
A list including the following elements:
  \item{sam.info}{
  the basic information about the two groups of samples, including the samples sizes and dimension.
}
  \item{pval}{
  the asymptotic p-value of the test proposed by Chen, Li, and Zhong (2014).
}
}
\references{
Chen SX, Li J, and Zhong PS (2014). "Two-Sample Tests for High Dimensional Means with Thresholding and Data Transformation." arXiv preprint arXiv:1410.2848.
}
\note{
This function assumes that the two groups of samples have a common covariance matrix. Also, this function does not transform the data with their precision matrix (see Chen, Li, and Zhong, 2014). To calculate the p-value of the test statisic with transformation, users can use transformed samples for \code{sam1} and \code{sam2}.
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
apval_Chen2014(sam1, sam2)
}
