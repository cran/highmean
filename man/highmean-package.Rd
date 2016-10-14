\name{highmean-package}
\alias{highmean-package}
\alias{highmean}
\docType{package}
\title{
Two-Sample Tests for High-Dimensional Mean Vectors
}
\description{
Provides various tests for comparing high-dimensional mean vectors in two sample populations.
}
\details{
Several two-sample tests for high-dimensional mean vectors have been proposed recently; see, e.g., Bai and Saranadasa (1996), Srivastava and Du (2008), Chen and Qin (2010), Cai et al (2014), and Chen et al (2014). However, these tests are powerful only against certain and limited alternative hypotheses. In practice, since the true alternative hypothesis is unknown, it is unclear how to choose one of these tests to yield high power. Accordingly, Pan et al (2014) and Xu et al (2016) proposed an adaptive test that may maintain high power across a wide range of situations; the asymptotic property of this test was also studied. This package provides functions to calculate p-values of the foregoing tests, using their asymptotic properties and the empirical (permutation or parametric bootstrap resampling) technique.
}
\author{
Lifeng Lin and Wei Pan

Maintainer: Lifeng Lin <linl@umn.edu>
}
\references{
Bai ZD and Saranadasa H (1996). "Effect of high dimension: by an example of a two sample problem." \emph{Statistica Sinica}, \bold{6}(2), 311--329.

Cai TT, Liu W, and Xia Y (2014). "Two-sample test of high dimensional means under dependence." \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, \bold{76}(2), 349--372.

Chen SX and Qin YL (2010). "A two-sample test for high-dimensional data with applications to gene-set testing." \emph{The Annals of Statistics}, \bold{38}(2), 808--835.

Chen SX, Li J, and Zhong PS (2014). "Two-Sample Tests for High Dimensional Means with Thresholding and Data Transformation." arXiv preprint arXiv:1410.2848.

Pan W, Kim J, Zhang Y, Shen X, and Wei P (2014). "A powerful and adaptive association test for rare variants." \emph{Genetics}, \bold{197}(4), 1081--1095.

Srivastava MS and Du M (2008). "A test for the mean vector with fewer observations than the dimension." \emph{Journal of Multivariate Analysis}, \bold{99}(3), 386--402.

Xu G, Lin L, Wei P, and Pan W (2016). "An adaptive two-sample test for high-dimensional means." \emph{Biometrika}, \bold{103}(3), 609--624.
}