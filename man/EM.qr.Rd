\name{EM.qr}
\alias{EM.qr}
\title{Quantile Regression Using Asymmetric Laplace Distribution}
\description{Return estimating the parameters in a quantile regression}

\usage{
EM.qr(y, x = NULL, tau = NULL, error = 0.000001 , iter = 2000, envelope=FALSE)
}

\arguments{
  \item{y}{vector of responses}
  \item{x}{the design matrix}
  \item{tau}{the quantile to be estimated, this is generally a number strictly between 0 and 1.}
  \item{error}{the covergence maximum error}
  \item{iter}{maximum iterations of the EM algorithm.}
  \item{envelope}{confidence envelopes for a curve based on bootstrap replicates }
}

\value{Estimated parameter for a quantile regression fit,standard error, log-likelihood.}

\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com}, Christian E. Galarza \email{cgalarza88@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br}}

\references{
[1] Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press.

[2] Yu, K. & Moyeed, R. (2001). Bayesian quantile regression. Statistics & Probability Letters, 54 (4), 437 to 447.

[3] Kotz, S., Kozubowski, T. & Podgorski, K. (2001). The laplace distribution and generalizations: A revisit with applications to communications, economics, engineering, and finance. Number 183. Birkhauser. 
}

\examples{
data(ais, package="sn")
attach(ais)
sexInd <- (sex=="female") + 0
x      <- cbind(1,LBM,sexInd)
y      <- BMI
tau    <- 0.5

## EM.qr
EM.qr(y,x,tau)


}

\keyword{regression}

