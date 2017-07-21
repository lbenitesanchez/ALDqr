\name{diag.qr}
\alias{diag.qr}
\title{Diagnostics for Quantile Regression Using Asymmetric Laplace Distribution}
\description{Return case-deletion estimating the parameters in a quantile regression}

\usage{
  diag.qr(y,x,tau,theta)
}

\arguments{
  \item{y}{vector of responses}
  \item{x}{the design matrix}
  \item{tau}{the quantile to be estimated, this is generally a number strictly between 0 and 1.}
  \item{theta}{parameter estimated}
}

\value{Hessian and gradient matrix. Also the generalized cook distance (GDi), approximation of the likelihood distance (QDi) }

\author{Luis Benites Sanchez \email{lbenitesanchez@gmail.com}, Christian E. Galarza \email{cgalarza88@gmail.com}, Victor Hugo Lachos \email{hlachos@ime.unicamp.br}}


\references{
  [1] Koenker, R. W. (2005). Quantile Regression, Cambridge U. Press.
  [2] Yu, K. & Moyeed, R. (2001). Bayesian quantile regression. Statistics & Probability Letters, 54 (4), 437 to 447.
  [3] Kotz, S., Kozubowski, T. & Podgorski, K. (2001). The laplace distribution and generalizations: A revisit with applications to communications, economics, engineering, and finance. Number 183. Birkhauser. 
}

\examples{\dontrun{

##############################################################
### Graphic of the generalized Cook distance for data(AIS) ###
##############################################################
#Dados 
data(ais, package="sn")
attach(ais)
sexInd <- (sex=="female") + 0
x      <- cbind(1,LBM,sexInd)
y      <- BMI


#Percentile
 perc         <- 0.5

res           <- EM.qr(y,x,perc)
diag          <- diag.qr(y,x,perc,res$theta)
HessianMatrix <- diag$MatrizQ
Gradiente     <- diag$mdelta
GDI           <- c()
for (i in 1:202) {
 GDI[i] <- t(Gradiente[,i])%*%solve(-HessianMatrix)%*%Gradiente[,i]
}

  
#Seccion de los graficos
 par(mfrow = c(1,1))
 plot(seq(1:202),GDI,xlab='Index',ylab=expression(paste(GD[i])),main='p=0.1')
 abline(h=2*(4+1)/202,lty=2)
 identify(GDI,n=1) 

 plot(seq(1:202),GDI,xlab='Index',ylab=expression(paste(GD[i])),main='p=0.5')
 abline(h=2*(4+1)/202,lty=2)
 identify(GDI,n=1) 

 plot(seq(1:202),GDI,xlab='Index',ylab=expression(paste(GD[i])),main='p=0.9')
 abline(h=2*(4+1)/202,lty=2)
 identify(GDI,n=4) 


#############################################################
### Graphic of the likelihood displacemente for data(AIS) ###
#############################################################
#Dados 
 data(ais, package="sn"); attach(ais); sexInd<-(sex=="female")+0; x=cbind(1,LBM,sexInd); y=BMI

#Percentile
 perc          <- 0.9
 n             <- nrow(x)

 res           <- EM.qr(y,x,perc)
 
 thetaest      <- res$theta
 sigmaest      <- thetaest[4]
 betaest       <- matrix(thetaest[1:3],3,1)

 taup2         <- (2/(perc*(1-perc)))
 thep          <- (1-2*peGraphic of the generalized Cook distance for data(AIS)rc)/(perc*(1-perc))

 diag          <- diag.qr(y,x,perc,thetaest)

 HessianMatrix <- diag$MatrizQ
 Gradiente     <- diag$mdelta

 sigma         <- sigmaest
 beta          <- betaest 

 muc           <- (y-x%*%beta) 
 delta2        <- (y-x%*%beta)^2/(taup2*sigma)
 gamma2        <- (2+thep^2/taup2)/sigma

 vchpN         <- besselK(sqrt(delta2*gamma2), 0.5-1)
                  /(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
 vchp1         <- besselK(sqrt(delta2*gamma2), 0.5+1)
                  /(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))
 
 Q             <- -0.5*n*log(sigmaest)-0.5*(sigmaest*taup2)^{-1}*
                  (sum(vchpN*muc^2 - 2*muc*thep + vchp1*(thep^2+2*taup2)))  
 ########################################################
 theta_i       <- thetaest%*%matrix(1,1,202) +(-solve(HessianMatrix))%*%Gradiente
 sigmaest      <- theta_i[4,]
 betaest       <- theta_i[1:3,]
 sigma         <- sigmaest
 beta          <- betaest
 muc           <- (y-x%*%beta) 
 
 delta2        <- (y-x%*%beta)^2/(taup2*sigma)
 gamma2        <- (2+thep^2/taup2)/sigma
 
 vchpN         <- besselK(sqrt(delta2*gamma2), 0.5-1)
                  /(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
 vchp1         <- besselK(sqrt(delta2*gamma2), 0.5+1)
                  /(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))

 Q1 <- c()
 for (i in 1:202)
 {
   Q1[i] <- -0.5*n*log(sigmaest[i])-sum(vchpN[,i]*muc[,i]^2 - 2*muc[,i]*thep
     + vchp1[,i]*(thep^2+2*taup2))/(2*(sigmaest[i]*taup2))
 }

 ######################################################## 
 QDi <- 2*(-Q+Q1)

 #Depois de escolger perc guardamos os valores de  QDi
 QDi0.1 <- QDi
 QDi0.5 <- QDi
 QDi0.9 <- QDi

#Seccion de los graficos
 par(mfrow = c(1,1))
 plot(seq(1:202),QDi0.1,xlab='Index',ylab=expression(paste(QD[i])),main='p=0.1')
 abline(h=mean(QDi0.1)+3.5*sd(QDi0.1),lty=2)
 identify(QDi0.1,n=3) 

 plot(seq(1:202),QDi0.5,xlab='Index',ylab=expression(paste(QD[i])),main='p=0.5')
 abline(h=mean(QDi0.5)+3.5*sd(QDi0.5),lty=2)
 identify(QDi0.5,n=3) 

 plot(seq(1:202),QDi0.9,xlab='Index',ylab=expression(paste(QD[i])),main='p=0.9')
 abline(h=mean(QDi0.9)+3.5*sd(QDi0.9),lty=2)
 identify(QDi0.9,n=4) 
} 
}

\keyword{diagnostic}

