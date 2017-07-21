################################################################################
### DIAGNOSTICO: Ponderacao de casos    
################################################################################


diag.qr       <- function(y,x,tau,theta)
{
  ## PC= Ponderacao de casos
  p          <- ncol(x)
  n          <- nrow(x)
  taup2      <- (2/(tau*(1-tau)))
  thep       <- (1-2*tau)/(tau*(1-tau))
  
  beta       <- theta[1:p]
  sigma      <- theta[p+1]
  mu         <- x%*%beta
  muc        <- y-mu
  
  B          <- thep/(taup2*sigma)
  E          <- (thep^2+2*taup2)/(2*taup2*sigma^2)
  
  
  delta2     <- (y-x%*%beta)^2/(taup2*sigma)
  gamma2     <- (2+thep^2/taup2)/sigma
  
  
  DerBB      <- matrix(0,p,p)
  DerSS      <- 0
  DerBS      <- matrix(0,1,p)
  
  MatrizQ    <- matrix(0,nrow=(p+1),ncol=(p+1))
  GradQ      <- matrix(0,nrow=(p+1),ncol=n)
  
  
  vchpN1     <- besselK(sqrt(delta2*gamma2), 0.5-1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
  vchpN2     <- besselK(sqrt(delta2*gamma2), 0.5-2)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-2)
  vchpP1     <- besselK(sqrt(delta2*gamma2), 0.5+1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))
  vchpP2     <- besselK(sqrt(delta2*gamma2), 0.5+2)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(2) 
  
  for(i in 1:n)
  {
   Ai        <- muc[i]/(taup2*sigma)
   Ci        <- -0.5*(3*taup2*sigma+2*muc[i]*thep)/(taup2*sigma^2)
   Di        <- muc[i]^2/(2*taup2*sigma^2)
    
   DerBB     <- DerBB+ (-vchpN1[i]/(taup2*sigma))*((x[i,])%*%t(x[i,]))
   DerBS     <- DerBS+(-(vchpN1[i]*muc[i]-thep)/(taup2*sigma^2))*(x[i,])
   DerSS     <- DerSS+(1.5/sigma^2-(vchpN1[i]*(muc[i])^2-2*thep*muc[i]+vchpP1[i]*(2*taup2+thep^2))/(taup2*sigma^3))
    
   GradQ[,i] <- as.matrix(c((vchpN1[i]*Ai-B)*(x[i,]),(Ci+vchpN1[i]*Di+vchpP1[i]*E)),p+1,1)
 }
  
  MatrizQ[1:p,1:p] <- DerBB
  MatrizQ[p+1,1:p] <- (DerBS)
  MatrizQ[1:p,p+1] <- t(DerBS)
  MatrizQ[p+1,p+1] <- DerSS
  
  MatrizQ          <- (MatrizQ+t(MatrizQ))/2
  
  #############################################################
  ### Graphic of the likelihood displacemente for data      ###
  #############################################################

  thetaest      <- theta
  sigmaest      <- thetaest[p+1]
  betaest       <- matrix(thetaest[1:p],p,1)
  
  taup2         <- (2/(tau*(1-tau)))
  thep          <- (1-2*tau)/(tau*(1-tau))
  
  HessianMatrix <- MatrizQ
  Gradiente     <- GradQ
  
  sigma         <- sigmaest
  beta          <- betaest 
  
  muc           <- (y-x%*%beta) 
  delta2        <- (y-x%*%beta)^2/(taup2*sigma)
  gamma2        <- (2+thep^2/taup2)/sigma
  
  vchpN         <- besselK(sqrt(delta2*gamma2), 0.5-1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
  vchp1         <- besselK(sqrt(delta2*gamma2), 0.5+1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))
  
  Q             <- -0.5*n*log(sigmaest)-0.5*(sigmaest*taup2)^{-1}*(sum(vchpN*muc^2 - 2*muc*thep + vchp1*(thep^2+2*taup2)))  
  ########################################################
  theta_i       <- thetaest%*%matrix(1,1,n) +(-solve(HessianMatrix))%*%Gradiente
  sigmaest      <- theta_i[p+1,]
  betaest       <- theta_i[1:p,]
  sigma         <- sigmaest
  beta          <- betaest
  muc           <- (y-x%*%beta) 
  
  delta2        <- (y-x%*%beta)^2/(taup2*sigma)
  gamma2        <- (2+thep^2/taup2)/sigma
  
  vchpN         <- besselK(sqrt(delta2*gamma2), 0.5-1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))^(-1)
  vchp1         <- besselK(sqrt(delta2*gamma2), 0.5+1)/(besselK(sqrt(delta2*gamma2), 0.5))*(sqrt(delta2/gamma2))
  
  Q1 <- c()
  for (i in 1:n){Q1[i] <- -0.5*n*log(sigmaest[i])-sum(vchpN[,i]*muc[,i]^2 - 2*muc[,i]*thep + vchp1[,i]*(thep^2+2*taup2))/(2*(sigmaest[i]*taup2))}
  
  ######################################################## 
  QDi <- 2*(-Q+Q1)
  #############################################################
  
  #############################################################
  ## Graphic of the generalized Cook distance for data      ###
  #############################################################
  
  HessianMatrix <- MatrizQ
  Gradiente     <- GradQ
  GDi           <- c()
  for (i in 1:n) {GDi[i] <- t(Gradiente[,i])%*%solve(-HessianMatrix)%*%Gradiente[,i]}
  
  
  obj.out          <- list(MatrizQ = MatrizQ, mdelta=GradQ,QDi=QDi,GDi=GDi)
  return(obj.out)    
}

#theta   <- EM.qr(y,x,tau)$theta
#diag_qr(y,x,tau,theta)
  

