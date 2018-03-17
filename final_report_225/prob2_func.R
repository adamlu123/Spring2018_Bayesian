update_beta <- function(sig2,gamma,   X=X,R=diag(1,n),beta.LS=beta.LS,tau=tau,c=c){
  diagnoal <- ( (gamma==rep(0,n)) +(gamma==rep(1,n))*c ) * tau
  inv_D <- diag(1/diagnoal)
  A = solve(sig2^(-1)*t(X)%*%X+ inv_D%*%diag(1,n)%*%inv_D )
  mu <- A%*%t(X)%*%X%*%beta.LS/sig2
  beta.i <- mvrnorm(n = 1, mu, A)
  return(beta.i)
}

# update_beta(sig2,gamma,   X=X,R=diag(1,n),beta.LS=beta.LS,tau=tau,c=c)
update_sig2 <- function(beta.i){
  SSE <- sum(resid(fit)^2)
  sig2.i <- rinvgamma(1,(50+0)/2, SSE/2)
  return(sig2.i)
}
# update_sig2(beta.i)

update_gamma <- function(gamma,beta.i,method='McCullogh',  c,tau){
  for(j in 1:11){
    if(j==1){
    temp <- gamma # temp is 01 value 
    }
    
    if(method=='McCullogh'){
    diagonal <- rep(0,11)
    diagonal[j] <- (c[j]^(-2)-1)*tau[j]^(-2)
    dif_Var <- diag(diagonal)
    b_a <- c[j]*exp(0.5*t(beta.i)%*%dif_Var%*%beta.i) # b_a = b/a
    pj <- 1/(1+b_a)
    temp[j] <- rbinom(1, 1, pj)
    }else{
          inv_Var <- diag(( (gamma==rep(0,n))+(gamma==rep(1,n)*c) *tau)^(-2))
          pj <- c[j]^(-1)*exp(0.5*(1-1/c[j]^2 * t(beta.i)%*% inv_Var %*%beta.i  ) )
          temp[j] <- rbinom(1, 1, pj)
      
    }
    
  }
  return(temp)
}
update_gamma(gamma.i,beta.i,method='McCullogh',  c,tau)
