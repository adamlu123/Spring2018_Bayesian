library("MCMCpack")
library("BoomSpikeSlab")
?lm.spike
dinvgamma(10,0.0001,0.0001)
hist(rinvgamma(100, 1, rate = 1))
rinvgamma(10,0.001,0.001)
# Load data 
dt <- read.csv("http://www.ics.uci.edu/~mguindan/teaching/stats225/prob1.csv",sep="", header = T)
y <- dt2[,1]
dt2 <- matrix(unlist(dt), ncol = 11, byrow = FALSE)
X <- cbind(rep(1,50),dt2[,2:11])
fit <- lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, dt)
summary(fit)
n <- 11
# Set value for parameters
tau = rep(0.33, n)
c = rep(10, n)
R = diag(1,n)
v=0.01; lambda = 0.01

# Init for beta, sig2, gamma
beta.LS <- coef(fit)
sig2 <- var(y)
gamma.i <- rep(1,n)

## 
epoch <- 1000
freq_gamma <- rep(0,11)
beta_placeholder <- matrix(0, epoch, 11)

for(i in seq(epoch)){
  beta.i <- update_beta(sig2,gamma.i,   X=X,R=diag(1,n), beta.LS=beta.LS,tau=tau,c=c)
  
  beta_placeholder[i,] <- beta.i
  beta.i <- as.matrix(beta.i)
  
  sig2.i <- update_sig2(beta.i)
  gamma.i <- update_gamma(gamma.i,beta.i,method='McCullogh',  c,tau)
  freq_gamma <- freq_gamma + gamma.i
  }
print(freq_gamma)
plot(beta_placeholder[,1], type='l')
acf(beta_placeholder[,1])
est.beta <- colMeans(beta_placeholder)
print(est.beta)
# part 3
model <- lm.spike(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, dt, niter=epoch)
plot.ts(model$beta)
plot(model)
summary(model)
summary(fit)
model$beta
