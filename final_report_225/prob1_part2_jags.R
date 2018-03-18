
## part 2
p=11
data=list(y=y,X=X,b_tau=rep(0.33, 11),p=11, c=rep(10, 11))
# init=list(tau=1,taub=1,alpha=0,betaT=rep(0,p),ind=rep(0,p))

c = rep(10, 11)

model_string <- "
model{

for (i in 1:50) { for(j in 1:p) {z[i,j]<-X[i,j]*b[j] }

eta[i] <-sum(z[i,])
y[i]~dnorm( eta[i], y_tau )

}

# prior
for(i in 1:p){
gamma[i]~dbern(0.5)
sig2_b[i] <- pow( gamma[i]*c[i]*b_tau[i] + (1-gamma[i])*b_tau[i], 2 )
b[i]~dnorm(0, 1/sig2_b[i])

}
y_tau~dgamma(0.01, 0.01)
# mu - dnorm(0, 1/10)
# y_pred~dnorm(mu, 1/4)

}
" 
model.spec<-textConnection(model_string)

jags <- jags.model(model.spec,
                   data = data,
                   n.chains=1,
                   n.adapt=100)

output <- coda.samples(jags,
                           c('b','gamma'	),
                           n.iter=10000,
                           thin=10
)
print(summary(output))
plot(output)

