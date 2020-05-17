CallOption <- function(Stock,tau, Strike, rate, sigma) {
    d1 <- (log(Stock/Strike) + (rate + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price <- Stock*pnorm(d1) - Strike * exp(-rate*tau)*pnorm(d2)
    return(price)
    
    
}
tau <- 30/252
S <- 100
k <- 100
r <- 0.05
sigma <- 0.2
q = 0
h <- 0.0001

Delta <- function(S,tau,k,r,sigma){
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    return(exp(-q*tau)*pnorm(d1))
}
Vega <- function(S,tau, k,r,sigma){
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    vega <-  S * exp(-q*tau) *sqrt(tau) * (1/sqrt(2*pi)) * exp(-d1^2/2)
    return(vega)
}
Gamma <- function(S,tau,k,r,sigma) {
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    gamma <- (exp(-d1^2/2) * exp(-q*tau)) / (S * sigma * sqrt(2*pi*tau))
    return (gamma)
}
Delta(S,tau ,k,r,sigma)
Vega(S,tau,k,r,sigma)
Gamma(S,tau,k,r,sigma)

DeltaApprox <- function(S,tau,k,r,sigma) {
    Delta_approx <- (CallOption(S+h,tau,k,r,sigma) - CallOption(S,tau,k,r,sigma)) / h
    return(Delta_approx)
}
VegaApprox <- function(S,tau,k,r,sigma){
    Vega_approx<- (CallOption(S,tau,k,r,sigma+h) - CallOption(S,tau,k,r,sigma))/h
    return(Vega_approx)
}
GammaApprox <- function(S,tau,k,r,sigma){ 
    Gamma_approx <-(CallOption(S+2*h,tau,k,r,sigma)- 2*CallOption(S+h,tau,k,r,sigma) 
                    + CallOption(S,tau,k,r,sigma) )/h^2
    return(Gamma_approx)
}
VegaApprox(S,tau,k,r,sigma)
delta1 <- matrix(nrow = 1, ncol = 20)
delta2 <- matrix(nrow = 1, ncol = 20)
delta3<- matrix(nrow = 1, ncol = 20)
vega1 <- matrix(nrow = 1, ncol = 20)
vega2 <- matrix(nrow = 1, ncol = 20)
vega3 <- matrix(nrow = 1, ncol = 20)
gamma1 <-  matrix(nrow = 1, ncol = 20)
gamma2 <-  matrix(nrow = 1, ncol = 20)
gamma3 <-  matrix(nrow = 1, ncol = 20)
for(i in 1:20){
    delta1[i] <-  DeltaApprox(Stock1$Last,tau,month1$calls.Strike[i],0.04,vol1[i])
    delta2[i] <- DeltaApprox(Stock1$Last,tau,month2$calls.Strike[i],0.04,vol2[i])
    delta3[i] <- DeltaApprox(Stock1$Last,tau,month3$calls.Strike[i],0.04,vol3[i])
}
for(i in 1:20){
    vega1[i] <- VegaApprox(Stock1$Last,tau,month1$calls.Strike[i],0.04,vol1[i])
    vega2[i] <- VegaApprox(Stock1$Last,tau,month2$calls.Strike[i],0.04,vol2[i])
    vega3[i] <- VegaApprox(Stock1$Last,tau,month3$calls.Strike[i],0.04,vol3[i])
}
for (i in 1:20) {
    gamma1[i] <- GammaApprox(Stock1$Last,tau,month1$calls.Strike[i],0.04,vol1[i])
    gamma2[i] <- GammaApprox(Stock1$Last,tau,month2$calls.Strike[i],0.04,vol2[i])
    gamma3[i] <- GammaApprox(Stock1$Last,tau,month3$calls.Strike[i],0.04,vol3[i])
}
gamma <- cbind(gamma1,gamma2,gamma3)
delta <- cbind(delta1,delta2,delta3)
vega <- cbind(vega1,vega2,vega3)
t(gamma)
t(delta)
t(vega)
DeltaApprox(S,tau,k,r,sigma)
GammaApprox(S,tau,k,r,sigma)
VegaApprox(S,tau,k,r,sigma)
