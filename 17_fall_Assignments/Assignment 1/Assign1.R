# Question 1 
# Black Scholes ----
CallOption <- function(Stock,tau, Strike, rate, sigma) {
    d1 <- (log(Stock/Strike) + (rate + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price <- Stock*pnorm(d1) - Strike * exp(-rate*tau)*pnorm(d2)
   return(price)
 
    
}

PutOption <- function(Stock, tau, Strike, rate, sigma){
    
    d1 <- (log(Stock/Strike) + (rate + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price <- Strike * exp(-rate*tau) * pnorm(-d2) - Stock*pnorm(-d1)
    print(price)
}

OptionPrice <- function(Stock,tau, Strike, rate, sigma ){
 
    
    d1 <- (log(Stock/Strike) + (rate + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    CallPrice <- Stock*pnorm(d1) - Strike * exp(-rate*tau)*pnorm(d2)
    PutPrice <- Strike * exp(-rate*tau) * pnorm(-d2) - Stock*pnorm(-d1)
    print("Put price is")
    print(PutPrice)
    print("Call price is")
    print(CallPrice)
    }    
    
#Call put parity ----
PutCallParity <- function(Stock,tau, Strike, rate, sigma) {
LHS <- CallOption(Stock,tau, Strike, rate, sigma ) - PutOption(Stock,tau, Strike, rate, sigma )
RHS <- Stock - Strike *exp(-rate*tau)
print(RHS)
print(LHS)
return(LHS-RHS)

}
PutCallParity(100,30/252,100,0.05,0.2)
# Option data ----
maturity1 <- getOptionChain("FB","2017-03-17")
maturity2 <- getOptionChain("FB","2017-04-21")
maturity3 <- getOptionChain("FB","2017-09-15")
maturity1 <- maturity1["calls"]
maturity2 <- maturity2["calls"]
maturity3 <- maturity3["calls"]
month1 <- data.frame(maturity1)
month2 <- data.frame(maturity2)
month3 <- data.frame(maturity3)
month1 <- month1[1:20,]
month2 <- month2[1:20,]
month3 <- month3 [1:20,]
avg1 <- (month1$calls.Bid + month1$calls.Ask )/2
avg2 <- (month2$calls.Bid + month2$calls.Ask) / 2
avg3 <- (month3$calls.Bid + month3$calls.Ask) / 2
Stock1 <- getQuote("FB")

# Bisection Method ----
BisectionMethod <-  function(S, tau, Strike, r, market){
    up <- 2
    down <- 0
    mid <- (up + down) / 2
    i<- 0
    tol <- CallOption(S, tau, Strike, r, mid) - market 
     while(abs(tol) > 1e-04 && i<10000){
        if(tol < 0){
            down <- mid
        }else{
            up <- mid
            }
        mid<- (up + down)/2
        tol <- CallOption(S, tau,Strike, r, mid) - market 
        i <- i + 1
    }
    return(mid)
}

vol1 <- matrix(nrow = 1, ncol = 20)
vol2 <- matrix(nrow = 1, ncol = 20)
vol3 <- matrix(nrow = 1, ncol = 20)

for(i in 1:20) {
    vol1[i] = BisectionMethod(Stock1$Last,26/360,month1$calls.Strike[i],0.04,t(avg1[i]))
   vol2[i] = BisectionMethod(Stock1$Last,58/360,month2$calls.Strike[i],0.04,t(avg2[i]))
    vol3[i] = BisectionMethod(Stock1$Last,203/360,month3$calls.Strike[i],0.04,t(avg3[i]))
}

t(volat)
# Plotting ----
plot( month1$calls.Strike, t(vol1) , type = 'l', col = 'green',xlab = 'Strikes' , ylab = 'Volatility')
lines(month2$calls.Strike , t(vol2),type = 'l' , col = 'purple')
lines(month3$calls.Strike, t(vol3) , type = 'l' , col = 'red')

# Secant Method ----
SecantMethod <- function(S,tau,K,r,market){
 x1 <- 0
 x2 <- CallOption(S,tau,k,0.04,1)
 i <- 0
 while( i < 100) {
     ans[i] = market - (CallOption(S,tau,K,r,x1) - Stock1$Last)*
         (x2-x1)/(CallOption(S,tau,K,r,x2) - CallOption(S,tau,K,r,x1))
     x1 = x2
     x2 = ans[i]
     i = i +1
     
 }
 return(x2)
    
    }
impvol1 <- matrix(nrow = 1, ncol = 20)
impvol2 <- matrix(nrow = 1, ncol = 20)
impvol3 <- matrix(nrow = 1, ncol = 20)
impvol1[1] = SecantMethod(Stock1$Last,30/360,month1$calls.Strike[1],0.04,avg1[1])
impvol1
for (c in 1:20) {
    impvol1[c] = SecantMethod(Stock1$Last,26/360,month1$calls.Strike[c],0.04,avg1[c])
    impvol2[c] = SecantMethod(Stock1$Last,58/360,month2$calls.Strike[c],0.04,avg2[c])
    impvol1[c] = SecantMethod(Stock1$Last,203/360,month3$calls.Strike[c],0.04,avg3[c])
    }

impvolat <- cbind(impvol1,impvol2,impvol3)
t(impvolat) / 1000
plot(t(impvolat) , type = 'l')


# Greeks ----
tau <- 30/252
r <- 0.05
sigma <- 0.2
h <- 0.0001
Delta <- function(S,tau,k,r,sigma){
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    return(pnorm(d1))
}
Vega <- function(S,tau, k,r,sigma){
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    vega <-  S * sqrt(tau) * (1/sqrt(2*pi)) * exp(-d1^2/2)
    return(vega)
}
Gamma <- function(S,tau,k,r,sigma) {
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    gamma <- exp(-d1^2/2) / (S * sigma * sqrt(2*pi*tau))
   return (gamma)
}
Delta(100,30/252,100,5/100,0.2)
Delta(S,tau,k,r,sigma)
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
Vega
t(gamma)
t(delta)
t(vega)
# Greeks test ----
deltatest <- matrix(nrow = 1, ncol = 20)
gammatest <- matrix(nrow = 1, ncol = 20)
vegatest <- matrix(nrow =1, ncol = 20)
for(i in 1:20) {
    vegatest[i] <-  Vega(Stock1$Last,tau,month1$calls.Strike[i],0.04,vol1[i])
}
for(i in 1:20) {
    deltatest[i] <-  Vega(Stock1$Last,tau,month1$calls.Strike[i],0.04,vol1[i])
}
for(i in 1:20) {
    gammatesttest[i] <-  Vega(Stock1$Last,tau,month1$calls.Strike[i],0.04,vol1[i])
}
