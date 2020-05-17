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
    return(price)
}
CallOption(100,30/252,100,5/100,0.2)
PutOption(100,30/252,100,5/100,0.2)

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
    vol1[i] = BisectionMethod(Stock1$Last,30/360,month1$calls.Strike[i],0.04,t(avg1[i]))
    vol2[i] = BisectionMethod(Stock1$Last,30/360,month2$calls.Strike[i],0.04,t(avg2[i]))
    vol3[i] = BisectionMethod(Stock1$Last,30/360,month3$calls.Strike[i],0.04,t(avg3[i]))
}

striketable1 <- matrix(nrow = 21,ncol=1)
voltable1 <- matrix(nrow = 20,ncol=1)
for (i in 1:20) {
    table1[i] = month1$calls.Strike[i]
}
for (i in 1:20){
    voltable1[i] <- vol1[i]
    
}
striketable1
voltable1
table1<- cbind(striketable1,voltable1)
table1
