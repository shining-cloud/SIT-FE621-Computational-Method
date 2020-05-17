PutCallParity <- function(Stock,tau, Strike, rate, sigma) {
    LHS <- CallOption(Stock,tau, Strike, rate, sigma ) - PutOption(Stock,tau, Strike, rate, sigma )
    RHS <- Stock - Strike *exp(-rate*tau)
    print(RHS)
    print(LHS)
    return(LHS-RHS)
    
}
volat
strikes <- cbind(month1$calls.Strike,month2$calls.Strike,month3$calls.Strike) 
t(volat)
plot(t(volat), type = 'l')
PutCallParity(100,30/252,100,0.05,0.2)
strikes
