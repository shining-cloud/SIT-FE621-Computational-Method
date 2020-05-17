# Black Scholes ----
CallOption <- function(S,k, tau, r, sigma) {
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price =  S*pnorm(d1) - k * exp(-r*tau)*pnorm(d2)
    return(price)
    
    
}

PutOption <- function(S,k, tau, r, sigma){
    
    d1 <- (log(S/k) + (r + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price <- k * exp(-r*tau) * pnorm(-d2) - S*pnorm(-d1)
    return(price)
}

# ############Question 3 Exotic Options##########
#part a ----
EUpandOutCall <- function(S0,k,tau,r,sigma,div,H,N){
    dt <- tau / N
    nu <- r - (0.5 * sigma^2)
    dxu <- sqrt(sigma^2 * dt + (nu * dt)^2)
    dxd <-  -dxu
    pu <- 0.5 + 0.5 * ((nu * dt) / dxu)
    pd <- 1 - pu
    disc <- exp(-r*dt)
    
    dpu <-  disc * pu
    dpd <- disc * pd
    edxud <- exp( dxu - dxd)
    edxd <- exp(dxd) 
    
    
    S <- c()
    S[1] = S0 * exp(N*dxd)
    S
    
    for(i in 2:(N+1) ){
        S[i] = S[i-1] * exp(dxu - dxd) 
        
    }
    S
    C <- c()
    for(i in 1:(N+1)) {
        if(S[i] < H){
            C[i] <- max(0 , (S[i] - k))
        } else C[i] = 0
    }
    S
    C
    for(i in N:1) {
        for(j in 1:i){
            S[j] = disc * (pu*S[j+1] + pd * S[j])
            # S[j] = S[j] / exp(dxd)
            if(S[j] < H){
                C[j] = disc * (pu * C[j+1] + pd * C[j])
            } else
                C[j] = 0
        }
    }
  #  print(S)
    return(C[1])
} 

EUpandOutCall(10,10,0.3,0.01,0.2,0,11,32)
# part b----


EUpandOutFormula <- function(S,k,tau,r,sigma,H,div){
    v = r-div-(sigma^2/2)
    dbs <- function(S,k) {
        v = r-div-(sigma^2/2)
        ans = (log(S/k) + v*tau )/ (sigma*sqrt(tau))
        ans
        return(ans)
    }
    part1 = CallOption(S,k,tau,r,sigma) -  CallOption(S,H,tau,r,sigma) -
        (H-k)*exp(-r*tau)*pnorm(dbs(S,H))
    part1
    part2 = (H/S)^((2*v)/sigma^2)*(CallOption((H^2/S),k,tau,r,sigma) -
                                       CallOption((H^2/S),H,tau,r,sigma)
                                   -(H-k)*exp(-r*tau)*pnorm(dbs(H,S)))
    part2
    answer = part1 - part2
    answer
    return(answer)
}
EUpandOutFormula(10,10,0.3,0.01,0.2,11,0)
EUpandOutCall(10,10,0.3,0.01,0.2,0,11,32)
