UpandInFormula <- function(S,k,tau,r,sigma,H,div){
    v = r-div-(sigma^2/2)
    dbs <- function(S,k) {
        v = r-div-(sigma^2/2)
        ans = (log(S/k) + v*tau )/ (sigma*sqrt(tau))
        return(ans)
    }
    part1 = (H/S)^((2*v)/sigma^2)*(PutOption((H^2/S),k,tau,r,sigma) -
                                       PutOption((H^2/S),H,tau,r,sigma)
                                   +(H-k)*exp(-r*tau)*pnorm(-dbs(H,S)))
    
    part2 = CallOption(S,H,tau,r,sigma) +
        (H-k)*exp(-r*tau)*pnorm(dbs(S,H))
    
    
    part2
    answer = part1 + part2
    answer
    return(answer)
}



EUpandInCall <- function(S0,k,tau,r,sigma,div,H,N){
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
        if(S[i] >= H){
            C[i] <- max(0 , (S[i] - k))
        } else C[i] = 0
    }
    S
    C
    for(i in N:1) {
        for(j in 1:i){
            
            
            C[j] = disc * (pu * C[j+1] + pd * C[j])
        }
    }
    print(S)
    return(C[1])
} 
EUpandInCall(10,10,0.3,0.01,0.2,0,11,5) 
UpandInFormula(10,10,0.3,0.01,0.2,11,0)

# In - Out Parity ----

N = 200
Parity = EUpandInCall(10,10,0.3,0.01,0.2,0,11,N) + EUpandOutCall(10,10,0.3,0.01,0.2,0,11,N) -
    CallOption(10,10,0.3,0.01,0.2)
abs(Parity)
