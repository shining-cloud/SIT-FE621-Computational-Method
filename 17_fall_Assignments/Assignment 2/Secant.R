########## Question 1##############
#Binomial Tree ----
Binomialtree <- function(isCall, isAmerican, S0, k, tau, r, sigma, N){
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
    
    
    cp = ifelse(isCall, 1, -1)
    S <- c()
    S[1] = S0 * exp(N*dxd)
    
    
    for(i in 2:(N+1) ){
        S[i] = S[i-1] * exp(dxu - dxd) 
        if(isAmerican) {
            S[i]  = S[i-1] * edxud
        }
    }
    
    C <- c()
    for(i in 1:(N+1)) {
        C[i] <- max(0 , cp*(S[i] - k))
    }
    
    for(i in N:1) {
        for(j in 1:i){
            if (isAmerican == F){
                C[j] = disc * (pu * C[j+1] + pd * C[j])
            } else  {
                C[j] = dpd * C[j]  + dpu * C[j+1]
                S[j] = S[j] / edxd
                C[j] = max(C[j] , cp*(S[j] - k))
            }
        }
    }
    
    return(C[1])
} 
Binomialtree(T,F, 100, 100, 1, 0.06, 0.2, 3)
PutOption(100,100,1,0.06,0.2)
CallOption(100,100,1,0.06,0.2)
# Option data 
# 1 month data ----
datacall1 <- data.frame()
dataput1 <- data.frame()
datacall1<- read.csv("1C.csv")
dataput1 <- read.csv("1P.csv")
kc1 <- datacall1$calls.Strike
kp1 <- dataput1$puts.Strike
bc1 <- datacall1$calls.Bid
bp1 <- dataput1$puts.Bid
ac1 <- datacall1$calls.Ask
ap1 <- dataput1$puts.Ask
ac1
bc1
mc1 <- c()
mp1 <- c()
for(i in 1:10){
    mc1[i] = abs((ac1[i] + bc1[i])/ 2)
    mp1[i] = abs((ap1[i] + bp1[i])/ 2)
}
# 2 month data ----
datacall2 <- data.frame()
dataput2 <- data.frame()
datacall2<- read.csv("2C.csv")
dataput2 <- read.csv("2P.csv")
kc2 <- datacall2$calls.Strike
kp2 <- dataput2$puts.Strike
bc2 <- datacall2$calls.Bid
bp2 <- dataput2$puts.Bid
ac2 <- datacall2$calls.Ask
ap2 <- dataput2$puts.Ask
mc2 <- c()
mp2 <- c()
for(i in 1:10){
    mc2[i] = abs((ac2[i] + bc2[i])/ 2)
    mp2[i] = abs((ap2[i] + bp2[i])/ 2)
} 

# 3 month data ----
datacall3 <- data.frame()
dataput3 <- data.frame()
datacall3<- read.csv("3C.csv")
dataput3 <- read.csv("3P.csv")
kc3 <- datacall3$calls.Strike
kp3 <- dataput3$puts.Strike
bc3 <- datacall3$calls.Bid
bp3 <- dataput3$puts.Bid
ac3 <- datacall3$calls.Ask
ap3 <- dataput3$puts.Ask
mc3 <- c()
mp3 <- c()
for(i in 1:10){
    mc3[i] = abs((ac3[i] + bc3[i])/ 2)
    mp3[i] = abs((ap3[i] + bp3[i])/ 2)
} 
mc3
kc3
kp3
#Binomialtree(T,F,100,100,1,0.06,0.2,3)

S = 49.5

# Volatility question ----
 
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

BS <- function(S,k,tau,r,market,sigma,isCall) {
    
    if(isCall) {
        ans = CallOption(S,k,tau,r,sigma) - market
    } 
    else {
        ans = PutOption(S,k,tau,r,sigma) - market
    }
    return(ans)
}
Secant <- function(S,k,tau,r,market,isCall){
    err = 0.0001
    u = 1
    d = 0.01
    for(i in 1:1000) {
        fu = BS(S,k,tau,r,market,u,isCall) 
        fd = BS(S,k,tau,r,market,d,isCall)
        ans = u - (fu) * ((u-d) / (fu - fd)) 
        if (abs(u-d) < err){
            return(ans)
        }
     
        if(is.na(ans)){
       return(0)
    }
     if (is.infinite(ans)){
    return(u)
    }
      if (u > 500){
          return(d)
      }
        d = u
        u = ans
    }
}

impcall1 <- c() 
impcall2 <- c() 
impcall3 <- c() 

impcall1[1] = SecantMethod(S,kc1[1],30/252,0.0075,mc1[1], T)
impcall1[1]
impcall1[1] = SecantMethod(S,kc1[1],30/252,0.0075,mc1[1], T)
impcall1[1]
# 1 month call implied vol ----
for(i in 1:10){
    impcall1[i] =   Secant(S,kc1[i],30/252,0.0075,mc1[i],T)
}
impcall1
S
kc1
mc1
# 2 month call implied vol----
kc2
mc2
for(i in 1:10){
    impcall2[i] =   Secant(S,kc2[i],60/252,0.0075,mc2[i],T)
}
impcall2
kc2
mc2
# 3 month call implied vol----
for(i in 1:10){
    impcall3[i] =   Secant(S,kc3[i],90/252,0.0075,mc3[i],T)
    if(impcall3[i] < 0 )
        impcall3[i] = 0
}

impcall3
impput1 <- c()
impput2 <- c()
impput3 <- c()
# 1 month put implied vol----
for(i in 1:10){
    impput1[i] =   Secant(S,kp1[i],30/252,0.0075,mp1[i],F)
    if(impput1[i] < 0 )
        impput1[i] = 0
}
impput1
# 2 month put implied vol ----
for(i in 1:10){
    impput2[i] =   Secant(S,kp2[i],60/252,0.0075,mp2[i],F)
    if(impput2[i] < 0 )
        impput2[i] = 0
}
impput2
# 3 month put implied vol ----
kp3
for(i in 1:10){
    impput3[i] =   Secant(S,kp3[i],90/252,0.0075,mp3[i],F)
    if(impput3[i] < 0 ) {
        impput3[i] = 0
    }
}
impput3
length(impcall1)
################ European Options ###########
treecall1 <- c()
treecall2 <- c()
treecall3 <- c()
BSC1 <- c()
BSC2 <- c()
BSC3 <- c()
for (i in 1:10) {
    treecall1[i] =  Binomialtree(T,F,S,kc1[i],30/252,0.0075,impcall1[i],200)
    BSC1[i] = CallOption(S,kc1[i],30/252,0.0075,impcall1[i])
}
treecall1
mc1
BSC1

# 2 month Call ----
for (i in 1:10) {
    treecall2[i] =  Binomialtree(T,F,S,kc2[i],60/252,0.0075,impcall2[i],200)
    BSC2[i] = CallOption(S,kc2[i],60/252,0.0075,impcall2[i])
}
treecall2
mc2
BSC2
# 3 month Call  ----
for (i in 1:10) {
    treecall3[i] =  Binomialtree(T,F,S,kc3[i],90/252,0.0075,impcall3[i],200)
    BSC3[i] = CallOption(S,kc3[i],90/252,0.0075,impcall3[i])
}
treecall3
mc3
BSC3

# 1 month put ----
treeput1 <- c()
treeput2 <- c()
treeput3 <- c()
BSP1 <- c()
BSP2 <- c()
BSP3 <- c()
for (i in 1:10) {
    treeput1[i] =  Binomialtree(F,F,S,kp1[i],30/252,0.0075,impput1[i],200)
    BSP1[i] = PutOption(S,kp1[i], 30/252,0.0075,impput1[i])
}
treeput1
BSP1
mp1
# 2 month put  ----
for (i in 1:10) {
    treeput2[i] =  Binomialtree(F,F,S,kp2[i],60/252,0.0075,impput2[i],200)
    BSP2[i] = PutOption(S,kp2[i], 60/252,0.0075,impput2[i])
}
treeput2
BSP2
mp2
# 3 month put ----
for (i in 1:10) {
    treeput3[i] =  Binomialtree(F,F,S,kp3[i],90/252,0.0075,impput3[i],200)
    BSP3[i] = PutOption(S,kp3[i], 90/252,0.0075,impput3[i])
}
treeput3
BSP3
mp3
############### American Options ##########
# 1 month call  ----
Atreecall1 <- c()
Atreecall2 <- c()
Atreecall3 <- c()
ABSC1 <- c()
ABSC2 <- c()
ABSC3 <- c()
for (i in 1:10) {
    Atreecall1[i] =  Binomialtree(T,T,S,kc1[i],30/252,0.0075,impcall1[i],200)
    ABSC1[i] = CallOption(S,kc1[i],30/252,0.0075,impcall1[i])
}
Atreecall1
ABSC1
mc1
# 2 month call  ----
for (i in 1:10) {
    Atreecall2[i] =  Binomialtree(T,T,S,kc2[i],60/252,0.0075,impcall2[i],200)
    ABSC2[i] = CallOption(S,kc2[i],60/252,0.0075,impcall2[i])
}
Atreecall2
ABSC2
mc2
# 3 month call  ----
for (i in 1:10) {
    Atreecall3[i] =  Binomialtree(T,T,S,kc3[i],90/252,0.0075,impcall3[i],200)
    ABSC3[i] = CallOption(S,kc3[i],30/252,0.0075,impcall3[i])
}
Atreecall3
ABSC3
mc3

# 1 month put  ----
Atreeput1 <- c()
Atreeput2 <- c()
Atreeput3 <- c()
ABSP1 <- c()
ABSP2 <- c()
ABSP3 <- c()
for (i in 1:10) {
    Atreeput1[i] =  Binomialtree(F,T,S,kp1[i],30/252,0.0075,impput1[i],200)
    ABSP1[i] = PutOption(S,kp1[i],30/252,0.0075,impput1[i])
}
Atreeput1
ABSP1
mp1
# 2 month put  ----
for (i in 1:10) {
    Atreeput2[i] =  Binomialtree(F,T,S,kp2[i],60/252,0.0075,impput2[i],200)
    ABSP2[i] = PutOption(S,kp2[i],60/252,0.0075,impput2[i])
}
Atreeput2
ABSP2
mp2
# 3 month put ----
for (i in 1:10) {
    Atreeput3[i] =  Binomialtree(F,T,S,kp3[i],90/252,0.0075,impput3[i],200)
    ABSP3[i] = PutOption(S,kp3[i],90/252,0.0075,impput3[i])
}
Atreeput3
ABSP3
mp3

# Part d Absolute error ----
AbsoluteError <- function(S,k,tau,r= 0.0075,vol) {
    N <- c(10,20,30,40,50,100,150,200,250,300,350,400)
    err <- matrix(0,nrow = length(N), ncol = length(impput1))
    length(N)
    
    for(i in 1:length(N)){
        for(j in 1:length(impput1)){
            err[i,j] <- abs(Binomialtree(F,F,S,k[j],tau,r,vol[j],N[i]) - 
                                PutOption(S,k[j],tau,r,vol[j]));
            
        }
    }
    return(err)
}
AbsoluteError(S,kp1,tau,r=0.0075,impput1)
AbsoluteError(S,kp2,tau,r=0.0075,impput2)
AbsoluteError(S,kp3,tau,r=0.0075,impput3)
### As the number of time steps increases, the error decreases
###############Question 2 ##################
#Trinomial Tree---- ----
TrinomialTree <- function(isCall, isAmerican,S0, K, tau, r, sigma, N, div, dx) {
    dt = tau/N 
    edx <- exp(dx)
    nu = r - div - 0.5 * sigma^2
    pu = 0.5 * ( (sigma^2*dt + nu^2 *dt^2)/dx^2 + nu*dt/dx )
    pm = 1.0 -   (sigma^2*dt + nu^2 *dt^2)/dx^2  
    pd = 0.5 * ( (sigma^2*dt + nu^2 *dt^2)/dx^2 - nu*dt/dx )
    disc = exp(-r*dt)
    cp = ifelse(isCall, 1, -1)
    
    S <- c()
    S[1] <- S0 * exp(-N*dx)
    
    for(i in 2: (2*N + 1)) {
        S[i] = S[i-1] * edx
    }
    
    S <- rev(S)
    
    C <- matrix(0, nrow = (2*N+1), ncol = (N+1))
    
    for(i in 1:(2*N+1)) {
        C[i,(N+1)]  = ( max(0 , cp*(S[i] - K)) )
        
        
    }
    
    for(j in (N):1){
        for(i in ((N+1)-j+1):((N+1)+j-1)) {
            C[i,j] = disc*(pu*C[i-1,j+1] + pm*C[i,j+1] + pd*C[i+1,j+1])
            if (isAmerican) {
                C[i,j] = ( max(C[i,j] , cp*(S[i] - K)) )
            }
        }
    }
    C
    price <- C[N+1,1]
    
    return(price)
}

TrinomialTree(T,F,100,100,1,0.06,0.2,3,0.03,0.2)

S0 = k =  100 
tau = 1
sigma = 0.25
r = 0.06
N = 200
div = 0.03
dx = 0.2
EUCall =  TrinomialTree(T,F, S0,k,tau,r,sigma,200,div,dx)
EUCall
EUPut = TrinomialTree(F,F, S0,k,tau,r,sigma,200,div,dx)
EUPut
ACall = TrinomialTree(T,T, S0,k,tau,r,sigma,200,div,dx)
ACall
APut = TrinomialTree(F,T, S0,k,tau,r,sigma,200,div,dx)
APut


# Part b ----
# dx   ----
dx1c <- c()
dx1p <- c()
dx2c <- c()
dx2p <- c()
dx3c <- c()
dx3p <- c()
for (i in 1:10) {
    dx1c[i] = impcall1[i] * sqrt(3*(30/252)/N)
    dx1p[i] = impput1[i] * sqrt(3*(30/252)/N)
    dx3c[i] = impcall2[i] * sqrt(3*(90/252)/N)
    dx3p[i] = impput2[i] * sqrt(3*(90/252)/N)
    dx2c[i] = impcall2[i] * sqrt(3*(60/252)/N)
    dx2p[i] = impput2[i] * sqrt(3*(60/252)/N)
}



# 1 month trinomial tree European ----
EUtricall1 <- c()
BSEUcall1 <- c()
for(i in 1:10) {
    EUtricall1[i] = TrinomialTree(T,F,100,100,30/252,0.06,impput1[i],200,0.03,dx1p[i])
    BSEUcall1[i] = CallOption(100,30/252,100,0.06,impcall1[i])
}
EUtricall1
BSEUcall1
EUtriput1 <- c()
BSEUput1 <- c()
for(i in 1:10) {
    EUtriput1[i] = TrinomialTree(F,F,100,100,30/252,0.06,impput1[i],200,0.03,dx1p[i])
    BSEUput1[i] = PutOption(100,30/252,100,0.06,impput1[i])
}
EUtriput1
BSEUput1

# 2 month trinomial tree European ----
EUtricall2 <- c()
BSEUcall2 <- c()
for(i in 1:10) {
    EUtricall2[i] = TrinomialTree(T,F,100,100,60/252,0.06,impcall2[i],200,0.03,dx2c[i])
    BSEUcall2[i] = CallOption(100,60/252,100,0.06,impcall2[i])
}
EUtricall2
BSEUcall2
EUtriput2 <- c()
BSEUput2 <- c()
for(i in 1:10) {
    EUtriput2[i] = TrinomialTree(F,F,100,100,30/252,0.06,impput2[i],200,0.03,dx2p[i])
    BSEUput2[i] = PutOption(100,30/252,100,0.06,impput2[i])
}
EUtriput2
BSEUput2

# 3 month trinomial tree European ----
EUtricall3 <- c()
BSEUcall3 <- c()
for(i in 1:10) {
    EUtricall3[i] = TrinomialTree(T,F,100,100,90/252,0.06,impcall3[i],200,0.03,dx3c[i])
    BSEUcall3[i] = CallOption(100,90/252,100,0.06,impcall3[i])
}
EUtricall3
BSEUcall3
EUtriput3 <- c()
BSEUput3 <- c()
for(i in 1:10) {
    EUtriput3[i] = TrinomialTree(F,F,100,100,90/252,0.06,impput3[i],200,0.03,dx3p[i])
    BSEUput3[i] = PutOption(100,30/252,100,0.06,impput3[i])
}
impput3
dx3p
EUtriput3
BSEUput3

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
    C
    return(C)
} 

EUpandOutCall(100,100,1,0.06,0.2,0,95,3)
# part b----

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

EUpandOutFormula <- function(S,k,tau,r,sigma,H,div){
    v = r-div-(sigma^2/2)
    dbs <- function(S,k) {
        v = r-div-(sigma^2/2)
        ans = (log(S/k) + v*tau )/ (sigma*sqrt(tau))
        ans
        return(ans)
    }
    part1 = CallOption(S,tau,k,r,sigma) -  CallOption(S,tau,H,r,sigma) -
        (H-k)*exp(-r*tau)*pnorm(dbs(S,H))
    part1
    part2 = (H/S)^((2*v)/sigma^2)*(CallOption((H^2/S),tau,k,r,sigma) -
                                       CallOption((H^2/S),tau,H,r,sigma)
                                   -(H-k)*exp(-r*tau)*pnorm(dbs(H,S)))
    part2
    answer = part1 - part2
    answer
    return(answer)
}
EUpandOutFormula(100,100,1,0.06,0.2,95,0)





#Up and in formula ----
UpandInFormula <- function(S,k,tau,r,sigma,H,div){
    v = r-div-(sigma^2/2)
    dbs <- function(S,k) {
        v = r-div-(sigma^2/2)
        ans = (log(S/k) + v*tau )/ (sigma*sqrt(tau))
        ans
        return(ans)
    }
    part1 = (H/S)^((2*v)/sigma^2)*(PutOption((H^2/S),tau,k,r,sigma) -
                                       PutOption((H^2/S),tau,H,r,sigma)
                                   +(H-k)*exp(-r*tau)*pnorm(dbs(H,S)))
    
    part2 = CallOption(S,tau,H,r,sigma) +
        (H-k)*exp(-r*tau)*pnorm(dbs(S,H))
    
    
    part2
    answer = part1 - part2
    answer
    return(answer)
}
UpandInFormula(10,10,0.3,0.02,0.2,11,0)

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
            
            S[j] = disc * (pu*S[j+1] + pd * S[j])
            if(S[j] >= H){
                C[j] = disc * (pu * C[j+1] + pd * C[j])
                #C[j] <- max(0,(S[j] - k))
            } else
                C[j] = 0
        }
    }
    C
    return(C)
} 
EUpandInCall(100,100,1,0.06,0.2,0,95,3)
EUpandInCall(100,100,1,0.06,0.2,11,0,30)
EUpandInCall(10,10,0.3,0.01,0.2,0,11,30)

# d) American Up and In ----
AUpandInPut <- function(S0,k,tau,r,sigma,div,H,N){
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
        S[i] = S[i-1] * edxud 
        
    }
    S
    C <- c()
    for(i in 1:(N+1)) {
        if(S[i] > H){
            C[i] <- max(0 , (k - S[i]))
        } else C[i] = 0
    }
    for(i in N:1) {
        for(j in 1:i){
            
            S[j] = S[j] / exp(dxd)
            if(S[j] > H){
                C[j] = disc * (pu * C[j+1] + pd * C[j])
                C[j] = max(C[j], k - S[j])
            } else
                C[j] = 0
        }
    }
    C
    return(C)
}
AUpandInPut(100,100,1,0.06,0.2,0,95,3)


## up and in = put - up and out