########## Question 1##############
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
Binomialtree(F,T, 100, 100, 1, 0.06, 0.2, 500)
CallOption(100,100,1,0.06,0.2)
PutOption(100,100,1,0.06,0.2)

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
    if(impcall2[i] < 0 )
    impcall2[i] = 0
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
################ European Options Binomial tree ########### 
treecall1 <- c()
treecall2 <- c()
treecall3 <- c()
BSC1 <- c()
BSC2 <- c()
BSC3 <- c()
# 1 month call ----
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
############### American Options Binomial tree ##########
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
# strike, vol, m ,bs, ebinomial, abinomail
#i,


# matrices for results Binomial results ----
call1 <- matrix(0,nrow = 10, ncol = 6)
call2 <- matrix(0,nrow = 10, ncol = 6)
call3 <- matrix(0,nrow = 10, ncol = 6)
put1 <- matrix(0,nrow = 10, ncol = 6)
put2 <- matrix(0,nrow = 10, ncol = 6)
put3 <- matrix(0,nrow = 10, ncol = 6)
for(i in 1:10) {
    call1[i,1] = kc1[i]
    call1[i,2] = impcall1[i]
    call1 [i,3] = mc1[i]
    call1[i,4] = BSC1[i]
    call1[i,5] = treecall1[i]
    call1[i,6] = Atreecall1[i]
    
    call2[i,1] = kc2[i]
    call2[i,2] = impcall2[i]
    call2 [i,3] = mc2[i]
    call2[i,4] = BSC2[i]
    call2[i,5] = treecall2[i]
    call2[i,6] = Atreecall2[i]
    
    call3[i,1] = kc3[i]
    call3[i,2] = impcall3[i]
    call3 [i,3] = mc3[i]
    call3[i,4] = BSC3[i]
    call3[i,5] = treecall3[i]
    call3[i,6] = Atreecall3[i]
    
    put1[i,1] = kp1[i]
    put1[i,2] = impput1[i]
    put1 [i,3] = mp1[i]
    put1[i,4] = BSP1[i]
    put1[i,5] = treeput1[i]
    put1[i,6] = Atreeput1[i]
    
    put2[i,1] = kp2[i]
    put2[i,2] = impput2[i]
    put2 [i,3] = mp2[i]
    put2[i,4] = BSP2[i]
    put2[i,5] = treeput2[i]
    put2[i,6] = Atreeput2[i]
    
    put3[i,1] = kp3[i]
    put3[i,2] = impput3[i]
    put3 [i,3] = mp3[i]
    put3[i,4] = BSP3[i]
    put3[i,5] = treeput3[i]
    put3[i,6] = Atreeput3[i]
    
   
}
write.csv(call1, file = "bcall1.csv")
write.csv(call2, file = "bcall2.csv")
write.csv(call3, file = "bcall3.csv")
write.csv(put1, file = "bput1.csv")
write.csv(put2, file = "bput2.csv")
write.csv(put3, file = "bput3.csv")
# Part d Absolute error ----
AbsoluteError <- function(S,k,tau,r= 0.0075,vol) {
    N <- c(10,20,30,40,50,100,150,200,250,300,350,400)
    err <- matrix(0,nrow = length(N), ncol = length(vol))
    length(N)
    
    for(i in 1:length(N)){
        for(j in 1:length(vol)){
            err[i,j] <- abs(Binomialtree(F,F,S,k[j],tau,r,vol[j],N[i]) - 
                                PutOption(S,k[j],tau,r,vol[j]));
            
        }
    }
    return(err)
}
err1 <- matrix(0,nrow = 12, ncol = 10)
err2 <- matrix(0,nrow = 12, ncol = 10)
err3 <- matrix(0,nrow = 12, ncol = 10)
err1 <- AbsoluteError(S,kp1,tau,r=0.0075,impput1) 
err2 = AbsoluteError(S,kp2,tau,r=0.0075,impput2)
err3 = AbsoluteError(S,kp3,tau,r=0.0075,impput3)
err1
write.csv(err1, file = "berr1.csv")
write.csv(err2, file = "berr2.csv")
write.csv(err3, file = "berr3.csv")

# Plotting error ----
N <- c(10,20,30,40,50,100,150,200,250,300,350,400)
plot(N,err1[,1],type = 'l', col = 'purple', ylab = 'error', xlab = "N")
for (i in 2:10) {
    lines(N,err1[,i], type = 'l' , col = 'purple')
}
err2
plot(N,err2[,1],type = 'l', col = 'red', ylab = 'error', xlab = "N")
for (i in 2:10) {
    plot(N,err2[,i] , type = 'l', col = 'red' , ylab = "Error", xlab = "N")
}
plot(N,err3[,1],type = 'l', col = 'purple', ylab = 'error', xlab = "N")
for (i in 2:10) {
    lines(N,err3[,i], type = 'l' , col = 'purple')
}
err3

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
    if(is.nan(price)){
        return(0)
    }
   return(price)
}
    
   
    TrinomialTree(F,T,10,10,0.3,0.06,0.25,350,0.03,0.2)
    PutOption(10,10,0.3,0.06,0.2)
   CallOption(10,10,0.3,0.06,0.2)
    S0 = k =  100 
    tau = 1
    sigma = 0.25
    r = 0.06
    N = 350
    div = 0.03
    dx = sigma * sqrt(3*tau/N)
    EUCall =  TrinomialTree(T,F, S0,k,tau,r,sigma,350,div,dx)
    EUCall
    EUPut = TrinomialTree(F,F, S0,k,tau,r,sigma,200,div,dx)
    EUPut
    ACall = TrinomialTree(T,T, S0,k,tau,r,sigma,200,div,dx)
    ACall
    APut = TrinomialTree(F,T, S0,k,tau,r,sigma,200,div,dx)
    APut

# Part b ---- 
#  Calculate dx   ----
    S0 = k =  100 
    tau = 1
    sigma = 0.25
    r = 0.06
    N = 350
    div = 0.03
    dx1c <- c()
    dx1p <- c()
    dx2c <- c()
    dx2p <- c()
    dx3c <- c()
    dx3p <- c()
    for (i in 1:10) {
        dx1c[i] = impcall1[i] * sqrt(3*(30/252)/N)
        dx1p[i] = impput1[i] * sqrt(3*(30/252)/N)
        dx3c[i] = impcall3[i] * sqrt(3*(90/252)/N)
        dx3p[i] = impput3[i] * sqrt(3*(90/252)/N)
        dx2c[i] = impcall2[i] * sqrt(3*(60/252)/N)
        dx2p[i] = impput2[i] * sqrt(3*(60/252)/N)
    }
   
############# European Options Trinomial ###########
# 1 month trinomial tree European ----
    
EUtricall1 <- c()
BSEUcall1 <- c()
for(i in 1:10) {
    EUtricall1[i] = TrinomialTree(T,F,100,100,30/252,0.06,impcall1[i],N,0.03,dx1c[i])
    BSEUcall1[i] = CallOption(100,100,30/252,0.06,impcall1[i])
}
EUtricall1
BSEUcall1
EUtriput1 <- c()
BSEUput1 <- c()
for(i in 1:10) {
    EUtriput1[i] = TrinomialTree(F,F,100,100,30/252,0.06,impput1[i],200,0.03,dx1p[i])
    BSEUput1[i] = PutOption(100,100,60/252,0.06,impput1[i])
}
EUtriput1
BSEUput1
    
# 2 month trinomial tree European ----
EUtricall2 <- c()
BSEUcall2 <- c()
for(i in 1:10) {
    EUtricall2[i] = TrinomialTree(T,F,100,100,60/252,0.06,impcall2[i],200,0.03,dx2c[i])
    BSEUcall2[i] = CallOption(100,100,60/252,0.06,impcall2[i])
}
EUtricall2
BSEUcall2
EUtriput2 <- c()
BSEUput2 <- c()
for(i in 1:10) {
    EUtriput2[i] = TrinomialTree(F,F,100,100,30/252,0.06,impput2[i],200,0.03,dx2p[i])
    BSEUput2[i] = PutOption(100,100, 30/252,0.06,impput2[i])
}
EUtriput2
BSEUput2

# 3 month trinomial tree European ----
EUtricall3 <- c()
BSEUcall3 <- c()
for(i in 1:10) {
    EUtricall3[i] = TrinomialTree(T,F,100,100,90/252,0.06,impcall3[i],N,0.03,dx3c[i])
    BSEUcall3[i] = CallOption(100,100,90/252,0.06,impcall3[i])
}
EUtricall3
BSEUcall3
EUtriput3 <- c()
BSEUput3 <- c()
for(i in 1:10) {
    EUtriput3[i] = TrinomialTree(F,F,100,100,90/252,0.06,impput3[i],200,0.03,dx3p[i])
    BSEUput3[i] = PutOption(100,100,90/252,0.06,impput3[i])
}
EUtriput3
BSEUput3

############ American Option Trinomial tree ############
# 1 month trinomial tree----
Atricall1 <- c()
BSAcall1 <- c()
for(i in 1:10) {
    Atricall1[i] = TrinomialTree(T,T,100,100,30/252,0.06,impcall1[i],N,0.03,dx1c[i])
    BSAcall1[i] = CallOption(100,100,30/252,0.06,impcall1[i])
}
Atricall1
BSAcall1
Atriput1 <- c()
BSAput1 <- c()
for(i in 1:10) {
    Atriput1[i] = TrinomialTree(F,T,100,100,30/252,0.06,impput1[i],200,0.03,dx1p[i])
    BSAput1[i] = PutOption(100,100,60/252,0.06,impput1[i])
}
Atriput1
BSAput1

# 2 month trinomial tree----
Atricall2 <- c()
BSAcall2 <- c()
for(i in 1:10) {
    Atricall2[i] = TrinomialTree(T,T,100,100,60/252,0.06,impcall2[i],200,0.03,dx2c[i])
    BSAcall2[i] = CallOption(100,100,60/252,0.06,impcall2[i])
}
Atricall2
BSAcall2
Atriput2 <- c()
BSAput2 <- c()
for(i in 1:10) {
    Atriput2[i] = TrinomialTree(F,T,100,100,30/252,0.06,impput2[i],200,0.03,dx2p[i])
    BSAput2[i] = PutOption(100,100, 30/252,0.06,impput2[i])
}
Atriput2
BSAput2

# 3 month trinomial tree ----
Atricall3 <- c()
BSAcall3 <- c()
for(i in 1:10) {
    Atricall3[i] = TrinomialTree(T,T,100,100,60/252,0.06,impcall3[i],N,0.03,dx3c[i])
    BSAcall3[i] = CallOption(100,100,60/252,0.06,impcall3[i])
}
Atricall3
BSAcall3
Atriput3 <- c()
BSAput3 <- c()
for(i in 1:10) {
    Atriput3[i] = TrinomialTree(F,T,100,100,90/252,0.06,impput3[i],200,0.03,dx3p[i])
    BSAput3[i] = PutOption(100,100,90/252,0.06,impput3[i])
}
dx3p
impput3
Atriput3
BSAput3
BSEUput3
tcall1 <- matrix(0,nrow = 10, ncol = 5)
tcall2 <- matrix(0,nrow = 10, ncol = 5)
tcall3 <- matrix(0,nrow = 10, ncol = 5)
tput1 <- matrix(0,nrow = 10, ncol = 5)
tput2 <- matrix(0,nrow = 10, ncol = 5)
tput3 <- matrix(0,nrow = 10, ncol = 5)

for(i in 1:10) {
    tcall1[i,1] = kc1[i]
    tcall1[i,2] = impcall1[i]
    tcall1[i,3] = BSEUcall1[i]
    tcall1[i,4] = EUtricall1[i]
    tcall1[i,5] = Atricall1[i]
    
    tcall2[i,1] = kc2[i]
    tcall2[i,2] = impcall2[i]
    tcall2[i,3] = BSEUcall2[i]
    tcall2[i,4] = EUtricall2[i]
    tcall2[i,5] = Atricall2[i]
    
    tcall3[i,1] = kc3[i]
    tcall3[i,2] = impcall3[i]
    tcall3[i,3] = BSEUcall3[i]
    tcall3[i,4] = EUtricall3[i]
    tcall3[i,5] = Atricall3[i]
    
    tput1[i,1] = kc1[i]
    tput1[i,2] = impput1[i]
    tput1[i,3] = BSEUput1[i]
    tput1[i,4] = EUtriput1[i]
    tput1[i,5] = Atriput1[i]
    
    tput2[i,1] = kc2[i]
    tput2[i,2] = impput2[i]
    tput2[i,3] = BSEUput2[i]
    tput2[i,4] = EUtriput2[i]
    tput2[i,5] = Atriput2[i]
    
    tput3[i,1] = kc3[i]
    tput3[i,2] = impput3[i]
    tput3[i,3] = BSEUput3[i]
    tput3[i,4] = EUtriput3[i]
    tput3[i,5] = Atriput3[i]
    
    
}
kc3
impput3
BSEUput3
EUtriput3
Atriput3
write.csv(tcall1, file = "tcall1.csv")
write.csv(tcall2, file = "tcall2.csv")
write.csv(tcall3, file = "tcall3.csv")
write.csv(tput1, file = "tput1.csv")
write.csv(tput2, file = "tput2.csv")
write.csv(tput3, file = "tput3.csv")

# Part d Absolute error----
TriAbsoluteError <- function(S,k,tau,r= 0.0075,vol, N, div, dx) {
    return(abs(TrinomialTree(F,T,S,k,tau,r,vol,N,div,dx) - 
                   PutOption(S,k,tau,r,vol)))
    
       
}
N <- c(10,20,30,40,50,100,150,200,250,300,350,400)
error1 <- matrix(0,nrow = length(N), ncol = 10)
error2 <- matrix(0,nrow = length(N), ncol = 10)
error3 <- matrix(0,nrow = length(N), ncol = 10)

for(i in 1:length(N)){
    for(j in 1:10){
        error1[i,j] = TriAbsoluteError(S,kp1[j],30/252,r=0.0075,impput1[j],N[i],div=0,dx1p[j])
        error2[i,j]=TriAbsoluteError(S,kp2[j],60/252,r=0.0075,impput2[j], N[i], div=0,dx2p[j])
  error3[i,j]=TriAbsoluteError(S,kp3[j],90/252,r = 0.0075, impput3[j], N[i], div=0, dx3p[j])
    }
}
error1
error2
error3
write.csv(error1, file = "terror1.csv")
write.csv(error2, file = "terror2.csv")
write.csv(error3, file = "terror3.csv")

N <- c(10,20,30,40,50,100,150,200,250,300,350,400)
plot(N,error1[,1],type = 'l', col = 'purple', ylab = 'error', xlab = "N")
for (i in 2:10) {
    lines(N,err1[,i], type = 'l' , col = 'purple')
}
err2
plot(N,err2[,1],type = 'l', col = 'red', ylab = 'error', xlab = "N")
for (i in 2:10) {
    plot(N,err2[,i] , type = 'l', col = 'red' , ylab = "Error", xlab = "N")
}
plot(N,err3[,1],type = 'l', col = 'purple', ylab = 'error', xlab = "N")
for (i in 2:10) {
    lines(N,err3[,i], type = 'l' , col = 'purple')
}
err3
# ############Question 3 Exotic Options##########
#part a) European up and out call option ----
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
    #print((S))
    return((C[1]))
} 
    
EUpandOutCall(10,10,0.3,0.01,0.2,0,11,5)
# part b) European up and out call formula ----
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

 (EUpandOutCall(10,10,0.3,0.01,0.2,0,11,500) -  EUpandOutFormula(10,10,0.3,0.01,0.2,11,0))




#Up and in formula ----
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

    return(C[1])
} 
EUpandInCall(10,10,0.3,0.01,0.2,0,11,20)
UpandInFormula(10,10,0.3,0.01,0.2,11,0)# mathces analytical at 20 steps

# In - Out Parity ----

N = 200
Parity = EUpandInCall(10,10,0.3,0.01,0.2,0,11,N) + EUpandOutCall(10,10,0.3,0.01,0.2,0,11,N) -
    CallOption(10,10,0.3,0.01,0.2)
abs(Parity)
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
        if(S[i] >= H){
            C[i] <- max(0,(k - S[i]))
        } else C[i] = 0
    }
    for(i in N:1) {
        for(j in 1:i){
            
            S[j] = S[j] / exp(dxd)
            #if(S[j] > H){
                C[j] = disc * (pu * C[j+1] + pd * C[j])
                C[j] = max(C[j], k - S[j])
            }# else
             #   C[j] = 0
    }
    C
    print(S)
    return(C)
}
AUpandInPut(10,10,0.3,0.01,0.2,0,11,5)
PutOption(10,10,0.3,0.01,0.2)


   AUpandInForm <- function(S,k,tau,r,sigma,div,H,N){
       part1 = ((S/H)^(1-(2*r))/sigma^2) *
           (Binomialtree(F,T,(H^2/S),k,tau,r,sigma,N) - Binomialtree(F,F,(H^2/S),k,tau,r,sigma,N))
       part2 = UpandInFormula(S,k,tau,r,sigma,H,div)
       ans = part2 + part1
       return(ans)
   }
   Binomialtree(F,T,12.1,10,0.3,0.01,0.2,10)
   AUpandInForm(10,10,0.3,0.01,0.2,0,11,200)

## up and in = put - up and out