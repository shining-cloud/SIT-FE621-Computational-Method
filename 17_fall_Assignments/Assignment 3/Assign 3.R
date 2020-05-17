############# Question 1#############
# Explicit finite difference ----

FiniteExplicit <- function(S0,k,r,tau,sig,N, Nj, div, dx, isCall) {
    
    dt = tau/N
    nu = r - div - 0.5 * sig^2
    edx = exp(dx)
    pu = 0.5 * dt * ( (sig/dx)^2 + nu/dx )
    pm = 1.0 - dt *   (sig/dx)^2 - r*dt 
    pd = 0.5 * dt * ( (sig/dx)^2 - nu/dx)
    
    cp = ifelse(isCall, 1, -1)
    
    
    S = c()
    S
    V = matrix(0, 2*Nj+1, N+1)
    
    V
    S[(2*Nj+1)] = S0 * exp(-Nj * dx)
    S
    
    
    for(i in (2*Nj) :1) {
        S[i] = S[i+1] * edx
        
    }
    S
    for (j in (2*Nj+1):1) {
        V[j, (N+1)] = max( 0, cp * (S[j] - k))
    }
    V
    
    for(i in N:1){
        for(j in (2*Nj):2){
            V[j,i] = pu*V[j-1,i+1] + pm*V[j,i+1] + pd*V[j+1,i+1] 
            
        }
        # Boundary conditions
        if(isCall) {
            bc = S[1] - S[2]
            V[(2*Nj+1),i] = V[(2*Nj), i] 
            V[1,i ] = V[2,i] + bc
        }
        else {
            bc = S[2] - S[2*Nj+1]
            V[(2*Nj+1),i] = V[(2*Nj), i] + bc
            V[1,i ] = V[2,i] 
        }
        
    }
    
    
    price = V[Nj+1,1]
    list(price = price, Stock = S, Value = V)
}

FiniteExplicit(100,100,0.06,1,0.2,3,3,0.03,0.2,T)
FiniteExplicit(100,100,0.06,1,0.2,3,3,0.03,0.2,F)



# Implicit finite difference
# Implicit finite difference----
FiniteImplicit <- function(S0,k,r,tau,sig,N, Nj, div, dx, isCall){
    
    dt = tau/N
    nu = r - div - 0.5 * sig^2
    edx = exp(dx)
    pu = -0.5 * dt * ( (sig/dx)^2 + nu/dx )
    pm = 1.0 + dt *   (sig/dx)^2 + r*dt 
    pd = -0.5 * dt * ( (sig/dx)^2 - nu/dx)
    
    cp = ifelse(isCall, 1, -1)
    
    
    S = c()
    
    V = matrix(0, 2*Nj+1, N+1)
    
    pp = matrix(0, 2*Nj+1, N+1)
    pmp = matrix(0, 2*Nj+1, N+1)
    
    pp
    pmp
    S[(2*Nj+1)] = S0 * exp(-Nj * dx)
    
    
    
    for(i in (2*Nj) :1) {
        S[i] = S[i+1] * edx
        
    }
    S
    for (j in (2*Nj+1):1) {
        V[j, (N+1)] = max( 0, cp * (S[j] - k))
    }
    V
    lambdaL = -1 * (S[2*Nj] - S[2*Nj+1])
    lambdaU = 0
    
    ImplicitTridiagnol <- function(V, pu, pm, pd, lambdaL, lambdaU, colI){
        pp = numeric(2*Nj+1)
        pmp = numeric(2*Nj+1)
        pmp[(2*Nj)] = pm + pd
        pp[(2*Nj)] = V[(2*Nj),(N+1)] + pd*lambdaL
        pmp
        pp
        #  colI = 3
        for(i in (2*Nj-1): 2){
            pmp[i] = pm - pu * pd /pmp[i+1]
            pp[i] = V[i,(colI+1)] - pp[i+1]*pd/pmp[i+1]
        }
        pmp
        pp
        V[1,colI] = (pp[2] + pmp[2] * lambdaU)/(pu + pmp[2])
        V[2,colI] = V[1,colI] - lambdaU
        V
        #i = 3
        for(i in 3:(2*Nj)){
            V[i,colI] = (pp[i] - pu*V[i-1,colI])/pmp[i]
        }
        V
        V[2*Nj+1,colI] = V[2*Nj,colI] - lambdaL
        V
        list(V=V, pmp = pmp, pp = pp)
    }
    
    for(i in N:1){
        tri = ImplicitTridiagnol(V,pu,pm,pd,lambdaL,lambdaU,i)
        V = tri$V
        pmp[,i] = tri$pmp
        pp[,i] = tri$pp
        
    }
    V
    price = V[Nj+1,1]
    list(price = price, Stock = S, Value = V)
}
FiniteImplicit(100,100,0.06,1,0.2,3,3,0.03,0.2,T)
FiniteImplicit(100,100,0.06,1,0.2,3,3,0.03,0.2,F)


# Part C  ----
parameters <- function(tau,sig,eps){
    delt = eps/(3*sig^2 + 1)
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    list(N = N, Delta_T = delt, Delta_x = delx, Nj = Nj)
}

para = parameters(1,0.25,0.001)
para
k=100
tau=1
S0=100
r=0.06
sig=0.25
N=para$N
div=0.03
dx=para$Delta_x
K = k
Tm = tau
# Part D ----

para = parameters(1,0.25,0.001)

FiniteExplicit(100, 100, 0.06, 1, 0.25, para$N, para$Nj, 0.03, para$Delta_x, T)  # 10.85769
FiniteExplicit(100, 100, 0.06, 1, 0.25, para$N, para$Nj, 0.03, para$Delta_x, F)  # 11.08906


FiniteImplicit(100, 100, 0.06, 1, 0.25, para$N, para$Nj, 0.03, para$Delta_x, T)  # 10.82705
FiniteImplicit(100, 100, 0.06, 1, 0.25, para$N, para$Nj, 0.03, para$Delta_x, F)  # 8.235997

# Part E----
S0 = k = 100
r = 0.06
tau = 1
sig = 0.25
div = 0.06
CallOption <- function(S,k, tau, r, sigma,q) {
    d1 <- (log(S/k) + (r - q + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price =  S* exp(-q*tau)*pnorm(d1) - k * exp(-r*tau)*pnorm(d2)
    return(price)
    
    
}
PutOption <- function(S,k, tau, r, sigma,q){
    
    d1 <- (log(S/k) + (r - q + sigma^2/2 ) * tau) / (sigma * sqrt(tau))
    d2 <- d1 - sigma*sqrt(tau)
    price <- k * exp(-r*tau) * pnorm(-d2) - S * exp(-q*tau)*pnorm(-d1)
    return(price)
}
BSCall = CallOption(100,100,1,0.06,0.25,0.03)
BSPut = PutOption(100,100,1,0.06,0.25,0.03)
eps = 0.001
delt = eps/(3*sig^2 + 1)
N = tau/delt
delx = sig * sqrt(3*delt)
Nj = (2*sqrt(3*N) - 1)/2


counter = 0 
while(abs(FiniteExplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,T)$price - BSCall) > eps){
    delt = delt + 0.00001
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    counter = counter + 1
}
N
Nj
delt
delx
counter # 90
BSCall
er = abs(FiniteExplicit(100,100,0.06,1,0.25,574.0181,40.99764,0.03,0.018077332,T)$price - BSCall)
FiniteExplicit(100,100,0.06,1,0.25,574,40,0.03,0.01807332,T) # 11.00988

counter = 0 
eps = 0.001
delt = eps/(3*sig^2 + 1)  
N = tau/delt
delx = sig * sqrt(3*delt)
Nj = (2*sqrt(3*N) - 1)/2
while(abs(FiniteExplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,F)$price - BSPut) > eps){
    delt = delt + 0.00001
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    counter = counter + 1
}
N
Nj
delt
delx
counter
counter # 1141
FiniteExplicit(100,100,0.06,1,0.25,1187,60,0.03,0.125,F)
FiniteExplicit(100,100,0.06,1,0.25,66.88021,13.6641,0.03,delx,F)
er1 = FiniteExplicit(100,100,0.06,1,0.25,66.88201,13.6641,0.03,0.05294,F) 
er1$price - BSPut

counter = 0 
eps = 0.001
delt = eps/(3*sig^2 + 1)
N = tau/delt
delx = sig * sqrt(3*delt)
Nj = (2*sqrt(3*N) - 1)/2
while(abs(FiniteImplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,T)$price - BSCall) > eps){
    delt = delt + 0.000001
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    counter = counter + 1
}
N
Nj
delt
delx
counter
er = abs(FiniteImplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,T)$price - BSCall)

counter # 201
counter = 0 
eps = 0.001
delt = eps/(3*sig^2 + 1) 
N = tau/delt
delx = sig * sqrt(3*delt)
Nj = (2*sqrt(3*N) - 1)/2
while(abs(FiniteImplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,F)$price - BSPut) > eps){
    delt = delt + 0.000001
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    counter = counter + 1
}
N
Nj
delt
delx
counter
er = abs(FiniteImplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,F)$price - BSPut)
er
counter # 458

counter = 0 
eps = 0.001
delt = eps/(3*sig^2 + 1) 
N = tau/delt
delx = sig * sqrt(3*delt)
Nj = (2*sqrt(3*N) - 1)/2
while(abs(CrankNicolson(100,100,0.06,1,0.25,N,Nj,0.03,delx,T)$price - BSCall) > eps){
    delt = delt + 0.00001
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    counter = counter + 1
}
N
Nj
delt
delx
counter # counter = 2216
er = abs(CrankNicolson(100,100,0.06,1,0.25,N,Nj,0.03,delx,T)$price - BSCall)
er
counter = 0 
eps = 0.001
delt = eps/(3*sig^2 + 1) # in terms of tau 
N = tau/delt
delx = sig * sqrt(3*delt)
Nj = (2*sqrt(3*N) - 1)/2

while(abs(CrankNicolson(100,100,0.06,1,0.25,N,Nj,0.03,delx,F)$price - BSPut) > eps){
    delt = delt + 0.00001
    N = tau/delt
    delx = sig * sqrt(3*delt)
    Nj = (2*sqrt(3*N) - 1)/2
    counter = counter + 1
}
N
Nj
delt
delx
counter # counter = 2216
er = abs(FiniteImplicit(100,100,0.06,1,0.25,N,Nj,0.03,delx,F)$price - BSPut)
er
#Part F----
sigma = seq(from = 0.05, to = 0.6, by = 0.05)
sigma
# keeping N,dx,and dt constant
 dt <- para$Delta_T
 dx <- para$Delta_x
 
 
 tau = 1
eps = 0.001
 r = 0.06
 div = 0.03
 pu = pd = pm = c()
 
N = tau/dt
 for(i in 1:length(sigma)){
 pu[i] = -0.5 * dt * ( (sigma[i]/dx)^2 + (r - div - 0.5 * sigma[i]^2)/dx )
 pm[i] =  1.0 + dt *   (sigma[i]/dx)^2 + r*dt 
 pd[i]= -0.5 * dt * ( (sigma[i]/dx)^2 - (r - div - 0.5 * sigma[i]^2)/dx)
 }
pu
pm    
pd
prob = matrix(0,nrow = 12, ncol = 3 )
    for(i in 1:12){
            prob[i,1] = pu[i]
            prob[i,2] = pm[i]
            prob[i,3] = pd[i]
    }
    prob
    write.csv(prob, file = "Part F.csv")
    plot(sigma,pm,type = 'l' , col = 'purple',  xlab = 'Volitility', xlim =c(0,0.6), ylim = c(-2.5,3.0),
        ylab ='Pu , Pm, Pd' )
    lines(sigma,pd,col = 'green')
    lines(sigma,pu, col = 'yellow')
   
    
# Crank-Nicolson finite difference----
    CrankNicolson <- function(S0,k,r,tau,sig,N,Nj,div,dx,isCall){
        
        dt = tau/N
        nu = r - div - 0.5 * sig^2
        edx = exp(dx)
        pu = -0.25 * dt * ( (sig/dx)^2 + nu/dx )
        pm = 1.0 + 0.5 * dt *   (sig/dx)^2 + 0.5 * r * dt 
        pd = -0.25 * dt * ( (sig/dx)^2 - nu/dx)
        
        cp = ifelse(isCall, 1, -1)
        S = c()
        
        V = matrix(0, 2*Nj+1, N+1)
        
        pp = matrix(0, 2*Nj+1, N+1)
        pmp = matrix(0, 2*Nj+1, N+1)
        
        pp
        pmp
        S[(2*Nj+1)] = S0 * exp(-Nj * dx)
        
        
        
        for(i in (2*Nj) :1) {
            S[i] = S[i+1] * edx
            
        }
        S
        for (j in (2*Nj+1):1) {
            V[j, (N+1)] = max( 0, cp * (S[j] - k))
        }
        V
        lambdaL = -1 * (S[2*Nj] - S[2*Nj+1])
        lambdaU = 0
        
        CrankNicolsontridagonal <- function(V, pu, pm, pd, lambdaL, lambdaU, colI){
            pp = numeric(2*Nj+1)
            pmp = numeric(2*Nj+1)
            pmp[(2*Nj)] = pm + pd
            pp[(2*Nj)] = -pu * V[2*Nj-1,(N+1)] - (pm - 2) * V[2*Nj,N+1] - pd*V[2*Nj+1, N+1]+  pd*lambdaL
            pmp
            pp
            #  colI = 3
            for(j in (2*Nj-1): 2){
                pmp[j] = pm - pu * pd /pmp[j+1]
                pp[j] = ( - pu   *V[j-1, colI+1] 
                          -(pm-2) *V[j  , colI+1]
                          - pd    *V[j+1, colI+1] 
                          -pp[j+1]*pd/pmp[j+1])
            }
            pmp
            pp
            V[1,colI] = (pp[2] + pmp[2] * lambdaU)/(pu + pmp[2])
            V[2,colI] = V[1,colI] - lambdaU
            V
            # i = 3
            for(i in 3:(2*Nj)){
                V[i,colI] = (pp[i] - pu*V[i-1,colI])/pmp[i]
            }
            V
            V[2*Nj+1,colI] = V[2*Nj,colI] - lambdaL
            V
            list(V=V, pmp = pmp, pp = pp)
        }
        for(i in N:1){
            tri = CrankNicolsontridagonal(V,pu,pm,pd,lambdaL,lambdaU,i)
            V = tri$V
            pmp[,i] = tri$pmp
            pp[,i] = tri$pp
            
        }
        V
        price = V[Nj+1,1]
        list(price = price)
    }
    CrankNicolson(100,100,0.06,1,0.2,3,3,0.03,0.2,T)
    CrankNicolson(100,100,0.06,1,0.25,para$N,para$Nj,0.03,para$Delta_x,T) #10.8286
    CrankNicolson(100,100,0.06,1,0.25,N,Nj,0.03,delx,F) # 8.237
    CrankNicolson(100,100,0.06,1,0.25,43.47428,10.92028,0.03,0.06567623,T)
    FiniteImplicit(100,100,0.06,1,0.25,331,31,0.03,0.02380036,F)
    FiniteExplicit(100,100,0.06,1,0.25,574.018,40.99764,0.03,0.01807,T)
 
# Hedge sensitivities ----
ans <- Explicitfindiff(100,100,0.06,1,0.2,3,0.03,0.2,T)

   HedgeSentivities <- function(S0,k,r,tau,sigma, N, Nj, div, dx , isCall){
      
   
    dt = tau/N
    expdiff <- FiniteExplicit(S0,k,r,tau,sigma,N,Nj, div, dx, isCall)
   
     del <- (expdiff$Value[N,1] - expdiff$Value[N+2,1]) / 
         (expdiff$Stock[N] - expdiff$Stock[N+2])
     del
     
     gamma1 <- (expdiff$Value[N,1] - expdiff$Value[N-1,1])/
                (expdiff$Stock[N] - expdiff$Stock[N-1])
     gamma2 <- (expdiff$Value[N-1,1] - expdiff$Value[N+2,1])/
         (expdiff$Stock[N-1] - expdiff$Stock[N+2])
     gamma3 = 0.5*(expdiff$Stock[N] - expdiff$Stock[N+2])
     gamma = (gamma1 - gamma2) / gamma3
     gamma
     expdiff$Value
    theta = (expdiff$Value[N+1,1] - expdiff$Value[N+1,2]) / dt
    theta
    delsig = 0.001*sigma
    vega = (FiniteExplicit(S0,k,r,tau,(sigma + delsig),N,Nj, div,dx,isCall)$price
    - FiniteExplicit(S0,k,r,tau,(sigma - delsig),N,Nj, div,dx,isCall)$price) / (2*delsig)
    vega
    list(Delta = del, Gamma = gamma, Theta = theta, Vega = vega)
   }
        
        
   HedgeSentivities(100,100,0.06,1,0.2,3,3,0.03,0.2,T)
   S0 = k =100
   r = 0.06
   tau = 1
   sigma = 0.2
   N = Nj = 3
   div = 0.03
   dx = 0.2
   isCall = T
   
ans
expdiff1 <- Explicitfindiff(100,100,0.06,1,0.2,3,0.03,0.2,T)
expdiff2 <- Explicitfindiff(95,100,0.06,1,0.2,3,0.03,0.2,T)
expdiff1$price - expdiff2$price
############## Question 2##########
# Part a) Implied Vol ----
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
mc1
mp1
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
S = 49.5

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
for(i in 1:10){
    impcall1[i] =   Secant(S,kc1[i],30/252,0.0075,mc1[i],T)
}
impcall1
S
kc1
mc1
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
for(i in 1:10){
    impcall3[i] =   Secant(S,kc3[i],90/252,0.0075,mc3[i],T)
    if(impcall3[i] < 0 )
        impcall3[i] = 0
}

impcall3
impput1 <- c()
impput2 <- c()
impput3 <- c()
for(i in 1:10){
    impput1[i] =   Secant(S,kp1[i],30/252,0.0075,mp1[i],F)
    if(impput1[i] < 0 )
        impput1[i] = 0
}
impput1
for(i in 1:10){
    impput2[i] =   Secant(S,kp2[i],60/252,0.0075,mp2[i],F)
    if(impput2[i] < 0 )
        impput2[i] = 0
}
impput2
kp3
for(i in 1:10){
    impput3[i] =   Secant(S,kp3[i],90/252,0.0075,mp3[i],F)
    if(impput3[i] < 0 ) {
        impput3[i] = 0
    }
}
impput3
#Part b) Convergence----
Convergence <- function(tau,eps = 0.001, isCall){
    if(isCall){
    if(tau == 30/252){
        
        for(i in 1:10){
            sig[i] = impcall1[i]
        }
    } 
    if (tau == 60/252){
        
        for(i in 1:10){
            sig[i] = impcall2[i]
        }
    }
    if(tau == 90/252) {
        
        for(i in 1:10){
            sig[i] = impcall3[i]
        }
    }
    }
    else{if(tau == 30/252){
        
        for(i in 1:10){
            sig[i] = impput1[i]
        }
    } 
        if (tau == 60/252){
            
            for(i in 1:10){
                sig[i] = impput2[i]
            }
        }
        if(tau == 90/252) {
            
            for(i in 1:10){
                sig[i] = impput3[i]
            }
        }
        
    }
dxc = dtc = Nc =  Njc = c()
for(i in 1:10){
    dtc[i] = eps/(3* sig[i]^2 + 1)
    Nc[i] = tau/dtc[i]
    dxc[i] = sig[i] * sqrt(3*dtc[i])
    Njc[i] = (2*sqrt(3*Nc[i]) - 1)/2
}

list( dt = dtc, N = Nc, dx = dxc, Nj = Njc)
}

# 1 month ----
    ConCall1 = Convergence(30/252,eps,T)
    ConPut1 = Convergence(30/252,eps,F)
IC1 = EC1 = CNC1 = c()
IP1 = EP1 = CNP1 = c()
for(i in 1:10) {
   EC1[i] =  FiniteExplicit(S,kc1[i], 0.0075, 30/252, impcall1[i], ConCall1$N[i], ConCall1$Nj[i], 0, ConCall1$dx[i], T )
  IC1[i] = FiniteImplicit(S,kc1[i], 0.0075, 30/252, impcall1[i], ConCall1$N[i], ConCall1$Nj[i], 0 , ConCall1$dx[i], T)
   CNC1[i] = CrankNicolson(S,kc1[i],0.0075,30/252,impcall1[i],ConCall1$N[i], ConCall1$Nj[i],0,ConCall1$dx[i], T)
   
   EP1[i] =  FiniteExplicit(S,kp1[i], 0.0075, 30/252, impput1[i], ConPut1$N[i], ConPut1$Nj[i], 0, ConPut1$dx[i], F )
   IP1[i] = FiniteImplicit(S,kp1[i], 0.0075, 30/252, impput1[i], ConPut1$N[i], ConPut1$Nj[i], 0 , ConPut1$dx[i], F)
   CNP1[i] = CrankNicolson(S,kp1[i],0.0075,30/252,impput1[i],ConPut1$N[i], ConPut1$Nj[i],0,ConPut1$dx[i], F)
}

write.csv(impcall1, file = "vc1.csv")
write.csv(impput1, file = "vp1.csv")
write.csv(ac1, file = "ac1.csv")
write.csv(ap1, file = "ap1.csv")
write.csv(bc1, file = "bc1.csv")
write.csv(bp1, file = "bp1.csv")
write.csv(kc1, file = "kc1.csv")
write.csv(kp1, file = "kp1.csv")
write.csv(mc1, file = "mc1.csv")
write.csv(mp1, file = "mp1.csv")
write.csv(t(t(EC1)), file = 'EC1.csv')
write.csv(t(t(IC1)), file = 'IC1.csv')
write.csv(t(t(CNC1)), file = 'CNC1.csv')

write.csv(t(t(EP1)), file = 'EP1.csv')
write.csv(t(t(IP1)), file = 'IP1.csv')
write.csv(t(t(CNP1)), file = 'CNP1.csv')
# 2 month ----
ConCall2 = Convergence(60/252,eps,T)
ConPut2 = Convergence(60/252,eps,F)
IC2 = EC2 = CNC2 = c()
IP2 = EP2 = CNP2 = c()
for(i in 1:10) {
    EC2[i] =  FiniteExplicit(S,kc2[i], 0.0075, 60/252, impcall2[i], ConCall2$N[i], ConCall2$Nj[i], 0, ConCall2$dx[i], T )
    IC2[i] = FiniteImplicit(S,kc2[i], 0.0075, 60/252, impcall2[i], ConCall2$N[i], ConCall2$Nj[i], 0 , ConCall2$dx[i], T)
    CNC2[i] = CrankNicolson(S,kc2[i],0.0075,60/252,impcall2[i],ConCall2$N[i], ConCall2$Nj[i],0, ConCall2$dx[i], T)
    
    EP2[i] =  FiniteExplicit(S,kp2[i], 0.0075, 60/252, impput2[i], ConPut2$N[i], ConPut2$Nj[i], 0, ConPut2$dx[i], F)
    IP2[i] = FiniteImplicit(S,kp2[i], 0.0075, 60/252, impput2[i], ConPut2$N[i], ConPut2$Nj[i], 0, ConPut2$dx[i], F)
    CNP2[i] = CrankNicolson(S,kp2[i],0.0075,60/252,impput2[i], ConPut2$N[i], ConPut2$Nj[i], 0, ConPut2$dx[i], F)
}
write.csv(impcall2, file = "vc2.csv")
write.csv(impput2, file = "vp2.csv")
write.csv(ac2, file = "ac2.csv")
write.csv(ap2, file = "ap2.csv")
write.csv(bc2, file = "bc2.csv")
write.csv(bp2, file = "bp2.csv")
write.csv(kc2, file = "kc2.csv")
write.csv(kp2, file = "kp2.csv")
write.csv(mc2, file = "mc2.csv")
write.csv(mp2, file = "mp2.csv")

write.csv(t(t(EC2)), file = 'EC2.csv')
write.csv(t(t(IC2)), file = 'IC2.csv')
write.csv(t(t(CNC2)), file = 'CNC2.csv')

write.csv(t(t(EP2)), file = 'EP2.csv')
write.csv(t(t(IP2)), file = 'IP2.csv')
write.csv(t(t(CNP2)), file = 'CNP2.csv')

# 3 month ----
ConCall3 = Convergence(90/252,eps,T)
ConPut3 = Convergence(90/252,eps,F)
IC3 = EC3 = CNC3 = c()
IP3 = EP3 = CNP3 = c()
for(i in 1:10) {
    EC3[i] =  FiniteExplicit(S,kc3[i], 0.0075, 90/252, impcall3[i], ConCall3$N[i], ConCall3$Nj[i], 0, ConCall3$dx[i], T )
    IC3[i] = FiniteImplicit(S,kc3[i], 0.0075, 90/252, impcall3[i], ConCall3$N[i], ConCall3$Nj[i], 0 , ConCall3$dx[i], T)
    CNC3[i] = CrankNicolson(S,kc3[i],0.0075,90/252,impcall3[i],ConCall3$N[i], ConCall3$Nj[i],0, ConCall3$dx[i], T)
    
    EP3[i] =  FiniteExplicit(S,kp3[i], 0.0075, 90/252, impput3[i], ConPut3$N[i], ConPut3$Nj[i], 0, ConPut3$dx[i], F)
    IP3[i] = FiniteImplicit(S,kp3[i], 0.0075, 90/252, impput3[i], ConPut3$N[i], ConPut3$Nj[i], 0, ConPut3$dx[i], F)
    CNP3[i] = CrankNicolson(S,kp3[i],0.0075,90/252,impput3[i], ConPut3$N[i], ConPut3$Nj[i], 0, ConPut3$dx[i], F)
}

write.csv(impcall3, file = "vc3.csv")
write.csv(impput3, file = "vp3.csv")
write.csv(ac3, file = "ac3.csv")
write.csv(ap3, file = "ap3.csv")
write.csv(bc3, file = "bc3.csv")
write.csv(bp3, file = "bp3.csv")
write.csv(kc3, file = "kc3.csv")
write.csv(kp3, file = "kp3.csv")
write.csv(mc3, file = "mc3.csv")
write.csv(mp3, file = "mp3.csv")

write.csv(t(t(EC3)), file = 'EC3.csv')
write.csv(t(t(IC3)), file = 'IC3.csv')
write.csv(t(t(CNC3)), file = 'CNC3.csv')

write.csv(t(t(EP3)), file = 'EP3.csv')
write.csv(t(t(IP3)), file = 'IP3.csv')
write.csv(t(t(CNP3)), file = 'CNP3.csv')

write.csv(ConCall1, file = 'cc1.csv')
write.csv(ConPut1, file = 'cp1.csv')

write.csv(ConCall2, file = 'cc2.csv')
write.csv(ConPut2, file = 'cp2.csv')

write.csv(ConCall3, file = 'cc3.csv')
write.csv(ConPut3, file = 'cp3.csv')
#Part C Hedge Sensitivites ----
d1c = g1c = v1c = t1c= c()
d2c = g2c = v2c = t2c= c()
d3c = g3c = v3c = t3c= c()

d1p = g1p = v1p = t1p= c()
d2p = g2p = v2p = t2p= c()
d3p = g3p = v3p = t3p= c()
i = 1
for(i in 1:10){
    d1c[i] = HedgeSentivities(49.5,kc1[i],0.0075,30/252,impcall1[i],(ConCall1$N[i]), ConCall1$N[i],0,ConCall1$dx[i],T)$Delta
    g1c[i] = HedgeSentivities(49.5,kc1[i],0.0075,30/252,impcall1[i],(ConCall1$N[i]), ConCall1$N[i],0,ConCall1$dx[i],T)$Gamma
    v1c[i] = HedgeSentivities(49.5,kc1[i],0.0075,30/252,impcall1[i],(ConCall1$N[i]), ConCall1$N[i],0,ConCall1$dx[i],T)$Vega
    t1c[i] = HedgeSentivities(49.5,kc1[i],0.0075,30/252,impcall1[i],(ConCall1$N[i]), ConCall1$N[i],0,ConCall1$dx[i],T)$Theta
   
    d2c[i] = HedgeSentivities(49.5,kc2[i],0.0075,60/252,impcall2[i],(ConCall2$N[i]), ConCall2$N[i],0,ConCall2$dx[i],T)$Delta
    g2c[i] = HedgeSentivities(49.5,kc2[i],0.0075,60/252,impcall2[i],(ConCall2$N[i]), ConCall2$N[i],0,ConCall2$dx[i],T)$Gamma
    v2c[i] = HedgeSentivities(49.5,kc2[i],0.0075,60/252,impcall2[i],(ConCall2$N[i]), ConCall2$N[i],0,ConCall2$dx[i],T)$Vega
    t2c[i] = HedgeSentivities(49.5,kc2[i],0.0075,60/252,impcall2[i],(ConCall2$N[i]), ConCall2$N[i],0,ConCall2$dx[i],T)$Theta
    
    d3c[i] = HedgeSentivities(49.5,kc3[i],0.0075,90/252,impcall3[i],(ConCall3$N[i]), ConCall3$N[i],0,ConCall3$dx[i],T)$Delta
    g3c[i] = HedgeSentivities(49.5,kc3[i],0.0075,90/252,impcall3[i],(ConCall3$N[i]), ConCall3$N[i],0,ConCall3$dx[i],T)$Gamma
    v3c[i] = HedgeSentivities(49.5,kc3[i],0.0075,90/252,impcall3[i],(ConCall3$N[i]), ConCall3$N[i],0,ConCall3$dx[i],T)$Vega
    t3c[i] = HedgeSentivities(49.5,kc3[i],0.0075,90/252,impcall3[i],(ConCall3$N[i]), ConCall3$N[i],0,ConCall3$dx[i],T)$Theta
    
     }
write.csv(d1c,file = 'd1c.csv')
write.csv(g1c, file = 'g1c.csv')
write.csv(v1c, file = 'v1c.csv')
write.csv(t1c, file = 't1c.csv')

write.csv(d2c,file = 'd2c.csv')
write.csv(g2c, file = 'g2c.csv')
write.csv(v2c, file = 'v2c.csv')
write.csv(t2c, file = 't2c.csv')

write.csv(d3c,file = 'd3c.csv')
write.csv(g3c, file = 'g3c.csv')
write.csv(v3c, file = 'v3c.csv')
write.csv(t3c, file = 't3c.csv')

for(i in 1:10){
    d1p[i] = HedgeSentivities(49.5,kp1[i],0.0075,30/252,impput1[i],(ConPut1$N[i]), ConPut1$N[i],0,ConPut1$dx[i],F)$Delta
    g1p[i] = HedgeSentivities(49.5,kp1[i],0.0075,30/252,impput1[i],(ConPut1$N[i]), ConPut1$N[i],0,ConPut1$dx[i],F)$Gamma
    v1p[i] = HedgeSentivities(49.5,kp1[i],0.0075,30/252,impput1[i],(ConPut1$N[i]), ConPut1$N[i],0,ConPut1$dx[i],F)$Vega
    t1p[i] = HedgeSentivities(49.5,kp1[i],0.0075,30/252,impput1[i],(ConPut1$N[i]), ConPut1$N[i],0,ConPut1$dx[i],F)$Theta
    
    d2p[i] = HedgeSentivities(49.5,kp2[i],0.0075,60/252,impput2[i],(ConPut2$N[i]), ConPut2$N[i],0,ConPut2$dx[i],F)$Delta
    g2p[i] = HedgeSentivities(49.5,kp2[i],0.0075,60/252,impput2[i],(ConPut2$N[i]), ConPut2$N[i],0,ConPut2$dx[i],F)$Gamma
    v2p[i] = HedgeSentivities(49.5,kp2[i],0.0075,60/252,impput2[i],(ConPut2$N[i]), ConPut2$N[i],0,ConPut2$dx[i],F)$Vega
    t2p[i] = HedgeSentivities(49.5,kp2[i],0.0075,60/252,impput2[i],(ConPut2$N[i]), ConPut2$N[i],0,ConPut2$dx[i],F)$Theta
    
    d3p[i] = HedgeSentivities(49.5,kp3[i],0.0075,90/252,impput3[i],(ConPut3$N[i]), ConPut3$N[i],0,ConPut3$dx[i],F)$Delta
    g3p[i] = HedgeSentivities(49.5,kp3[i],0.0075,90/252,impput3[i],(ConPut3$N[i]), ConPut3$N[i],0,ConPut3$dx[i],F)$Gamma
    v3p[i] = HedgeSentivities(49.5,kp3[i],0.0075,90/252,impput3[i],(ConPut3$N[i]), ConPut3$N[i],0,ConPut3$dx[i],F)$Vega
    t3p[i] = HedgeSentivities(49.5,kp3[i],0.0075,90/252,impput3[i],(ConPut3$N[i]), ConPut3$N[i],0,ConPut3$dx[i],F)$Theta
    
}


write.csv(d1p,file = 'd1p.csv')
write.csv(g1p, file = 'g1p.csv')
write.csv(v1p, file = 'v1p.csv')
write.csv(t1p, file = 't1p.csv')

write.csv(d2p,file = 'd2p.csv')
write.csv(g2p, file = 'g2p.csv')
write.csv(v2p, file = 'v2p.csv')
write.csv(t2p, file = 't2p.csv')

write.csv(d3p,file = 'd3p.csv')
write.csv(g3p, file = 'g3p.csv')
write.csv(v3p, file = 'v3p.csv')
write.csv(t3p, file = 't3p.csv')

Call1
FiniteExplicit(49.5,52,0.0075,30/252,0.12667,124,18,0,0.00667,T)
HedgeSentivities(49.5,52,0.0075,30/252,0.12667,124,124,0,0.00667,T)
Call1del
# Question 3 ----
HestonExplicit <- function(param, S,k,r,tau,div,V){
    kappa = param[1]
    theta = param[2]
    sigma = param[3]
    V0 = param [4]
    rho =  param[5]
    lambda = param[6]
    
    Ns = Nv = 40
    Nt = 1000
    
    Smin = S[1]
    Vmin = V[1]
    Tmin = tau[1]
    Smax = S[Ns]
    Vmax = V[Nv]
    Tmax = tau[Nt]
    
    ds = (Smax - Smin)/ (Ns - 1)
    dv = (Smax - Smin)/ (Nv - 1)
    dt = (Tmax - Tmin) / (Nt -1)
    
    U = matrix(0 , nrow = Ns, ncol = Nv)
    u = matrix(0, nrow = Ns, ncol = Nv )
    
    for(i in 1: Ns) {
        for(j in 1:Nv){
            U[i,j] = max(0 , k - S[i])
        }
    }
    for (tm in 1:(Nt-1)){
        for(j in 1:(Nv-1)){
            U[1,j] = 0
            U[Ns,j] = max(0,  k -Smax )
        }
        for(s in 1:Ns){
            U[s,Nv] = max(0, k - S[s])
        }
        u = U
        
        for(s in 2:(Ns-1)){
            derV = (u[s,2] - u[s,1])/dv
            derS = (u[s,2] - u[s-1,1])/2/ds
            U[s,1] = u[s,1]*(1 - r*dt - kappa*theta*dt/dv) +
                dt*(0.5*(r-div)*s*u[s+1,1] - u[s-1,1]) + 
                kappa*theta*dt/dv*u[s,2];
        }
        u = U
        for (s in 2:(NS-1)) {
            for (v in 2:(NV-1)) {
                A = (1 - dt*(s-1)^2*(v-1)*dv - sigma^2*(v-1)*dt/dv - r*dt);
                B = (1/2*dt*(s-1)^2*(v-1)*dv - 1/2*dt*(r-div)*(s-1));
                CC = (1/2*dt*(s-1)^2*(v-1)*dv + 1/2*dt*(r-div)*(s-1));
                D = (1/2*dt*sigma^2*(v-1)/dv - 1/2*dt*kappa*(theta-(v-1)*dv)/dv);
                E = (1/2*dt*sigma^2*(v-1)/dv + 1/2*dt*kappa*(theta-(v-1)*dv)/dv);
                FF = 1/4*dt*sigma*(s-1)*(v-1);
                U[s,v] = A*u[s,v] +   B*u[s-1,v] + CC*u[s+1,v] +
                    D*u[s,v-1] + E*u[s,v+1] +
                    FF*(u[s+1,v+1]+u[s-1,v-1]-u[s-1,v+1]-u[s+1,v-1])
            }
        }
    }
    price = U[N+1,1]
    
    
    
    }
