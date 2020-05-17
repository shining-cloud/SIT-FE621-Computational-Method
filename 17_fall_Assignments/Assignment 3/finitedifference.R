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
    list(price = price)
}
FiniteImplicit <- function(S0,k,r,tau,sig,N, Nj, div, dx, isCall){
    
   isCall = T
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
        i = 3
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
    list(price = price)
}
k=100
tau=1
S0=100
r=0.06
sig=0.25
N=3
Nj = 3
div=0.03
dx=0.2
K = k
Tm = tau
FiniteExplicit(100,100,0.06,1,0.2,3,3,0.03,0.2,T)
FiniteExplicit(100,100,0.06,1,0.25,para$N,para$Nj,0.03,para$Delta_x,T)
FiniteImplicit(100,100,0.06,1,0.25,para$N,para$Nj,0.03,para$Delta_x,T)
FiniteImplicit(100,100,0.06,1,0.2,3,3,0.03,0.2,T)

Crank <- function(S0,k,r,tau,sig,N,Nj,div,dx,isCall){
    
    dt = tau/N
    nu = r - div - 0.5 * sig^2
    edx = exp(dx)
    pu = -0.25 * dt * ( (sig/dx)^2 + nu/dx )
    pm = 1.0 + 0.5 * dt *   (sig/dx)^2 + 0.5 * r * dt 
    pd = -0.25 * dt * ( (sig/dx)^2 - nu/dx)
    
    cp = ifelse(isCall, 1, -1)
    
    
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

Crank(100,100,0.06,1,0.2,30,30,0.03,0.2,T)
Crank(100,100,0.06,1,0.25,1187.5,60,0.03,0.2,F)
FiniteExplicit(100,100,0.06,1,0.25,1187.5,60,0.03,0.2,T)
FiniteImplicit(100,100,1,0.06,0.2,30,30,0.03,0.2,T)
