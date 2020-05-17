# Assignment 4
#Question 1 ----
setwd("~/R/Courses/FE 621/Assignments/Assignment 4")
Sample_Data = read.csv("sample_data.csv", header = T)
Sample_Data
D1 = Sample_Data$stock1
D2 = Sample_Data$stock2
D3 = Sample_Data$stock3
D4 = Sample_Data$stock4
D5 = Sample_Data$stock5
#install.packages("Sim.DiffProc")
library(Sim.DiffProc)

ModelCheck <- function(S){
    X = ts(S)
    fx <- expression(theta[1]*x)
    gx<- expression(theta[2]*x^theta[3])
    mod <- fitsde(data = X, drift = fx, diffusion = gx, start = list(theta1=1, theta2=1,theta3=1),pmle="euler")
    
    fx1 <- expression( theta[1]+theta[2]*x )
    gx1 <- expression( theta[3]*x^theta[4] )
    mod1 <- fitsde(data=X,drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="euler")
    
    
    fx2 <- expression( theta[1]+theta[2]*x)
    gx2 <- expression( theta[3]*sqrt(x) )
    mod2 <- fitsde(data=X,drift=fx2,diffusion=gx2,start = list(theta1=1, theta2=1,theta3=1),pmle="euler")
    
    fx3 <- expression( theta[1] )
    gx3 <- expression( theta[2]*x^theta[3] )
    mod3 <- fitsde(data=X,drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1),pmle="euler")
    
    fx4 = expression(theta[1]*x)
    gx4 = expression(theta[2] + theta[3]*x^theta[4])
    mod4 <- fitsde(data=X,drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="euler")
    
    aic <- c(AIC(mod),AIC(mod1),AIC(mod2),AIC(mod3),AIC(mod4))
    bic <- c(BIC(mod),BIC(mod1),BIC(mod2),BIC(mod3),BIC(mod4))
    Testaic<- data.frame(aic,row.names = c("mod1","mod2","mod3","mod4", "mod5"))
    Testbic <- data.frame(bic,row.names = c("mod1","mod2","mod3","mod4", "mod5"))
    TrueModAIC = rownames(Testaic)[which.min(Testaic[,1])]
    TrueModBIC = rownames(Testbic)[which.min(Testbic[,1])]
 
     list(AICModel = Testaic,BICModel = Testbic, BestfitAIC = TrueModAIC, BestfitBIC= TrueModBIC)
    
    
}
ModelCheck(D1) # Best Model: Model 2
ModelCheck(D2) # Best Model: Model 4
ModelCheck(D3) # Best Model: Model 1,3,4
ModelCheck(D4) # Best Model: Model 2
ModelCheck(D5) # Best Model: Model 4

# part 2 
# Data 1 Best Model is model 2----
# Euler
fx1 <- expression( theta[1]+theta[2]*x )
gx1 <- expression( theta[3]*x^theta[4] )
Eumod1 <- fitsde(data=ts(D1),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="euler")
coef(Eumod1)
summary(Eumod1)
vcov(Eumod1)
confint(Eumod1, level=0.95)
# Ozaki
Ozmod1 <- fitsde(data=ts(D1),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="ozaki")
summary(Ozmod1)
#Shoji-Ozaki
Shmod1 <-fitsde(data=ts(D1),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="shoji")
coef(Shmod1)
Summary(Shmod1) 
# Kesler
Kesmod1 <-fitsde(data=ts(D1),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="kessler")
coef(Kesmod1)
Summary(Kesmod1) 
# Data 2: Best Model is model 4 ----
fx3 <- expression( theta[1] )
gx3 <- expression( theta[2]*x^theta[3] )
Eumod2 <- fitsde(data=ts(D2),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1, theta4=1),pmle="euler")
Summary(Eumod2)
coef(Eumod2)
# Ozaki
Ozmod2 <- fitsde(data=ts(D2),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="ozaki")
summary(Ozmod2)
coef(Ozmod2)
#Shoji-Ozaki
Shmod2 <-fitsde(data=ts(D2),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="shoji")
coef(Shmod2)
Summary(Shmod2) 
# Kesler
Kesmod2 <-fitsde(data=ts(D1),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="kessler")
coef(Kesmod2)
Summary(Kesmod2) 

# Data 4: Best Model is model 2----
fx1 <- expression( theta[1]+theta[2]*x )
gx1 <- expression( theta[3]*x^theta[4] )
Eumod4 <- fitsde(data=ts(D4),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="euler")
coef(Eumod4)
summary(Eumod4)
# Ozaki
Ozmod4 <- fitsde(data=ts(D4),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="ozaki")
coef(Ozmod4)
summary(Ozmod)
#Shoji-Ozaki
Shmod4 <-fitsde(data=ts(D4),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="shoji")
coef(Shmod4)
Summary(Shmod4) 
# Kesler
Kesmod4 <-fitsde(data=ts(D4),drift=fx1,diffusion=gx1,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="kessler")
coef(Kesmod4)
Summary(Kesmod4) 

# Data 5: Best Model is Model 4 ----
fx3 <- expression( theta[1] )
gx3 <- expression( theta[2]*x^theta[3] )
Eumod5 <- fitsde(data=ts(D5),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="euler")
Summary(Eumod5)
coef(Eumod5)
# Ozaki
Ozmod5 <- fitsde(data=ts(D5),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="ozaki")
summary(Ozmod5)
coef(Ozmod5)
#Shoji-Ozaki
Shmod5 <-fitsde(data=ts(D5),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="shoji")
coef(Shmod5)
Summary(Shmod5) 
# Kesler
Kesmod5 <-fitsde(data=ts(D5),drift=fx3,diffusion=gx3,start = list(theta1=1, theta2=1,theta3=1,theta4=1),pmle="kessler")
coef(Kesmod5)
Summary(Kesmod5) 

# Generating data values for model 2
theta1 = coef(Eumod1)

fx1 <- expression( theta1[1]+theta1[2]*x )
gx1 <- expression( theta1[3]*x^theta1[4] )
sim <- snssde1d(drift=fx1,diffusion=gx1,x0=100,M=1,N=100000,Dt=1/365) 
data1=sim$X
error1 = mean((data1-D1)^2)
error1


# Question 2 ----
# Assume values for rho, alpha, v execute equation 4
#
SigmaB <- function(f,al,beta,rho,v,tau){ # Equation 4
 ans =   (al*(1+((1-beta)^2/24)*(al^2)/(f^(2-2*beta))+0.25*(rho*beta*v*al)
         /(f^(1-beta))+(2-3*rho^2)/(24)*v^2)*tau)/(f^(1-beta))
 return(ans)
}

FindAlpha <- function(f,beta,ATM,tau,rho,v){
  
    Model<- function(x){
  c( p1= ((1-beta^2)*tau)/(24*f^(2-2*beta))*x[1]^3 + (rho*beta*v*tau)/(4*f^(1-beta)) * x[1]^2 + (1+(2-3*rho^2)/(24)*v^2*tau) * x[1] - ATM*f^(1-beta))
  }
    ss <- multiroot(f = Model, start = c(1))
   
    return(ss$root)
  
}

FindAlpha(1.47,0.5,0.4291,1,-1,0)
SigmaB(f,FindAlpha(1.47,0.4291,1,-1,0),beta,rho,v,tau)
f = 1.37
beta = 0.5
rho = -1
v = 0
ATM = 0.2879
tau = 1
install.packages("rJava")
library(rJava)
table = read.csv("2017_2_15_mid.csv",header = T)
table

table$X1Yr
j = 1
ATM = c()
f = c()
for(i in 1:(length(table$X1Yr)/2)){
    ATM[i] = table$X1Yr[j]
    f[i] = table$X1Yr[j+1]
    j = j + 2
    }
ATM
f
FindAlpha(f[1],ATM[1],1)
for(i in 1:length(ATM)){}
ss <- multiroot(f = FindAlpha(f,ATM,tau,x), start = c(1))
ss
SigmaB(f,ss$root,beta,rho,v,tau)


f = 2.68
beta = 0.5
rho = 1
v = 1.1
SigmaB(f,0.9857506,beta,rho,v,tau)
tau = 1
SigmaB(f,1.100208,beta,rho,v,tau)

ss <- multiroot(f = FindAlpha, start = c(1))
ss
ss <- multiroot(f = findal, start = c(1))


# Question 3 ----
# Assume alpha, N, eta
alp = 1.75
i = sqrt(as.complex(-1))
N = 512
sigma = 0.15
eta = 0.03
lambda = (2*pi)/(N*eta)
#eta = (2*pi)/(N*lambda)
a = N*eta # upper limit
b = 0.5*N*lambda

CallPrice <- function(S,k,r,tau,sigma,N,i,eta,lambda,alp){
    ans= c()
    k = c()
    
    for(u in 1:N){
        Sum = 0
        for(j in 1:N){
            v = V(j,eta)
             Sum = Sum + exp((-i*2*pi/N)*(j-1)*(u-1))* exp(i* b * V(j,eta))*Si(v,i,S,r,sigma,tau,alp,j,eta) * (eta/3)*(3+(-1)^j-del(j-1))
           
        }
        
    ans[u] = (exp(-alp*Ku(b,u,lambda))/pi)*Re(Sum)
    k[u] = Ku(b,u,lambda)
    }
    list(Price = ans,Strike = k)
}

Ku <- function(b,u,lambda){
   ans = -b + lambda*(u-1)  
   return(ans)
}
phi <- function(v,i,S,r,sigma,tau,j,eta){
    ans = exp(i*(log(S)+((r-sigma^2)/2)*tau)*V(j,eta) - sigma^2*tau*V(j,eta)^2/2)
    return(ans)
}
Si <- function(v,i,S,r,sigma,tau,alp,j,eta){
    v = V(j,eta)
    ans = exp(-r*tau)*phi(v-(alp+1)*i,i,S,r,sigma,tau,j,eta)/(alp^2 + alp - v^2 + i*(2*alp +1)*v)
}

V<- function(j,eta){
    ans = eta*(j-1)
    return(ans)
}

del <- function(N){
    if(N == 0) {return(1)}
    else {return(0)}
}
r = 0.05
tau = 1
S = 100

CallPrice(S,k,r,tau,sigma,N,i,eta,lambda,alp)

# test ----
model2 <- function(x) {
    c(F1 = -3*x[1]^2 - x[1] +6)
}
# first solution

one = multiroot(model2, c(1), useFortran = FALSE)

two = multiroot(model2, c(100), useFortran = FALSE)

