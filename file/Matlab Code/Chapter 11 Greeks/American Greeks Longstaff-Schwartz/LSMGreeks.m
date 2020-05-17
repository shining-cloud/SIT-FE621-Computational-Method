function [Euro Amer] = LSMGreeks(S,K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,Greek)

kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = 0;

ds = 0.01*S;
dv = 0.01*v0;
dt = 0.01*T;
dr = 0.01*r;

if strcmp(Greek,'price')
    [Spaths  V] = MMSimGreeks(params,S,T,r,q,NT,NS,Zv,Zs);
    [EuroPrice AmerPrice]  = LSM(Spaths',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    Euro = EuroPrice;
    Amer = AmerPrice;
end

if strcmp(Greek,'delta') || strcmp(Greek,'gamma')
    [Spathsp V] = MMSimGreeks(params,S+ds,T,r,q,NT,NS,Zv,Zs);
    [Spathsm V] = MMSimGreeks(params,S-ds,T,r,q,NT,NS,Zv,Zs);
    [EuroPricep AmerPricep] = LSM(Spathsp',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricem AmerPricem] = LSM(Spathsm',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    if strcmp(Greek,'delta')
        Euro = (EuroPricep - EuroPricem)/2/ds;
        Amer = (AmerPricep - AmerPricem)/2/ds;
    else
        [Spaths V] = MMSimGreeks(params,S,T,r,q,NT,NS,Zv,Zs);
        [EuroPrice AmerPrice] = LSM(Spaths',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
        Amer = (AmerPricep - 2*AmerPrice + AmerPricem)/ds^2;
        Euro = (EuroPricep - 2*EuroPrice + EuroPricem)/ds^2;
    end
end

if strcmp(Greek,'vega1')
    paramsp = [kappa theta sigma v0+dv rho lambda];
    paramsm = [kappa theta sigma v0-dv rho lambda];
    [Spathsp V] = MMSimGreeks(paramsp,S,T,r,q,NT,NS,Zv,Zs);
    [Spathsm V] = MMSimGreeks(paramsm,S,T,r,q,NT,NS,Zv,Zs);
    [EuroPricep AmerPricep] = LSM(Spathsp',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricem AmerPricem] = LSM(Spathsm',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    Euro = (EuroPricep - EuroPricem)/2/dv*2*sqrt(v0);
    Amer = (AmerPricep - AmerPricem)/2/dv*2*sqrt(v0);
end

if strcmp(Greek,'vanna')
    paramsp = [kappa theta sigma v0+dv rho lambda];
    paramsm = [kappa theta sigma v0-dv rho lambda];
    [Spathspp V] = MMSimGreeks(paramsp,S+ds,T,r,q,NT,NS,Zv,Zs);
    [Spathspm V] = MMSimGreeks(paramsm,S+ds,T,r,q,NT,NS,Zv,Zs);
    [Spathsmp V] = MMSimGreeks(paramsp,S-ds,T,r,q,NT,NS,Zv,Zs);
    [Spathsmm V] = MMSimGreeks(paramsm,S-ds,T,r,q,NT,NS,Zv,Zs);
    [EuroPricepp AmerPricepp] = LSM(Spathspp',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricepm AmerPricepm] = LSM(Spathspm',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricemp AmerPricemp] = LSM(Spathsmp',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricemm AmerPricemm] = LSM(Spathsmm',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    Euro = (EuroPricepp - EuroPricepm - EuroPricemp + EuroPricemm)/4/dv/ds*2*sqrt(v0);
    Amer = (AmerPricepp - AmerPricepm - AmerPricemp + AmerPricemm)/4/dv/ds*2*sqrt(v0);
end

if strcmp(Greek,'theta')
    [Spathsp V] = MMSimGreeks(params,S,T+dt,r,q,NT,NS,Zv,Zs);
    [Spathsm V] = MMSimGreeks(params,S,T-dt,r,q,NT,NS,Zv,Zs);
    [EuroPricep AmerPricep] = LSM(Spathsp',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricem AmerPricem] = LSM(Spathsm',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    Euro = -(EuroPricep - EuroPricem)/2/dt;
    Amer = -(AmerPricep - AmerPricem)/2/dt;
end    

if strcmp(Greek,'rho')
    [Spathsp V] = MMSimGreeks(params,S,T,r+dr,q,NT,NS,Zv,Zs);
    [Spathsm V] = MMSimGreeks(params,S,T,r-dr,q,NT,NS,Zv,Zs);
    [EuroPricep AmerPricep] = LSM(Spathsp',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    [EuroPricem AmerPricem] = LSM(Spathsm',K,r,q,T,NT,NS,PutCall,XmatrixHandle);
    Euro = (EuroPricep - EuroPricem)/2/dr;
    Amer = (AmerPricep - AmerPricem)/2/dr;
end    
