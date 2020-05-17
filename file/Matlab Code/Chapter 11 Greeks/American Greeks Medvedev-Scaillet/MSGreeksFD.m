function y = MSGreeksFD(params,S,K,r,q,T,method,A,B,N,yinf,hi,NumTerms,Greek)

kappaV = params(1);
thetaV = params(2);
sigmaV = params(3);
v0     = params(4);
rho    = params(5);
lambda = 0;
trap = 1;

% Define the finite difference increments
ds = 0.005 * S;
dt = 0.005 * T;
dv = 0.005 * v0;
dr = 0.005 * r;

if strcmp(Greek,'price')
    [x1 x2 AmerPut x3 x4 x5] = MSPrice(S,K,T,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    y = AmerPut;
end

if (strcmp(Greek,'delta') || strcmp(Greek,'gamma'))
    [x1 x2 AmerPutp x3 x4 x5] = MSPrice(S+ds,K,T,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutm x3 x4 x5] = MSPrice(S-ds,K,T,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    if strcmp(Greek,'gamma')
        [x1 x2 AmerPut x3 x4 x5] = MSPrice(S,K,T,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
        y = (AmerPutp - 2*AmerPut + AmerPutm)/ds^2;
    else
        y = (AmerPutp - AmerPutm)/2/ds;
    end
end

if strcmp(Greek,'theta')
    [x1 x2 AmerPutp x3 x4 x5] = MSPrice(S,K,T+dt,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutm x3 x4 x5] = MSPrice(S,K,T-dt,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    y = -(AmerPutp - AmerPutm)/2/dt;
end

if strcmp(Greek,'vega1') | strcmp(Greek,'volga')
    paramsp = [kappaV thetaV sigmaV v0+dv rho];
    paramsm = [kappaV thetaV sigmaV v0-dv rho];
    [x1 x2 AmerPutp x3 x4 x5] = MSPrice(S,K,T,r,q,paramsp,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutm x3 x4 x5] = MSPrice(S,K,T,r,q,paramsm,trap,method,A,B,N,NumTerms,yinf,hi);
    Vega1 = (AmerPutp - AmerPutm)/2/dv*2*sqrt(v0);
    if strcmp(Greek,'volga')
        [x1 x2 AmerPut x3 x4 x5] = MSPrice(S,K,T,r,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
        dC2 = (AmerPutp - 2*AmerPut + AmerPutm)/(dv^2);
        y = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1/4/v0);
    else
        y = Vega1;
    end
end

if strcmp(Greek,'vanna')
    paramsp = [kappaV thetaV sigmaV v0+dv rho];
    paramsm = [kappaV thetaV sigmaV v0-dv rho];
    [x1 x2 AmerPutpp x3 x4 x5] = MSPrice(S+ds,K,T,r,q,paramsp,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutpm x3 x4 x5] = MSPrice(S+ds,K,T,r,q,paramsm,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutmp x3 x4 x5] = MSPrice(S-ds,K,T,r,q,paramsp,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutmm x3 x4 x5] = MSPrice(S-ds,K,T,r,q,paramsm,trap,method,A,B,N,NumTerms,yinf,hi);
    y = (AmerPutpp - AmerPutpm - AmerPutmp + AmerPutmm)/4/dv/ds*2*sqrt(v0);
end

if strcmp(Greek,'rho')
    [x1 x2 AmerPutp x3 x4 x5] = MSPrice(S,K,T,r+dr,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    [x1 x2 AmerPutm x3 x4 x5] = MSPrice(S,K,T,r-dr,q,params,trap,method,A,B,N,NumTerms,yinf,hi);
    y = -(AmerPutp - AmerPutm)/2/dr;
end


