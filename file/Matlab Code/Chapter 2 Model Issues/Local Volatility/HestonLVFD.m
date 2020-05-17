function y = HestonLVFD(S,K,T,rf,kappa,theta,sigma,lambda,v0,rho,Trap,x,w,dt,dK)

% Heston Local Volatility using Dupire's formula.
% Uses finite differences for the required derivatives

% dC/dT by central finite difference
CT_1 = HestonCallGaussLaguerre(S,K,T-dt,rf,kappa,theta,sigma,lambda,v0,rho,Trap,x,w);
CT1  = HestonCallGaussLaguerre(S,K,T+dt,rf,kappa,theta,sigma,lambda,v0,rho,Trap,x,w);
dCdT =  (CT1 - CT_1) / (2*dt);

% dC2/dK2 by central finite differences
CK_1 = HestonCallGaussLaguerre(S,K-dK,T,rf,kappa,theta,sigma,lambda,v0,rho,Trap,x,w);
CK0  = HestonCallGaussLaguerre(S,K   ,T,rf,kappa,theta,sigma,lambda,v0,rho,Trap,x,w);
CK1  = HestonCallGaussLaguerre(S,K+dK,T,rf,kappa,theta,sigma,lambda,v0,rho,Trap,x,w);
dC2dK2 = (CK_1 - 2*CK0 + CK1) / (dK)^2;

% Local variance
LocalVar = 2*dCdT / (K^2*dC2dK2);

% Local volatility
y = sqrt(LocalVar);


