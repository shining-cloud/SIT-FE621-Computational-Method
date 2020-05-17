% Heston and Black Scholes PDE
clc; clear;

%% Settings and parameters
S0 = 100;
K = 100;
r = 0.05;
q = 0.03;
rho = -0.8;
kappa = 5;
lambda = 0;
sigma = 0.4;
v0 = 0.07;
theta = v0;
PutCall = 'C';
trap = 1;
T = 1;

%% Heston Greeks
[x w] = GenerateGaussLaguerre(32);
Price = HestonPriceGaussLaguerre(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
Delta = HestonGreeks(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Delta');
Gamma = HestonGreeks(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Gamma');
Theta = HestonGreeks(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Theta');
Vega1 = HestonGreeks(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vega1') / (2*sqrt(v0));
Vanna = HestonGreeks(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vanna') / (2*sqrt(v0));
Volga = HestonGreeks(PutCall,S0,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Volga');
Volga = (1/4/sqrt(v0)*Volga - Vega1/2/sqrt(v0))/sqrt(v0);



%% Heston PDE
 Heston = Theta + 0.5*v0*S0^2*Gamma + (r-q)*S0*Delta - r*Price ...
         + rho*sigma*v0*S0*Vanna + 0.5*sigma^2*v0*Volga + kappa*(theta-v0)*Vega1;

%% Black Scholes Greeks
Theta = BSGreeks(PutCall,S0,K,r,q,T,sigma,'Theta');
Gamma = BSGreeks(PutCall,S0,K,r,q,T,sigma,'Gamma');
Delta = BSGreeks(PutCall,S0,K,r,q,T,sigma,'Delta');
Price = BSPrice(PutCall,S0,K,r,q,T,sigma);

%% Black Scholes PDE
BS = Theta + 0.5*sigma^2*S0^2*Gamma + (r-q)*S0*Delta - r*Price;

%% Output the result
fprintf('The Heston PDE equals         %12.2e \n',Heston);
fprintf('The Black-Scholes PDE equals  %12.2e \n',BS);

