% Greeks by exact and finite difference formulas
clc; clear;

%% Settings and parameters
S = 100;
K = 100;
r = 0.05;
q = 0.03;
T = 0.5 ;
rho = -0.8;
kappa = 0;
lambda = 0;
sigma = 1e-5;
v0 = 0.07;
theta = v0;
PutCall = 'C';
trap = 1;


%% Closed form Greeks
[x w] = GenerateGaussLaguerre(32);
Price = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
Delta = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Delta');
Gamma = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Gamma');
Rho   = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Rho');
Vega1 = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vega1');
Theta = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Theta');
Vanna = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vanna');
Volga = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Volga');

%% Finite difference Greeks
% Delta
dS = 1 ;
D1 = HestonPriceGaussLaguerre(PutCall,S+dS,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
D2 = HestonPriceGaussLaguerre(PutCall,S-dS,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
DeltaFD = (D1-D2)/2/dS;

% Gamma
GammaFD = (D1 - 2*Price + D2)/dS^2;

% Rho
dr = 1e-5;
R1 = HestonPriceGaussLaguerre(PutCall,S,K,T,r+dr,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
R2 = HestonPriceGaussLaguerre(PutCall,S,K,T,r-dr,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
RhoFD = (R1-R2)/2/dr;

% Vega #1  
dv = 1e-5;
V1 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0+dv,rho,trap,x,w);
V2 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0-dv,rho,trap,x,w);
Vega1FD = (V1-V2)/2/dv*2*sqrt(v0);

% Theta
dt = 1e-2;
T1 = HestonPriceGaussLaguerre(PutCall,S,K,T-dt,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
T2 = HestonPriceGaussLaguerre(PutCall,S,K,T+dt,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
ThetaFD = (T1-T2)/2/dt;
	
% Vanna
dv = 1e-5;
ds = 1e-1;
C1 = HestonPriceGaussLaguerre(PutCall,S+ds,K,T,r,q,kappa,theta,sigma,lambda,v0+dv,rho,trap,x,w);
C2 = HestonPriceGaussLaguerre(PutCall,S+ds,K,T,r,q,kappa,theta,sigma,lambda,v0-dv,rho,trap,x,w);
C3 = HestonPriceGaussLaguerre(PutCall,S-ds,K,T,r,q,kappa,theta,sigma,lambda,v0+dv,rho,trap,x,w);
C4 = HestonPriceGaussLaguerre(PutCall,S-ds,K,T,r,q,kappa,theta,sigma,lambda,v0-dv,rho,trap,x,w);
VannaFD = (C1 - C2 - C3 + C4)/4/dv/ds*2*sqrt(v0);
	
% Volga
dC2 = (V1 - 2*Price + V2)/(dv^2);
VolgaFD = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD/4/v0);

	
%% BlackScholes price and Greeks
vol  = sqrt(v0);
BS   = BSPrice(PutCall,S,K,r,q,T,vol);
BSD  = BSGreeks(PutCall,S,K,r,q,T,vol,'Delta');
BSG  = BSGreeks(PutCall,S,K,r,q,T,vol,'Gamma');
BSV  = BSGreeks(PutCall,S,K,r,q,T,vol,'Vega');
BSR  = BSGreeks(PutCall,S,K,r,q,T,vol,'Rho');
BST  = BSGreeks(PutCall,S,K,r,q,T,vol,'Theta');
BSVa = BSGreeks(PutCall,S,K,r,q,T,vol,'Vanna');
BSVo = BSGreeks(PutCall,S,K,r,q,T,vol,'Volga');

   
%% Display the results
fprintf('Exact Price         %10.5f\n', Price)
fprintf('Black Scholes Price %10.5f\n', BS)
fprintf('------------------------------------------------------\n')
fprintf('Greek          Exact        FD       Black-Scholes\n')
fprintf('------------------------------------------------------\n')
fprintf('Delta      %10.4f  %10.4f  %10.4f \n', Delta,DeltaFD,BSD)
fprintf('Gamma      %10.4f  %10.4f  %10.4f \n', Gamma,GammaFD,BSG)
fprintf('Rho        %10.4f  %10.4f  %10.4f \n', Rho,RhoFD,BSR)
fprintf('Theta      %10.4f  %10.4f  %10.4f \n', Theta,ThetaFD,BST)
fprintf('Vega #1    %10.4f  %10.4f  %10.4f \n', Vega1,Vega1FD,BSV)
fprintf('Vanna      %10.4f  %10.4f  %10.4f \n', Vanna,VannaFD,BSVa)
fprintf('Volga      %10.4f  %10.4f  %10.4f \n', Volga,VolgaFD,BSVo)
fprintf('\n')

%% Calculate and display the other Heston sensitivities

% Change parameter settings
kappa = 5;
sigma = 0.35;

% Closed form
Vega2 = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vega2');
Corr  = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Corr');
Sigma = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Sigma');
Kappa = HestonGreeks(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Kappa');

% Finite differences
% Vega #2
dh = 1e-4;
H1 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta+dh,sigma,lambda,v0,rho,trap,x,w); 
H2 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta-dh,sigma,lambda,v0,rho,trap,x,w); 
Vega2FD = (H1-H2)/2/dh*2*sqrt(theta);

% Correlation rho
dp = 1e-2;	
C1 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho+dp,trap,x,w);
C2 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho-dp,trap,x,w);
CorrFD = (C1-C2)/2/dp;

% Volatility of variance sigma
ds = 1e-2;	
S1 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma+ds,lambda,v0,rho,trap,x,w);
S2 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma-ds,lambda,v0,rho,trap,x,w);
SigmaFD = (S1-S2)/2/ds;

% Mean reversion speed kappa
dk = 1e-5;
K1 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa+dk,theta,sigma,lambda,v0,rho,trap,x,w);
K2 = HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa-dk,theta,sigma,lambda,v0,rho,trap,x,w);
KappaFD = (K1-K2)/2/dk;

fprintf('Other Sensitivities\n')
fprintf('------------------------------------------------------\n')
fprintf('Greek                 Exact             FD  \n')
fprintf('------------------------------------------------------\n')
fprintf('Vega #2         %15.5d  %15.5d  \n', Vega2,Vega2FD)
fprintf('Correlation     %15.5d  %15.5d  \n', Corr ,CorrFD )
fprintf('Kappa           %15.5d  %15.5d  \n', Kappa,KappaFD)
fprintf('Sigma           %15.5d  %15.5d  \n', Sigma,SigmaFD)
