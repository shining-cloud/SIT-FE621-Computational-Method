% Greeks by exact formula
clc; clear;

%% Settings and parameters
S = 100;         % Spot price.
K = 100;         % Strike
T = 0.25;        % Time to maturity.
r = 0.05;        % Risk free rate.
q = 0.0;         % Dividend
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.1;     % Heston parameter: volatility of vol
v0    = 0.05;    % Heston parameter: initial variance.
rho   = -0.9;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
PutCall = 'C';

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

% Volga
dC2 = (V1 - 2*Price + V2)/(dv^2);
VolgaFD = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD/4/v0);

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
VannaFD = (C1 - C2 - C3 + C4)/2/dv/ds*sqrt(v0);


%% Display the results
fprintf('Greek    ClosedForm    FinDiff \n')
fprintf('---------------------------- \n')
fprintf('Price   %10.4f \n',Price)
fprintf('Delta   %10.4f %10.4f \n',Delta,DeltaFD)
fprintf('Gamma   %10.4f %10.4f \n',Gamma,GammaFD)
fprintf('Theta   %10.4f %10.4f \n',Theta,ThetaFD)
fprintf('Rho     %10.4f %10.4f \n',Rho  ,RhoFD)
fprintf('Vega1   %10.4f %10.4f \n',Vega1,Vega1FD)
fprintf('Vanna   %10.4f %10.4f \n',Vanna,VannaFD)
fprintf('Volga   %10.4f %10.4f \n',Volga,VolgaFD)
fprintf('---------------------------- \n')
