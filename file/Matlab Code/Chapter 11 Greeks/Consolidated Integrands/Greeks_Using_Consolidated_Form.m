clc; clear;

%% Option features
S = 100;
K = 103;
r = 0.05;
q = 0.03;
tau = 0.5;

% Heston parameters
rho = -0.8;
kappa  = 5;
lambda = 0;
sigma  = 0.5;
v0 = 0.07;
theta = 0.05;
trap = 1;

                 
%% Obtain the Greeks
[x w] = GenerateGaussLaguerre(32);

ds = 0.1;
C1  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S+ds,r,q,v0,trap,x,w,'Price');
C0  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S   ,r,q,v0,trap,x,w,'Price');
C1_ = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S-ds,r,q,v0,trap,x,w,'Price');
DeltaFD = (C1 - C1_)/2/ds;
GammaFD = (C1 - 2*C0 + C1_)/ds^2;

dt = 0.01;
T1  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau+dt,K,S,r,q,v0,trap,x,w,'Price');
T1_ = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau-dt,K,S,r,q,v0,trap,x,w,'Price');
ThetaFD = -(T1 - T1_)/2/dt;

dr = 0.01;
R1  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r+dr,q,v0,trap,x,w,'Price');
R1_ = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r-dr,q,v0,trap,x,w,'Price');
RhoFD = (R1 - R1_)/2/dr;

dv = 0.001;
V1  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0+dv,trap,x,w,'Price');
V0  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0   ,trap,x,w,'Price');
V1_ = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0-dv,trap,x,w,'Price');
Vega1FD = (V1 - V1_)/2/dv*2*sqrt(v0);
dC2 = (V1 - 2*V0 + V1_)/(dv^2);
VolgaFD = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD/4/v0);

Va1  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S+ds,r,q,v0+dv,trap,x,w,'Price');
Va2  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S+ds,r,q,v0-dv,trap,x,w,'Price');
Va3  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S-ds,r,q,v0+dv,trap,x,w,'Price');
Va4  = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S-ds,r,q,v0-dv,trap,x,w,'Price');
VannaFD = (Va1 - Va2 - Va3 + Va4)/4/dv/ds*2*sqrt(v0);

Price = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Price');
Delta = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Delta');
Gamma = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Gamma');
Theta = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Theta');
Rho   = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Rho');
Vega1 = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Vega1');
Vanna = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Vanna');
Volga = HestonGreeksConsolidated(kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap,x,w,'Volga');


%% Output the results
fprintf('Greek          Analytic     Finite Diff \n');
fprintf('----------------------------------------\n');
fprintf('Price   %15.4f  \n',Price);
fprintf('Delta   %15.4f %12.4f \n',Delta,DeltaFD);
fprintf('Gamma   %15.4f %12.4f \n',Gamma,GammaFD);
fprintf('Rho     %15.4f %12.4f \n',Rho,RhoFD);
fprintf('Theta   %15.4f %12.4f \n',Theta,ThetaFD);
fprintf('Vega1   %15.4f %12.4f \n',Vega1,Vega1FD);
fprintf('Vanna   %15.4f %12.4f \n',Vanna,VannaFD);
fprintf('Volga   %15.4f %12.4f \n',Volga,VolgaFD);
fprintf('----------------------------------------\n');

