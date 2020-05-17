% Heston (1993) prices using the characteristic function described in Alan
% Lewis' paper "A Simple Option Formula for General Jump-Diffusion and other
% Exponential Levy Processes"
% Finds the call price in Equation (3.11) of his paper

clc; clear;

%% Spot price, strike price, risk free rate, dividend yield, and maturity
S = 100;
K = 100;
r = 0.05;
q = 0.0;
T = 0.25 ;

% Heston parameters
kappa = 2;
theta = 0.05;
sigma = 0.1;
v0    = 0.05;
rho   = -0.9;
lambda = 0;

% Absicissas and weights
[x w] = GenerateGaussLaguerre(32);
trap = 1;

%% Finite differences
Price = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Price');

% Delta and Gamma
ds = 0.1;
C1  = LewisGreeks311(S+ds,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
C0  = LewisGreeks311(S   ,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
C1_ = LewisGreeks311(S-ds,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
DeltaFD = (C1 - C1_)/2/ds;
GammaFD = (C1 - 2*C0 + C1_)/ds^2;

% Rho
dr = 1e-4;
R1  = LewisGreeks311(S,K,r+dr,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
R1_ = LewisGreeks311(S,K,r-dr,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
RhoFD = (R1 - R1_)/2/dr;

% Theta
dt = 1e-5;
T1  = LewisGreeks311(S,K,r,q,T+dt,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
T1_ = LewisGreeks311(S,K,r,q,T-dt,theta,kappa,sigma,rho,v0,trap,x,w,'Price');
ThetaFD = -(T1 - T1_)/2/dt;

% Vega and Volga
dv = 1e-4;
V1  = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0+dv,trap,x,w,'Price');
V0  = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,   trap,x,w,'Price');
V1_ = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0-dv,trap,x,w,'Price');
Vega1FD = (V1 - V1_)/2/dv * 2*sqrt(v0);
dC2 = (V1 - 2*V0 + V1_)/(dv^2);
VolgaFD = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD/4/v0);

Va1 = LewisGreeks311(S+ds,K,r,q,T,theta,kappa,sigma,rho,v0+dv,trap,x,w,'Price');
Va2 = LewisGreeks311(S+ds,K,r,q,T,theta,kappa,sigma,rho,v0-dv,trap,x,w,'Price');
Va3 = LewisGreeks311(S-ds,K,r,q,T,theta,kappa,sigma,rho,v0+dv,trap,x,w,'Price');
Va4 = LewisGreeks311(S-ds,K,r,q,T,theta,kappa,sigma,rho,v0-dv,trap,x,w,'Price');
VannaFD = (Va1 - Va2 - Va3 + Va4)/4/dv/ds*2*sqrt(v0);

% Closed form analytic Greeks
Delta = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Delta');
Gamma = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Gamma');
Rho   = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Rho');
Theta = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Theta');
Vega1 = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Vega1');
Volga = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Volga');
Vanna = LewisGreeks311(S,K,r,q,T,theta,kappa,sigma,rho,v0,trap,x,w,'Vanna');


%% Output the results;
fprintf('Greeks under the Lewis (2001) formulation \n');
fprintf('----------------------------------------\n');
fprintf('Greek            Analytic   Finite Diff \n');
fprintf('----------------------------------------\n');
fprintf('Price   %15.4f  \n',Price);
fprintf('Delta   %15.4f %12.4f \n',Delta,DeltaFD);
fprintf('Gamma   %15.4f %12.4f \n',Gamma,GammaFD);
fprintf('Theta   %15.4f %12.4f \n',Theta,ThetaFD);
fprintf('Rho     %15.4f %12.4f \n',Rho,RhoFD);
fprintf('Vega1   %15.4f %12.4f \n',Vega1,Vega1FD);
fprintf('Vanna   %15.4f %12.4f \n',Vanna,VannaFD);
fprintf('Volga   %15.4f %12.4f \n',Volga,VolgaFD);
fprintf('----------------------------------------\n');





