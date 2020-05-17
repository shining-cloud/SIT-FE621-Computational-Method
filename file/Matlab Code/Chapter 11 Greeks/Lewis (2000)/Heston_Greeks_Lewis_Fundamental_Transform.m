% Heston (1993) prices using the Fundamental transform described in Alan
% Lewis' book "Option Valuation Under Stochastic Volatility: With Mathematica Code"
% By Fabrice Douglas Rouah
% Finds the call price based on C1(S,K,t) and C2(S,K,t).

clc; clear;

%% Spot price, strike price, risk free rate, dividend yield, and maturity
S = 100;
K = 100;
rf = 0.05;
q  = 0.0;
tau = 0.25 ;

% Heston parameters
kappa = 2;
theta = 0.05;
sigma = 0.1;
v0    = 0.05;
rho   = -0.9;
lambda = 0;

% Gauss Laguerre 32-point rule
[x w] = GenerateGaussLaguerre(32);


%% First expression C1(K) for the Lewis Greeks
ki = 1.5;
form = 1;

% Price
PriceFD1 = HestonLewisCallPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);

% Delta and Gamma
ds = 0.1;
C1  = HestonLewisCallPrice(S+ds,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
C0  = HestonLewisCallPrice(S   ,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
C1_ = HestonLewisCallPrice(S-ds,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
DeltaFD1 = (C1 - C1_)/2/ds;
GammaFD1 = (C1 - 2*C0 + C1_)/ds^2;

% Rho
dr = 1e-10;
R1  = HestonLewisCallPrice(S,K,rf+dr,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
R1_ = HestonLewisCallPrice(S,K,rf-dr,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
RhoFD1 = (R1 - R1_)/2/dr;

% Theta
dt = 1e-2;
T1  = HestonLewisCallPrice(S,K,rf,q,v0,tau+dt,ki,theta,kappa,sigma,rho,form,x,w);
T1_ = HestonLewisCallPrice(S,K,rf,q,v0,tau-dt,ki,theta,kappa,sigma,rho,form,x,w);
ThetaFD1 = -(T1 - T1_)/2/dt;

% Vega and Volga
dv = 1e-5;
Ve1  = HestonLewisCallPrice(S,K,rf,q,v0+dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Ve0  = HestonLewisCallPrice(S,K,rf,q,v0   ,tau,ki,theta,kappa,sigma,rho,form,x,w);
Ve1_ = HestonLewisCallPrice(S,K,rf,q,v0-dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Vega1FD1 = (Ve1 - Ve1_)/2/dv*2*sqrt(v0);
dC2 = (Ve1 - 2*Ve0 + Ve1_)/(dv^2);
VolgaFD1 = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD1/4/v0);

% Vanna
Va1  = HestonLewisCallPrice(S+ds,K,rf,q,v0+dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Va2  = HestonLewisCallPrice(S+ds,K,rf,q,v0-dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Va3  = HestonLewisCallPrice(S-ds,K,rf,q,v0+dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Va4  = HestonLewisCallPrice(S-ds,K,rf,q,v0-dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
VannaFD1 = (Va1 - Va2 - Va3 + Va4)/4/dv/ds*2*sqrt(v0);

% Closed form Greeks from C1(K)
Price1 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Price');
Delta1 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Delta');
Gamma1 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Gamma');
Rho1   = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Rho');
Theta1 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Theta');
Vega11 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Vega1');
Vanna1 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Vanna');
Volga1 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Volga');

%% Second expression C1(K) for the Lewis Greeks
ki = 0.5;
form = 2;

% Price
PriceFD2 = HestonLewisCallPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);

% Delta and Gamma
C1  = HestonLewisCallPrice(S+ds,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
C0  = HestonLewisCallPrice(S   ,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
C1_ = HestonLewisCallPrice(S-ds,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
DeltaFD2 = (C1 - C1_)/2/ds;
GammaFD2 = (C1 - 2*C0 + C1_)/ds^2;

% Rho
R1  = HestonLewisCallPrice(S,K,rf+dr,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
R0  = HestonLewisCallPrice(S,K,rf   ,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
R1_ = HestonLewisCallPrice(S,K,rf-dr,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w);
RhoFD2 = (R1 - R1_)/2/dr;

% Theta
T1  = HestonLewisCallPrice(S,K,rf,q,v0,tau+dt,ki,theta,kappa,sigma,rho,form,x,w);
T0  = HestonLewisCallPrice(S,K,rf,q,v0,tau   ,ki,theta,kappa,sigma,rho,form,x,w);
T1_ = HestonLewisCallPrice(S,K,rf,q,v0,tau-dt,ki,theta,kappa,sigma,rho,form,x,w);
ThetaFD2 = -(T1 - T1_)/2/dt;

% Vega and Volga
Ve1  = HestonLewisCallPrice(S,K,rf,q,v0+dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Ve0  = HestonLewisCallPrice(S,K,rf,q,v0   ,tau,ki,theta,kappa,sigma,rho,form,x,w);
Ve1_ = HestonLewisCallPrice(S,K,rf,q,v0-dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Vega1FD2 = (Ve1 - Ve1_)/2/dv*2*sqrt(v0);
dC2 = (Ve1 - 2*Ve0 + Ve1_)/(dv^2);
VolgaFD2 = 4*sqrt(v0)*(dC2*sqrt(v0) + Vega1FD2/4/v0);

% Vanna
Va1  = HestonLewisCallPrice(S+ds,K,rf,q,v0+dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Va2  = HestonLewisCallPrice(S+ds,K,rf,q,v0-dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Va3  = HestonLewisCallPrice(S-ds,K,rf,q,v0+dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
Va4  = HestonLewisCallPrice(S-ds,K,rf,q,v0-dv,tau,ki,theta,kappa,sigma,rho,form,x,w);
VannaFD2 = (Va1 - Va2 - Va3 + Va4)/4/dv/ds*2*sqrt(v0);

% Closed form Greeks from C1(K)
Price2 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Price');
Delta2 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Delta');
Gamma2 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Gamma');
Rho2   = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Rho');
Theta2 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Theta');
Vega12 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Vega1');
Vanna2 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Vanna');
Volga2 = HestonLewisGreekPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,'Volga');


%% Output the results
fprintf('Comparison of Lewis (2000) C1(K) and C2(K) Greeks \n');
fprintf(' \n');
fprintf('              ------- C1(K) ---------   ------ C2(K) --------\n');
fprintf('Greek         Analytic    Finite Diff   Analytic  Finite Diff\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Price        %8.4f %12.4f %12.4f %12.4f\n',Price1,PriceFD1,Price2,PriceFD2);
fprintf('Delta        %8.4f %12.4f %12.4f %12.4f\n',Delta1,DeltaFD1,Delta2,DeltaFD2);
fprintf('Gamma        %8.4f %12.4f %12.4f %12.4f\n',Gamma1,GammaFD1,Gamma2,GammaFD2);
fprintf('Theta        %8.4f %12.4f %12.4f %12.4f\n',Theta1,ThetaFD1,Theta2,ThetaFD2);
fprintf('Rho          %8.4f %12.4f %12.4f %12.4f\n',Rho1,RhoFD1,Rho2,RhoFD2);
fprintf('Vega1        %8.4f %12.4f %12.4f %12.4f\n',Vega11,Vega1FD1,Vega12,Vega1FD2);
fprintf('Vanna        %8.4f %12.4f %12.4f %12.4f\n',Vanna1,VannaFD1,Vanna2,VannaFD2);
fprintf('Volga        %8.4f %12.4f %12.4f %12.4f\n',Volga1,VolgaFD1,Volga2,VolgaFD2);
fprintf('--------------------------------------------------------------\n');


