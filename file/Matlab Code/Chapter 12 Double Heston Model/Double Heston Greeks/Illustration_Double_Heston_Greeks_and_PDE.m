% Greeks from the Double Heston model.
% Using settings from from P. Gauthier and D. Possamai (2010);

clc; clear;

% Spot price, strike, risk free rate, dividend yield, maturity
S = 61.90;
K = 61.90;
rf = 0.03;
q = 0.0;
T = 1-.25;

% Parameters values
v01 = 0.6^2; 
v02 = 0.7^2;
sigma1 = 0.10;
sigma2 = 0.20;
kappa1 = 0.90;
kappa2 = 1.20;
rho1 = -0.5;
rho2 = -0.5;
theta1 = 0.10;
theta2 = 0.15;

% Stack the parameters into a single vector
param(1) = kappa1;
param(2) = theta1;
param(3) = sigma1;
param(4) =    v01;
param(5) =   rho1;
param(6) = kappa2;
param(7) = theta2;
param(8) = sigma2;
param(9) =    v02;
param(10)=   rho2;

% Define the flavor
PutCall = 'C';

%% Obtain call prices 
% Using the original Christoffersen, Heston, Jacobs characteristic function
% and using the "Little Trap" c.f.
a = 1e-20;
b = 100;
N = 500;
[x w] = GenerateGaussLaguerre(32);
trap = 0;
Price = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);

%% Greeks by finite differences
ds = 0.1;
C1  = DoubleHestonPriceGaussLaguerre(PutCall,S+ds,K,T,rf,q,param,trap,x,w);
C0  = DoubleHestonPriceGaussLaguerre(PutCall,S,   K,T,rf,q,param,trap,x,w);
C1_ = DoubleHestonPriceGaussLaguerre(PutCall,S-ds,K,T,rf,q,param,trap,x,w);
DeltaFD = (C1 - C1_)/2/ds;
GammaFD = (C1 - 2*C0 + C1_)/ds^2;

dr = 0.001;
R1  = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf+dr,q,param,trap,x,w);
R1_ = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf-dr,q,param,trap,x,w);
RhoFD = (R1 - R1_)/2/dr;

dt = 0.001;
T1  = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T+dt,rf,q,param,trap,x,w);
T1_ = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T-dt,rf,q,param,trap,x,w);
ThetaFD = -(T1 - T1_)/2/dt;

dv = 0.001;
param(4) = v01+dv;
V1  = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);
param(4) = v01-dv;
V1_ = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);
param(4) = v01;
V0  = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);
Vega11FD = (V1 - V1_)/2/dv*2*sqrt(v01);
dC1 = (V1 - 2*V0 + V1_)/(dv^2);
Volga1FD = 4*sqrt(v01)*(dC1*sqrt(v01) + Vega11FD/4/v01);
clear V0 V1 V1_

param(9) = v02+dv;
V1  = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);
param(9) = v02-dv;
V1_ = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);
param(9) = v02;
V0  = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w);
Vega12FD = (V1 - V1_)/2/dv*2*sqrt(v02);
dC2 = (V1 - 2*V0 + V1_)/(dv^2);
Volga2FD = 4*sqrt(v02)*(dC2*sqrt(v02) + Vega12FD/4/v02);

param(4) = v01+dv;
Va1  = DoubleHestonPriceGaussLaguerre(PutCall,S+ds,K,T,rf,q,param,trap,x,w);
Va3  = DoubleHestonPriceGaussLaguerre(PutCall,S-ds,K,T,rf,q,param,trap,x,w);
param(4) = v01-dv;
Va2  = DoubleHestonPriceGaussLaguerre(PutCall,S+ds,K,T,rf,q,param,trap,x,w);
Va4  = DoubleHestonPriceGaussLaguerre(PutCall,S-ds,K,T,rf,q,param,trap,x,w);
param(4) = v01;
Vanna1FD = (Va1 - Va2 - Va3 + Va4)/4/dv/ds*2*sqrt(v01);
clear Va1 Va2 Va3 Va4

param(9) = v02+dv;
Va1  = DoubleHestonPriceGaussLaguerre(PutCall,S+ds,K,T,rf,q,param,trap,x,w);
Va3  = DoubleHestonPriceGaussLaguerre(PutCall,S-ds,K,T,rf,q,param,trap,x,w);
param(9) = v02-dv;
Va2  = DoubleHestonPriceGaussLaguerre(PutCall,S+ds,K,T,rf,q,param,trap,x,w);
Va4  = DoubleHestonPriceGaussLaguerre(PutCall,S-ds,K,T,rf,q,param,trap,x,w);
param(9) = v02;
Vanna2FD = (Va1 - Va2 - Va3 + Va4)/4/dv/ds*2*sqrt(v02);


%% Greeks in closed form
Delta = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Delta');
Gamma = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Gamma');
Rho   = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Rho');
Theta = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Theta');
Vega11 = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Vega11');
Vega12 = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Vega12');
Volga1 = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Volga1');
Volga2 = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Volga2');
Vanna1 = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Vanna1');
Vanna2 = DoubleHestonGreeks(S,K,T,rf,q,param,trap,x,w,'Vanna2');


%% Output the results
fprintf('-----------------------------------------------------\n');
fprintf('Greek          Analytic    Finite Diff      Error \n');
fprintf('-----------------------------------------------------\n');
fprintf('Price   %15.4f  \n',Price);
fprintf('Delta   %15.4f %12.4f %15.4e\n',Delta,DeltaFD,abs(Delta-DeltaFD));
fprintf('Gamma   %15.4f %12.4f %15.4e\n',Gamma,GammaFD,abs(Gamma-GammaFD));
fprintf('Rho     %15.4f %12.4f %15.4e\n',Rho,RhoFD,abs(Rho-RhoFD));
fprintf('Theta   %15.4f %12.4f %15.4e\n',Theta,ThetaFD,abs(Theta-ThetaFD));
fprintf('Vega11  %15.4f %12.4f %15.4e\n',Vega11,Vega11FD,abs(Vega11-Vega11FD));
fprintf('Vega12  %15.4f %12.4f %15.4e\n',Vega12,Vega12FD,abs(Vega12-Vega12FD));
fprintf('Vanna1  %15.4f %12.4f %15.4e\n',Vanna1,Vanna1FD,abs(Vanna1-Vanna1FD));
fprintf('Vanna2  %15.4f %12.4f %15.4e\n',Vanna2,Vanna2FD,abs(Vanna2-Vanna2FD));
fprintf('Volga1  %15.4f %12.4f %15.4e\n',Volga1,Volga1FD,abs(Volga1-Volga1FD));
fprintf('Volga2  %15.4f %12.4f %15.4e\n',Volga2,Volga2FD,abs(Volga2-Volga2FD));
fprintf('-----------------------------------------------------\n');


%% Illustration of the double Heston PDE
Vega11p = Vega11 / (2*sqrt(v01));
Vega12p = Vega12 / (2*sqrt(v02));
Vanna1p = Vanna1 / (2*sqrt(v01));
Vanna2p = Vanna2 / (2*sqrt(v02));
Volga1p = (1/4/sqrt(v01)*Volga1 - Vega11p/2/sqrt(v01))/sqrt(v01);
Volga2p = (1/4/sqrt(v02)*Volga2 - Vega12p/2/sqrt(v02))/sqrt(v02);

%% The double Heston PDE
A = (rf-q)*S*Delta ...
    + kappa1*(theta1-v01)*Vega11p ...
    + kappa2*(theta2-v02)*Vega12p ...
    + 0.5*(v01+v02)*S^2*Gamma ...
    + rho1*sigma1*v01*S*Vanna1p ...
    + rho2*sigma2*v02*S*Vanna2p ...
    + 0.5*sigma1^2*v01*Volga1p ...
    + 0.5*sigma2^2*v02*Volga2p;
PDE = Theta + A - rf*Price;

%% Output the results
fprintf('Double Heston PDE evaluated at the Greeks\n');
fprintf('-----------------------------------------\n');
fprintf('Call price     %10.4f \n',Price);
fprintf('PDE value      %12.2e \n',PDE);
fprintf('-----------------------------------------\n');

