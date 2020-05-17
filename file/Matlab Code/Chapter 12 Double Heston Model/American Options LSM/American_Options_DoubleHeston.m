clc; clear;

% Spot price, risk free rate, dividend yield
S0 = 61.90;
K  = S0;
tau = 1;
rf = 0.03;
q  = 0;
PutCall = 'P';
trap = 1;

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

%% Obtain European call price using closed form 
[x w] = GenerateGaussLaguerre(32);
EuroClosed = DoubleHestonPriceGaussLaguerre(PutCall,S0,K,tau,rf,q,param,x,w,trap);

%% Simulation settings
XmatrixHandle = {@(y)ones(length(y),1), @(y)(1-y),@(y)1./2.*(2-4.*y-y.^2)};
NS = 50000;
NT = 100;

%% Select the simulation scheme and run the simulation
scheme = 'ZhuTV';

tic;
if strcmp(scheme,'Euler') || strcmp(scheme,'Alfonsi')
    [S V1 V2 EuroSim] = DHEulerAlfonsiSim(scheme,param,S0,K,tau,rf,q,NT,NS,PutCall);
elseif strcmp(scheme,'ZhuEuler') || strcmp(scheme,'ZhuTV')
    [S v1 v2 EuroSim] = DHTransVolSim(scheme,param,S0,K,tau,rf,q,NT,NS,PutCall);
elseif strcmp(scheme,'QE')
    [S V1 V2 EuroSim] = DHQuadExpSim(param,S0,K,tau,rf,q,NT,NS,PutCall);
end
SimTime = toc;

% Longstaff-Schwartz Algorithm and control variate price
[EuroLSM AmerLSM] = LSM(S',K,rf,q,tau,NT,NS,PutCall,XmatrixHandle);
AmerCV  = AmerLSM + (EuroClosed - EuroLSM);
AmerCV = EuroClosed + (AmerLSM - EuroLSM);

%% Output the results
fprintf('American Double Heston price using LSM\n');
fprintf('Number of simulations:%5.0f\n',NS);
fprintf('Number of time steps :%5.0f\n',NT);
fprintf('-----------------------------------------------------------------------\n');
fprintf('Method              Euro       Amer     ControlVar  SimTime\n');
fprintf('-----------------------------------------------------------------------\n');
fprintf('Closed          %10.4f\n',EuroClosed);
fprintf('Simulation      %10.4f %10.4f %10.4f %10.4f\n',EuroLSM,AmerLSM,AmerCV,SimTime);
fprintf('-----------------------------------------------------------------------\n');

