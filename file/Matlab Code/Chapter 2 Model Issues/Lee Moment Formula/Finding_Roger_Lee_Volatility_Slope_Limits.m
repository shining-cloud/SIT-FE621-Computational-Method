clc; clear;

% Underlying settings
tau = 37/365;
S = 129.14;
r = 0.00;
q = 0.00;
trap = 1;

% Parameter settings from S&P estimation
kappa =  4.5150;
theta =  0.1250;
sigma =  2.1149;
v0    =  0.0200;
rho   = -0.2741;
param = [kappa theta sigma v0 rho];

% Moment explosion setting
HiLimit =  10;
LoLimit = -10;

%% Find the moment bounds and the Lee slope bounds
[bR bL LowerAP UpperAP LowerCF UpperCF] = FindLeeBounds(S,r,q,tau,param,trap,LoLimit,HiLimit);

%% Output the moment bounds results
fprintf('Moment bound method       Lower Bound        Upper Bound\n')
fprintf('----------------------------------------------------------\n')
fprintf('Andersen & Piterbarg    %10.4f          %10.4f\n', LowerAP,UpperAP)
fprintf('Using the Heston C.F.   %10.4f          %10.4f\n', LowerCF,UpperCF)
fprintf('----------------------------------------------------------\n')
fprintf('\n');

%% Output the slope limits
fprintf('Roger Lee LHS slope limit  %10.4f\n', bL)
fprintf('Roger Lee RHS slope limit  %10.4f\n', bR)
fprintf('\n');
