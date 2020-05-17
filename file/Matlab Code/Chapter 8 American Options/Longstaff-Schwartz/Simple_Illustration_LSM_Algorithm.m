% Evaluation of American put options fromt the Heston model.
% Benchmark against values in Clarke and Parrot (1999)

clc; clear;

% Settings from Clarke and Parrot (1999)
S = 10;
K = 10;
T = 1/4;
r = -0.1;
q = 0.0;
PutCall = 'P';
params = [10 .05 .4 .06 -0.6 0];
rho = params(5);

% Number of time steps and number of stock price paths
NT = 7;
NS = 6;

% Design matrix for the LSM algorithm
XmatrixHandle = {@(y)ones(length(y),1), @(y)y};

%% The LSM algorithm
% Zv = randn(NT,NS);
% Zs = rho.*Zv + sqrt(1-rho.^2).*randn(NT,NS);
% [Spaths V] = MMSim(params,S,T,r,q,NT,NS,Zv,Zs);

Spaths = [...
   10.0000   10.0000   10.0000   10.0000   10.0000   10.0000;
    9.7649   10.1610    9.6139    9.5743   10.1996   10.1734;
    9.5249   10.2491    9.3050    9.8990   10.5908    9.2473;
    9.4070   10.8710    9.2125    9.3057   10.3304    8.8187;
    9.9986   11.7387    8.9473    9.5383   10.4294    9.4736;
   10.5262   11.6362    9.2696    9.8500   10.0660    9.0414;
   11.2665   11.3212    9.6290    9.6288   10.0228    9.3346]

[LSMEuro LSMAmer] = LSM_Illustrative(Spaths',K,r,q,T,NT,NS,PutCall,XmatrixHandle);

