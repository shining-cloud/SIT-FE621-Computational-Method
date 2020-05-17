% Evaluation of American put options fromt the Heston model.
% Benchmark against values in Clarke and Parrot (1999)

clc; clear;

% Import the stock prices and true values
S = [8 9 10 11 12];
TruePrice = [2 1.107641 0.520030 0.213668 0.082036];

% Settings from Clarke and Parrot (1999)
K = 10;
kappa = 5;
theta = 0.16;
sigma = 0.9;
V0 = 0.0625;
rho = 0.1;
lambda = 0;
T = 1/4;
rf = 0.1;
q  = 0.0;
param = [kappa theta sigma V0 rho lambda];

% Settings for the European prices
PutCall = 'P';
trap = 1;
[x w] = GenerateGaussLaguerre(32);

% Settings for the Belieava-Nawalkha tree
NT = 10;
threshold = 0.01;

%% Loop through each stock price, calculating the prices
for s=1:length(S)
    tic;
    [EuroTree(s) BNAmerTree(s) Euro Amer Yt V X Prob Branch] = BuildBivariateTree3(S(s),PutCall,K,T,rf,NT,kappa,theta,sigma,V0,rho,threshold);
    time = toc;
    EuroClosed = HestonPriceGaussLaguerre(PutCall,S(s),K,T,rf,q,param,trap,x,w);
    AmerCV(s) = BNAmerTree(s) + (EuroClosed - EuroTree(s));
    TreeError(s) = TruePrice(s) - BNAmerTree(s);
    CVError(s)   = TruePrice(s) - AmerCV(s);
    fprintf('Obtained option price %1.0f in %5.4f seconds\n',s,time);
end

%% Output the results
fprintf('----------------------------------------------------------------------\n');
fprintf('Reproduction of Clarke and Parrott (1999) prices\n');
fprintf('Using %1.0f time steps\n',NT);
fprintf('----------------------------------------------------------------------\n');
fprintf('  Spot    TruePrice   BNTreePrice   TreeError   CVPrice      CVError \n');
fprintf('----------------------------------------------------------------------\n');
for s=1:length(S)
    fprintf('%5.0f %12.6f %12.6f %10.3f %12.6f %10.3f\n',...
        S(s),TruePrice(s),BNAmerTree(s),TreeError(s),AmerCV(s),CVError(s));
end
fprintf('----------------------------------------------------------------------\n');

