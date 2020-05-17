% Evaluation of American put options fromt the Heston model.
% Benchmark against values in Clarke and Parrot (1999)

clc; clear;

% Import the stock prices and true values
S = [8 9 10 11 12];
TruePrice = [2.0000, 1.107641, 0.520030, 0.213668, 0.082036];

% Settings from Clarke and Parrot (1999)
K = 10;
kappa = 5;
theta = 0.16;
sigma = 0.9;
v0 = 0.0625;
rho = 0.1;
lambda = 0;
T = 1/4;
r = 0.1;
q = 0.0;
params = [kappa theta sigma v0 rho lambda];

% Settings for the European prices
PutCall = 'P';
[x w] = GenerateGaussLaguerre(32);
trap = 1;

% Settings for the LSM algorithm
% Number of time steps and number of stock price paths
NT =  1000;
NS = 50000;

% Design matrix for the LSM algorithm
XmatrixHandle = {@(x)ones(length(x),1), @(x)(1-x),@(x)1./2.*(2-4.*x-x.^2)};

% Matrices for the correlated N(0,1) random variables
Zv = randn(NT,NS);
Zs = rho.*Zv + sqrt(1-rho.^2).*randn(NT,NS);

%% Loop through each stock price, calculating the LSM prices and Greeks
for s=1:length(S)
    [EuroPrice(s) AmerPrice(s)] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'price');
    [EuroDelta(s) AmerDelta(s)] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'delta');
    [EuroGamma(s) AmerGamma(s)] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'gamma');
    [EuroTheta(s) AmerTheta(s)] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'theta');
    [EuroRho(s)   AmerRho(s)  ] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'rho');
    [EuroVega1(s) AmerVega1(s)] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'vega1');
    [EuroVanna(s) AmerVanna(s)] = LSMGreeks(S(s),K,params,T,r,q,NT,NS,Zv,Zs,PutCall,XmatrixHandle,'vanna');
end

%% Loop through each stock price, calculating the Closed form European prices and Greeks
for s=1:length(S)
    % Generate the paths for the stock price
    ClosedPrice(s) = HestonPriceGaussLaguerre(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
    ClosedDelta(s) = HestonGreeks(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Delta');
    ClosedGamma(s) = HestonGreeks(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Gamma');
    ClosedTheta(s) = HestonGreeks(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Theta');
    ClosedRho(s)   = HestonGreeks(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Rho');
    ClosedVega1(s) = HestonGreeks(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vega1');
    ClosedVanna(s) = HestonGreeks(PutCall,S(s),K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,'Vanna');
end


%% Display the results for LSM American prices
clc;
fprintf('-------------------------------------------------------------------- \n')
fprintf('Clarke and Parrott American Greeks with Least-Squares Monte-Carlo \n')
fprintf(['LSM uses ' num2str(NT) ' time steps, and ' num2str(NS) ' stock paths  \n'])
fprintf('-------------------------------------------------------------------- \n')
fprintf(' S(0)    Price   Delta    Gamma    Theta     Rho     Vega1    Vanna \n')
fprintf('-------------------------------------------------------------------- \n')
for s=1:length(S)
    fprintf('%3.0f %10.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
        S(s),AmerPrice(s),AmerDelta(s),AmerGamma(s),AmerTheta(s),AmerRho(s),AmerVega1(s),AmerVanna(s));
end
fprintf('-------------------------------------------------------------------- \n')
fprintf('\n')

%% Display the results for European puts
fprintf('LSM European prices and Greeks \n');
fprintf('-------------------------------------------------------------------- \n')
fprintf(' S(0)    Price   Delta    Gamma    Theta     Rho     Vega1    Vanna \n')
fprintf('-------------------------------------------------------------------- \n')
for s=1:length(S)
    fprintf('%3.0f %10.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
        S(s),EuroPrice(s),EuroDelta(s),EuroGamma(s),EuroTheta(s),EuroRho(s),EuroVega1(s),EuroVanna(s));
end
fprintf('-------------------------------------------------------------------- \n')
fprintf('\n')

%% Display the results for Closed-form European puts
fprintf('Closed-form European prices and Greeks \n');
fprintf('-------------------------------------------------------------------- \n')
fprintf(' S(0)    Price   Delta    Gamma    Theta     Rho     Vega1    Vanna \n')
fprintf('-------------------------------------------------------------------- \n')
for s=1:length(S)
    fprintf('%3.0f %10.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n',...
        S(s),ClosedPrice(s),ClosedDelta(s),ClosedGamma(s),ClosedTheta(s),ClosedRho(s),ClosedVega1(s),ClosedVanna(s));
end
fprintf('-------------------------------------------------------------------- \n')

