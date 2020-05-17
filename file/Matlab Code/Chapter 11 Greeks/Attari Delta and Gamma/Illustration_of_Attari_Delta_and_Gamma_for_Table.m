% Delta and Gamma using Attari (2004) formula
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.

clc; clear;

% Define the parameters and inputs
S = 100;         % Spot price.
K = 100;         % Strike price
T = 0.25;        % Time to maturity.
r = 0.05;        % Risk free rate.
q = 0.0;         % Dividend yield
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.1;     % Heston parameter: volatility of vol
v0    = .05;     % Heston parameter: initial variance.
rho   = -0.9;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
PutCall = 'C';

% Generate the abscissas and weights
[x w] = GenerateGaussLaguerre(32);


%% Analytic Greeks
Price = AttariPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
Delta = AttariGreeks(PutCall,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,'Delta',x,w);
Gamma = AttariGreeks(PutCall,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,'Gamma',x,w);


%% Output the results
fprintf('Attari for book table \n')
fprintf('---------------------- \n');
fprintf('Price  %10.4f\n',Price)
fprintf('Delta  %10.4f\n',Delta)
fprintf('Gamma  %10.4f\n',Gamma)
fprintf('---------------------- \n');


