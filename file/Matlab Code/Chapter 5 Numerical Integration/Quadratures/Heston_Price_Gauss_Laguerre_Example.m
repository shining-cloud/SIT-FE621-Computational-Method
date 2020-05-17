% Heston (1993) Call price by 32-point Gauss-Laguerre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% By Fabrice Douglas Rouah

clc; clear;

% Weights and abscissas
n = 32;
[x w] = GenerateGaussLaguerre(n);

% Define the parameters and inputs
S = 100;         % Spot price.
K = 100;         % Strike
T = 1.5;         % Time to maturity.
r = 0.05;        % Risk free rate.
q = 0.01;        % Dividend yield
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.05;    % Heston parameter: mean reversion level.
sigma = 0.3;     % Heston parameter: volatility of vol
lambda = 0;      % Heston parameter: risk.
v0 = .05;        % Heston parameter: initial variance.
rho = 0.45;      % Heston parameter: correlation
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation

CallPrice = HestonPriceGaussLaguerre('C',S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);
PutPrice  = HestonPriceGaussLaguerre('P',S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);

fprintf('Using %1.0f-point Gauss Laguerre integration\n',n)
fprintf('--------------------------------------------\n')
fprintf('Call price %10.5f\n',CallPrice)
fprintf('Put price  %10.5f\n',PutPrice)
fprintf('\n')
fprintf('  Abscissas   Weights\n')
fprintf('----------------------\n')
disp([x w])

