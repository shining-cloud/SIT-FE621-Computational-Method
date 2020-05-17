% Heston (1993) Call price by 32-point Gauss-Legendre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% By Fabrice Douglas Rouah

clc; clear;

% Generate the Gauss Legendre weights
n = 32;
[x w] = GenerateGaussLegendre(n);

% Define the parameters and inputs
S = 100;         % Spot price.
K = 95;          % Strike
T = .25;         % Time to maturity.
r = 0.03;        % Risk free rate.
q = 0.02;        % Dividend yield
kappa = 4;       % Heston parameter: mean reversion speed.
theta = 0.09;    % Heston parameter: mean reversion level.
sigma = 0.1;     % Heston parameter: volatility of vol
lambda = 0;      % Heston parameter: risk.
v0 = .04;        % Heston parameter: initial variance.
rho = -0.7;      % Heston parameter: correlation
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
b = 1000;
a = 0;
CallPrice = HestonPriceGaussLegendre('C',S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,a,b);
PutPrice  = HestonPriceGaussLegendre('P',S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,a,b);

fprintf('Using %1.0f-point Gauss Legendre integration\n',n)
fprintf('--------------------------------------------\n')
fprintf('Call price %10.5f\n',CallPrice)
fprintf('Put price  %10.5f\n',PutPrice)
fprintf('\n')
fprintf('  Abscissas   Weights\n')
fprintf('----------------------\n')
disp([x w])
