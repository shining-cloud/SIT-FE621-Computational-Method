% Heston (1993) Call or Put price by 32-point Gauss-Laguerre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% Uses the original Heston integrand or the Carr-Madan FFT integrand
% By Fabrice Douglas Rouah

clc; clear;

% Weights and abscissas
[x w] = GenerateGaussLaguerre(32);

% Define the parameters and inputs
S = 100;         % Spot price.
K = 100;         % Strike
T = 0.5;         % Time to maturity.
r = 0.10;        % Risk free rate.
q = 0.07;        % Dividend yield
kappa = 2;       % Heston parameter: mean reversion speed.
theta = 0.06;    % Heston parameter: mean reversion level.
sigma = 0.1;     % Heston parameter: volatility of vol
v0    = 0.06;    % Heston parameter: initial variance.
rho   = -0.7;    % Heston parameter: correlation
lambda = 0;      % Heston parameter: risk.
trap = 1;        % 1 = "Little Trap" formulation
                 % 0 = Original Heston formulation
alpha = 1.75;    % Carr-Madam dampening factor
PutCall = 'P';   % 'P'ut or 'C'all

% The price using the original Heston integrand
PriceHeston = HestonPriceGaussLaguerre('Heston',alpha,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall);

% The price using the Carr-Madan integrand with damping factor alpha
S = S*exp(-q*T);
PriceCarrMadan = HestonPriceGaussLaguerre('CarrMadan',alpha,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,PutCall);

% Display the results
fprintf('Method                    Price\n');
fprintf('------------------------------------\n');
fprintf('Heston integrand      %10.4f \n',PriceHeston);
fprintf('Carr-Madan integrand  %10.4f \n',PriceCarrMadan);
fprintf('------------------------------------\n');



