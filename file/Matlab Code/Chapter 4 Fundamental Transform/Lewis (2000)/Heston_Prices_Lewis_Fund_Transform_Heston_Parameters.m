% Heston (1993) prices using the Fundamental transform described in Alan
% Lewis' book "Option Valuation Under Stochastic Volatility: With Mathematica Code"
% By Fabrice Douglas Rouah
% Finds the call price based on C1(S,K,t) and C2(S,K,t).

clc; clear;

%% Spot price, strike price, risk free rate, dividend yield, and maturity
S   = 100;
K   = 100;
rf  = 0.05;
q   = 0.01;
tau = 0.25;

% Heston model parameters.
kappa = 2;       % Volatility reversion speed
theta = 0.05;    % Volatility reversion level
sigma = 0.1;     % Volatility of variance
rho   = -0.9;    % Correlation
v0    = 0.05;    % Initial variance
lambda = 0;      % Risk

%% Integration settings
% Trapezoidal rule
a = 1e-50;  % Lower limit of integration range
b = 100;    % Upper limit of integration range
N = 10000;  % Number of integration points

% Gauss Laguerre 32-point rule
[x w] = GenerateGaussLaguerre(32);

%% Second expression C2(S,K,t) for the Lewis call price
ki = 0.5;
form = 2;
IntRule = 1;
C2Trapz = HestonLewisCallPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,IntRule,a,b,N,x,w);
IntRule = 2;
C2Gauss = HestonLewisCallPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,IntRule,a,b,N,x,w);

%% First expression C1(S,K,t) for the Lewis call price
ki = 1.5;
form = 1;
IntRule = 1;
C1Trapz = HestonLewisCallPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,IntRule,a,b,N,x,w);
IntRule = 2;
C1Gauss = HestonLewisCallPrice(S,K,rf,q,v0,tau,ki,theta,kappa,sigma,rho,form,IntRule,a,b,N,x,w);

%% Heston price with 32-point Gauss-Laguerre and "Little Trap" formulation
trap = 1;
CallGL = HestonCallGaussLaguerre(S,K,tau,rf,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w);


%% Output the results
clc;
fprintf('Comparison of Lewis (2000) Call prices\n');
fprintf('Form             Integration        Call Price\n');
fprintf('----------------------------------------------\n');
fprintf('Heston Form      Gauss Laguerre  %10.4f\n',CallGL);
fprintf('Lewis C1         Trapezoidal     %10.4f\n',C1Trapz);
fprintf('Lewis C1         Gauss Laguerre  %10.4f\n',C1Gauss);
fprintf('Lewis C2         Trapezoidal     %10.4f\n',C2Trapz);
fprintf('Lewis C2         Gauss Laguerre  %10.4f\n',C2Gauss);
fprintf('----------------------------------------------\n');


