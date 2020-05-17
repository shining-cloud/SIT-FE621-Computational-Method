% Heston (1993) prices using the characteristic function described in Alan
% Lewis' paper "A Simple Option Formula for General Jump-Diffusion and other
% Exponential Levy Processes"
% Finds the call price in Equation (3.11) of his paper

clc; clear;

% Spot price, strike price, risk free rate, dividend yield, and maturity
S   = 100;
K   = 100;
r   = 0.05;
q   = 0.01;
T   = 0.25;

% Heston model parameters.
kappa  = 2;       % Volatility reversion speed
theta  = 0.05;    % Volatility reversion level
sigma  = 0.1;     % Volatility of variance
rho    = -0.9;    % Correlation
v0     = 0.05;    % Initial variance
lambda = 0;       % Vol risk parameter
trap   = 1;       % Little trap formulation
param = [kappa theta sigma v0 rho];

% Absicissas and weights
[x w] = GenerateGaussLaguerre(32);
trap = 1;
K = [95:105];

% Find the Lewis (2001) prices from his Equation (3.11) and the orignal
% and the original Heston (1993) prices
fprintf('Lewis (2001) and Heston (1993) call prices\n')
fprintf('---------------------------------------------\n')
fprintf('Strike   LewisPrice   HestonPrice     Error\n')
fprintf('---------------------------------------------\n')
for k=1:length(K)
	% Lewis (2001) price using his Equation 3.11
	LewisPrice(k) = LewisPrice311(S,K(k),r,q,T,theta,kappa,sigma,rho,v0,trap,x,w);
	% Find the Heston price using Simpsons' 3/8 rule and the "Little Heston Trap" formulation
    HestonPrice(k) = HestonPriceGaussLaguerre('C',S,K(k),T,r,q,param,trap,x,w);
	% Find the error
	e(k) = (HestonPrice(k)-LewisPrice(k))/HestonPrice(k)*100;
	% Display the results
    fprintf(' %3.0f %12.4f %12.4f  %12.4f\n', K(k),LewisPrice(k),HestonPrice(k),e(k));
end
fprintf('---------------------------------------------\n')

