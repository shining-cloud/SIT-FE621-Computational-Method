function y = d2P1dK2(phi,kappa,theta,lambda,rho,sigma,T,K,S,rf,v0)

% Returns the integrand for the second-order partial derivatives of P1 
% with respect to strike, namely d2P1/dK2.
% phi = integration variable
% Pnum = 1 or 2 (for the probabilities)
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v      = initial variance
% Option features.
%    T  = maturity
%    K  = strike price
%    S  = spot price
%    rf = risk free rate
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% Log of the stock price.
x = log(S);

% Parameters "u" and "b" for P1
u = 0.5;
b = kappa + lambda - rho*sigma;

d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

% Original Heston formulation.
G = (1 - g*exp(d*T))/(1-g);
C = rf*i*phi*T + kappa*sigma/sigma^2*((b - rho*sigma*i*phi + d)*T - 2*log(G));
D = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*T))/(1-g*exp(d*T)));

% The characteristic function
f1 = exp(C + D*v0 + i*phi*x);

% Return the real part of the integrand.
y = real((i*phi+1) * K^(-i*phi-2) * f1);

