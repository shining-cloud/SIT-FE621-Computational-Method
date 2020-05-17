function y = HestonProb(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v,Pnum,Trap);

% Returns the integrand for the risk neutral probabilities P1 and P2.
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
%    PutCall = 'C'all or 'P'ut
%    K = strike price
%    S = spot price
%    r = risk free rate
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% Log of the stock price.
x = log(S);

% Parameter "a" is the same for P1 and P2.
a = kappa*theta;

% Parameters "u" and "b" are different for P1 and P2.
if Pnum==1
	u = 0.5;
	b = kappa + lambda - rho*sigma;
else
	u = -0.5;
	b = kappa + lambda;
end

d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

if Trap==1
	% "Little Heston Trap" formulation
	c = 1/g;
	D = (b - rho*sigma*i*phi - d)/sigma^2*((1-exp(-d*tau))/(1-c*exp(-d*tau)));
	G = (1 - c*exp(-d*tau))/(1-c);
	C = (r-q)*i*phi*tau + a/sigma^2*((b - rho*sigma*i*phi - d)*tau - 2*log(G));
elseif Trap==0
	% Original Heston formulation.
	G = (1 - g*exp(d*tau))/(1-g);
	C = (r-q)*i*phi*tau + a/sigma^2*((b - rho*sigma*i*phi + d)*tau - 2*log(G));
	D = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));
end

% The characteristic function.
f = exp(C + D*v + i*phi*x);

% Return the real part of the integrand.
y = real(exp(-i*phi*log(K))*f/i/phi);

