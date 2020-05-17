function A = AttariProb(u,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Trap)

% Returns the integrand for the Attari integrand
% u = integration variable
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

% Parameters a, u, and b for second Heston c.f. (f2)
a = kappa*theta;
b = kappa + lambda;

d = sqrt((rho*sigma*i*u - b)^2 + sigma^2*(i*u + u^2));
g = (b - rho*sigma*i*u + d) / (b - rho*sigma*i*u - d);

if Trap==1
	% "Little Heston Trap" formulation
	c = 1/g;
	D = (b - rho*sigma*i*u - d)/sigma^2*((1-exp(-d*tau))/(1-c*exp(-d*tau)));
	G = (1 - c*exp(-d*tau))/(1-c);
	C = (r-q)*i*u*tau + a/sigma^2*((b - rho*sigma*i*u - d)*tau - 2*log(G));
elseif Trap==0
	% Original Heston formulation.
	G = (1 - g*exp(d*tau))/(1-g);
	C = (r-q)*i*u*tau + a/sigma^2*((b - rho*sigma*i*u + d)*tau - 2*log(G));
	D = (b - rho*sigma*i*u + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));
end

% The characteristic function for Attari.
f = exp(C + D*v0 - i*u*r*tau);

L = log(exp(-r*tau)*K/S);

% The Attari (2004) integrand
y = (real(f) + imag(f)/u)*cos(L*u) + (imag(f) - real(f)/u)*sin(L*u);
A = y/(1+u^2);

