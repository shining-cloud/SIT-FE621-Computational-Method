function y = HestonCF(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);

% Returns the integrand for the risk neutral probabilities P1 and P2.
% phi = integration variable
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance
% Option features.
%    PutCall = 'C'all or 'P'ut
%    K = strike price
%    S = spot price
%    r = risk free rate

% Log of the stock price.
x0 = log(S);

% Parameter "a"
a = kappa*theta;

% Parameters for P2
u = -0.5;
b = kappa + lambda;
d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

if trap==1
	% "Little Heston Trap" formulation
	c = 1/g;
	G = (1 - c*exp(-d*tau))/(1-c);
	D = (b - rho*sigma*i*phi - d)/sigma^2*((1-exp(-d*tau))/(1-c*exp(-d*tau)));
	C = (r-q)*i*phi*tau + a/sigma^2*((b - rho*sigma*i*phi - d)*tau - 2*log(G));
elseif trap==0
	% Original Heston formulation
	G = (1 - g*exp(d*tau))/(1-g);
	C = (r-q)*i*phi*tau + a/sigma^2*((b - rho*sigma*i*phi + d)*tau - 2*log(G));
	D = (b - rho*sigma*i*phi + d)/sigma^2*((1-exp(d*tau))/(1-g*exp(d*tau)));
end

% The characteristic function.
y = exp(C + D*v0 + i*phi*x0);

