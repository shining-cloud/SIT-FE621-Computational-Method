function y = HestonBivariateCF(phi1,phi2,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,trap);

% Returns the Heston bivariate characteristic function
% phi1,phi2 = integration variables
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

% Parameter "a" is the same for P1 and P2.
a = kappa*theta;

% Parameters "u" and "b" for P2.
u = -0.5;
b = kappa + lambda;

d = sqrt((rho*sigma*i*phi1 - b)^2 - sigma^2*(2*u*i*phi1 - phi1^2));
g = (b - rho*sigma*i*phi1 + d - sigma^2*i*phi2) ...
  / (b - rho*sigma*i*phi1 - d - sigma^2*i*phi2);
c = 1/g;
B = i*phi1;

if trap==1
	% Little Trap formulation in Kahl (2008)
	G = (c*exp(-d*tau)-1)/(c-1);
	A = (r-q)*i*phi1*tau + a/sigma^2*((b - rho*sigma*i*phi1 - d)*tau - 2*log(G));
	C = ((b - rho*sigma*i*phi1 - d) - (b - rho*sigma*i*phi1 + d)*c*exp(-d*tau)) ...
		/ sigma^2 / (1-c*exp(-d*tau));
elseif trap==0
	% Original Heston formulation.
	G = (1-g*exp(d*tau))/(1-g);
	A = (r-q)*i*phi1*tau + a/sigma^2*((b - rho*sigma*i*phi1 + d)*tau - 2*log(G));
	C = ((b - rho*sigma*i*phi1 + d) - (b - rho*sigma*i*phi1 - d)*g*exp(d*tau)) ...
		/ sigma^2 / (1-g*exp(d*tau));
end

% The characteristic function.
y = exp(A + B*x0 + C*v0);

