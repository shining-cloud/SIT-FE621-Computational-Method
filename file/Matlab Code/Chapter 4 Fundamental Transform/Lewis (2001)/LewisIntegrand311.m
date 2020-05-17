function y = LewisIntegrand311(u,kappa,theta,lambda,rho,sigma,tau,S,K,r,q,v0,Trap)

% Returns the integrand for the formulation in Equation (3.11) of Lewis (2001).
% u = integration variable
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance
% Option features.
%    tau = maturity
%    K = strike price
%    S = spot price
%    r = risk free rate
%    q = dividend yield
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% Log of the stock price.
x = log(S);

% Heston a and b parameters
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

% The Heston characteristic function f2 for ln S(T)
CFlnST = exp(C + D*v0 + i*u*x);

% The cf for the Levy process XT
Y = log(S) + (r-q)*tau;
CFXT = exp(-i*u*Y)*CFlnST;
 
% The integrand in Equation (3.11) of Lewis (2001)
u = u + i/2;
W = Y - log(K);
y = real(exp(i*u*W)*CFXT/(u^2 + 1/4));



