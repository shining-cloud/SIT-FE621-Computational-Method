function y = HestonProbConsol(phi,kappa,theta,lambda,rho,sigma,tau,K,S,r,q,v0,Pnum,Trap);

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
%    q = dividend yield
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% Log of the stock price.
x = log(S);

% Parameter "a" is the same for P1 and P2.
a = kappa*theta;

% First characteristic function f1
u1 = 0.5;
b1 = kappa + lambda - rho*sigma;
d1 = sqrt((rho*sigma*i*phi - b1)^2 - sigma^2*(2*u1*i*phi - phi^2));
g1 = (b1 - rho*sigma*i*phi + d1) / (b1 - rho*sigma*i*phi - d1);

if Trap==1
	% "Little Heston Trap" formulation
	c1 = 1/g1;
	D1 = (b1 - rho*sigma*i*phi - d1)/sigma^2*((1-exp(-d1*tau))/(1-c1*exp(-d1*tau)));
	G1 = (1 - c1*exp(-d1*tau))/(1-c1);
	C1 = (r-q)*i*phi*tau + a/sigma^2*((b1 - rho*sigma*i*phi - d1)*tau - 2*log(G1));
elseif Trap==0
	% Original Heston formulation.
	G1 = (1 - g1*exp(d1*tau))/(1-g1);
	C1 = (r-q)*i*phi*tau + a/sigma^2*((b1 - rho*sigma*i*phi + d1)*tau - 2*log(G1));
	D1 = (b1 - rho*sigma*i*phi + d1)/sigma^2*((1-exp(d1*tau))/(1-g1*exp(d1*tau)));
end
f1 = exp(C1 + D1*v0 + i*phi*x);

% Second characteristic function f2
u2 = -0.5;
b2 = kappa + lambda;
d2 = sqrt((rho*sigma*i*phi - b2)^2 - sigma^2*(2*u2*i*phi - phi^2));
g2 = (b2 - rho*sigma*i*phi + d2) / (b2 - rho*sigma*i*phi - d2);

if Trap==1
	% "Little Heston Trap" formulation
	c2 = 1/g2;
	D2 = (b2 - rho*sigma*i*phi - d2)/sigma^2*((1-exp(-d2*tau))/(1-c2*exp(-d2*tau)));
	G2 = (1 - c2*exp(-d2*tau))/(1-c2);
	C2 = (r-q)*i*phi*tau + a/sigma^2*((b2 - rho*sigma*i*phi - d2)*tau - 2*log(G2));
elseif Trap==0
	% Original Heston formulation.
	G2 = (1 - g2*exp(d2*tau))/(1-g2);
	C2 = (r-q)*i*phi*tau + a/sigma^2*((b2 - rho*sigma*i*phi + d2)*tau - 2*log(G2));
	D2 = (b2 - rho*sigma*i*phi + d2)/sigma^2*((1-exp(d2*tau))/(1-g2*exp(d2*tau)));
end
f2 = exp(C2 + D2*v0 + i*phi*x);

% Return the real part of the integrand.
y = real(exp(-i*phi*log(K))/i/phi*(S*exp(-q*tau)*f1 - K*exp(-r*tau)*f2));

