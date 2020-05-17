function y = LewisIntegrand311Greeks(u,kappa,theta,lambda,rho,sigma,tau,S,K,r,q,v0,Trap,Greek)

% Returns the integrand for the formulation in Equation (3.11) of Lewis (2001).
% u = integration variable
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0      = initial variance
% Option features.
%    K = strike price
%    S = spot price
%    r = risk free rate
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

% Change u to u - i/2
u = u - i/2;

% Log of the stock price.
x0 = log(S);

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
f2 = exp(C + D*v0 + i*u*x0);

% Derivative of f2 w.r.t. tau
dD = d*exp(d*tau)*(b-rho*sigma*u*i+d)*(g-1)/sigma^2/(1-g*exp(d*tau))^2;
dC = (r-q)*u*i + kappa*theta/sigma^2 ...
   * ((b-rho*sigma*u*i+d) + 2*g*d*exp(d*tau)/(1-g*exp(d*tau)));
dt = (dC + dD*v0);

% Derivative of f2 w.r.t. r
dr = i*u*tau;

% Change u - i/2 back to u
u = u + i/2;

% The integrand
if strcmp(Greek,'Price')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2);
elseif strcmp(Greek,'Delta')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*(i*u+1/2)/S);
elseif strcmp(Greek,'Gamma')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*(i*u+1/2)*(i*u-1/2)/S^2);
elseif strcmp(Greek,'Rho')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*(-tau+dr));
elseif strcmp(Greek,'Theta')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*(-r+dt));
elseif strcmp(Greek,'Vega1')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*D);
elseif strcmp(Greek,'Volga')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*(4*D*D*v0 + 2*D));
elseif strcmp(Greek,'Vanna')
    y = real(K^(-i*u)/(u^2+1/4)*exp(-r*tau)*f2*D*(i*u+1/2)/S);
end

