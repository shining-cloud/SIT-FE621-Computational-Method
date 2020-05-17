function A = AttariProbGreeks(u,kappa,theta,lambda,rho,sigma,tau,K,S0,r,q,v0,Trap,Greek)

% Returns the integrand for the Attari integrand for the Greeks
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
%    K  = strike price
%    S0 = spot price
%    r  = risk free rate
%    q = dividend yield
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation
% Greek = 'Delta' or 'Gamma'

% Parameters a and b for second Heston c.f. (f2)
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

% The L function
L  = log(exp(-r*tau)*K/S0);

% The coefficients
F = (real(f) + imag(f)/u);
G = (imag(f) - real(f)/u);

% Return the integrand for the chosen Greek
if strcmp(Greek,'Delta');
    A = (F*sin(L*u) - G*cos(L*u)) * u/S0/(1+u^2);
elseif strcmp(Greek,'Gamma');
    y = (-F*(cos(L*u)*u + sin(L*u)) + G*(-sin(L*u)*u + cos(L*u)));
    A = y * u/S0^2/(1+u^2);
end


