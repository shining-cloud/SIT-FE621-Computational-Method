function y = HestonLewisCallPrice(S,K,r,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w)
% Implements the Heston model using the Fundamental Transform described in
% the book by Alan Lewis (2000) "Option Valuation Under Stochastic
% Volatility: With Mathematica Code"

% INPUTS
%   S = Spot price
%   K = Strike price
%   r = Risk free rate
%   q = Dividend Yield
%   v0 = Heston parameter, initial variance
%   tau = Parameter
%   ki = Im(k) strip cutoff
%   theta = Heston parameter, mean reversion level
%   kappa = Heston parameter, mean reversion speed
%   sigma = Heston parameter, volatility of variance
%   form  = 2 produces C2(S,K,t) requires 0 < Im[k] < 1
%         = 1 produces C1(S,K,t) requires 1 < Im[k] < B
%   x = Gauss Laguerre rule: Abscissas
%   w = Gauss Laguerre rule: Weights

% Lewis Parameters
kmax = floor(max(1000,10/sqrt(v0*tau)));
X = log(S/K) + (r-q)*tau;

% Compute the integral
% Gauss Laguerre
for k = 1:length(x);
    u = x(k) + i*ki;
    int(k) = w(k) * LewisIntegrand(u,X,v0,tau,theta,kappa,sigma,rho);
end
Integral = sum(int);

% Compute the call price
if form==2
	% C2(S,K,t) Equation (2.10) page 41 of Lewis (2000)
	y = S*exp(-q*tau) - (1/pi)*K*exp(-r*tau)*Integral;
else
	% C1(S,K,t) Equation (2.8) page 40 of Lewis (2000)
	y = - (1/pi)*K*exp(-r*tau)*Integral;
end
