function y = HestonLewisGreekPrice(S,K,r,q,v0,tau,ki,theta,kappa,sigma,rho,form,x,w,Greek)
% Implements the Greeks on the Fundamental Transform described in
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
%   Greek = 'Price','Delta','Gamma','Rho','Theta','Vega1','Vanna', 'Volga'

% Lewis Parameters
kmax = floor(max(1000,10/sqrt(v0*tau)));
X = log(S/K) + (r-q)*tau;

% Compute the integral
% IntRule = 2 32-point Gauss Laguerre
% Gauss Laguerre
for k = 1:length(x);
    u = x(k) + i*ki;
    int(k) = w(k) * LewisIntegrandGreek(u,S,r,q,X,v0,tau,theta,kappa,sigma,rho,Greek);
end
Integral = sum(int);

% Compute the Greeks
switch Greek
    case 'Price'
    if form==2
        % Price based on C2(K) Equation (2.10) page 41 of Lewis (2000)
        y = S*exp(-q*tau) - (K/pi)*Integral;
    else
        % Price based on C1(K) Equatio+n (2.8) page 40 of Lewis (2000)
        y = - (K/pi)*Integral;
    end
    case 'Delta'
        if form==2
            % Delta based on C2(K) Equation (2.10) page 41 of Lewis (2000)
            y = exp(-q*tau) - (K/pi)*Integral;
        else
            % Delta based on C1(K) Equation (2.8) page 40 of Lewis (2000)
            y = - (K/pi)*Integral;
        end
    case 'Theta'
        if form==2
            y = q*S*exp(-q*tau) + (K/pi)*Integral;
        else
            y = (K/pi)*Integral;
        end
    otherwise
        % Remaining Greeks based on C1(K) or C2(K) have the same integrand
        y = - (K/pi)*Integral;
end

