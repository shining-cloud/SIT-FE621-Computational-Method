function y = HestonLVApprox(S,K,T,kappa,theta,sigma,v0,rho);

% Returns the Heston Local volatility in approximate form
% FOR RHO APPROXIMATELY 1 or -1
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance

% Modified parameters kappa' and theta'
kappa_ = kappa - rho*sigma/2;
theta_ = theta*kappa/kappa_;

% wT and vT'
xT = log(K/S);
wT = (v0-theta)*(1-exp(-kappa*T))/kappa + theta*T;
vT = (v0-theta_)*exp(-kappa_*T) + theta_;

% Integral
F1 = (v0-theta)/(kappa_-kappa);
E1 = exp((kappa_-kappa)*T) - 1;
F2 = theta/kappa_;
E2 = exp(kappa_*T) - 1;
Integral = exp(-kappa_*T)*(F1*E1 + F2*E2);

% Local Variance
uT = vT + rho*sigma*xT/wT*Integral;

% Local volatility
y = sqrt(uT);

% Gatheral formula (4.1) when rho = -1
% uT2 = vT - sigma*log(K/S)*(1-exp(-kappa_*T))/(kappa_*T);
% y2 = sqrt(uT2);

