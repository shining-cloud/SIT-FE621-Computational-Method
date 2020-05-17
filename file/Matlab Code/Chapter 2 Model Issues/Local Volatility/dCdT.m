function y = dCdT(kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,x,w);

% Returns the partial derivative of the call price
% with respect to maturity, namely dC/dT.
% Integration by Gauss Laguerre
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
%    T  = maturity
%    K  = strike price
%    S  = spot price
%    rf = risk free rate
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation
% Integration inputs
%    x = Gauss Laguerre abscissas
%    w = Gauss Laguerre weights

for k=1:length(x)
	int1(k) = w(k)*dPjdT(1,x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0);
	int2(k) = w(k)*dPjdT(2,x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0);
	int3(k) = w(k)*HestonProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,2,0);
end


dP1dT = (1/pi)*sum(int1);
dP2dT = (1/pi)*sum(int2);
P2    = (1/2) + (1/pi)*sum(int3);

y = S*dP1dT - K*exp(-rf*T)*(-rf*P2 + dP2dT);
