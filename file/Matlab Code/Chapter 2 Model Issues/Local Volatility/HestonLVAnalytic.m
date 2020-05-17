function y = HestonLVAnalytic(S,K,T,rf,kappa,theta,sigma,lambda,v0,rho,x,w,Trap);

% Returns the Heston Local volatility in analytic form
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
%    T = maturity
%    K = strike price
%    S = spot price
%    r = risk free rate
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation

for k=1:length(x)
	int1(k) = w(k)*dPjdT(1,x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,Trap);
	int2(k) = w(k)*dPjdT(2,x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,Trap);
	int3(k) = w(k)*HestonProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0,2,Trap);
	int4(k) = w(k)*d2P1dK2(x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0);
	int5(k) = w(k)*dP2dK2_2(x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,v0);
end

dP1dT = (1/pi)*sum(int1);
dP2dT = (1/pi)*sum(int2);
P2    = (1/2) + (1/pi)*sum(int3);

% dC/dT : derivative of call with respect to T
dCdT = S*dP1dT - K*exp(-rf*T)*(-rf*P2 + dP2dT);

dP1dK2 = (1/pi)*sum(int4);
TwodP2dK2 = (1/pi)*sum(int5);

% d2C/dK2 : second derivative of call with respect to K^2
d2CdK2 = S*dP1dK2 - exp(-rf*T)*TwodP2dK2;

% Local Variance
LocalVar = 2*dCdT / (K^2*d2CdK2);

% Local Volatility
y = sqrt(LocalVar);
