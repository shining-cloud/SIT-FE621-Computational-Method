function y = HestonPriceTrapezoidal(S,K,T,rf,q,kappa,theta,sigma,lambda,v0,rho,trap,PutCall,dphi,Uphi,Lphi)

% Heston (1993) option price by Trapezoidal rule
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% By Fabrice Douglas Rouah
% Returns the call price
% INPUTS -------------------------------------------------------
% S = Spot price.
% K = Strike
% T = Time to maturity.
% r = Risk free rate.
% kappa = Heston parameter: mean reversion speed.
% theta = Heston parameter: mean reversion level.
% sigma = Heston parameter: volatility of vol
% lambda = Heston parameter: risk.
% v0 = Heston parameter: initial variance.
% rho = Heston parameter: correlation
% trap:  1 = "Little Trap" formulation
%        0 = Original Heston formulation
% dphi = integration range increment
% Uphi = integration range upper limit
% Lphi = integration range lower limit

% OUTPUT -------------------------------------------------------
% The Heston call price

% Built the integration grid
PHI = [Lphi:dphi:Uphi];
N = length(PHI);
h = (PHI(N) - PHI(1))/N;

% Trapezoidal rule --------------------------------------------------
% Weights
wt = [1/2 ones(1,N-2) 1/2];
for k=1:N;
	phi = PHI(k);
	int1(k) = h*wt(k)*HestonProb(phi,kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,1,trap);
	int2(k) = h*wt(k)*HestonProb(phi,kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,2,trap);
end

% The integrals
I1 = sum(int1);
I2 = sum(int2);

% The probabilities P1 and P2
P1 = 1/2 + 1/pi*I1;
P2 = 1/2 + 1/pi*I2;

% The Call price
CallPrice = S*exp(-q*T)*P1 - K*exp(-rf*T)*P2;

% Output the option price
if strcmp(PutCall,'C')
	y = CallPrice;
else
	y = CallPrice - S*exp(-q*T) + K*exp(-rf*T);
end

