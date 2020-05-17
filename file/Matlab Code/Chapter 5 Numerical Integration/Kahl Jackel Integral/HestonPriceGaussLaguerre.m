function y = HestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,trap,x,w)

% Heston (1993) call or put price by Gauss-Laguerre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   r = Risk free rate.
%   param = vector of Heston parameters
%     kappa  = mean reversion speed.
%     theta  = mean reversion level.
%     sigma  = volatility of vol
%     lambda = risk.
%     v0     = initial variance.
%     rho    = Heston parameter: correlation
%   trap:  1 = "Little Trap" formulation
%          0 = Original Heston formulation
%   x = Gauss Laguerre abscissas
%   w = Gauss Laguerre weights
% OUTPUT -------------------------------------------------------
%   The Heston call or put price

kappa  = param(1); 
theta  = param(2); 
sigma  = param(3); 
v0     = param(4); 
rho    = param(5);
lambda = 0;

% Numerical integration
for k=1:length(x);
	int1(k) = w(k)*HestonProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,1,trap);
	int2(k) = w(k)*HestonProb(x(k),kappa,theta,lambda,rho,sigma,T,K,S,rf,q,v0,2,trap);
end

% Define P1 and P2
P1 = 1/2 + 1/pi*sum(int1);
P2 = 1/2 + 1/pi*sum(int2);

% The call price
HestonC = S*exp(-q*T)*P1 - K*exp(-rf*T)*P2;

% The put price by put-call parity
HestonP = HestonC - S*exp(-q*T) + K*exp(-rf*T);

% Output the option price
if strcmp(PutCall,'C')
	y = HestonC;
else
	y = HestonP;
end

