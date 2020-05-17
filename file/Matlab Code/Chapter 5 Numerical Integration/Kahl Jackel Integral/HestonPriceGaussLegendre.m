function y = HestonPriceGaussLegendre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,a,b)

% Heston (1993) call or put  price by 32-point Gauss-Legendre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   r = Risk free rate.
%   kappa  = Heston parameter: mean reversion speed.
%   theta  = Heston parameter: mean reversion level.
%   sigma  = Heston parameter: volatility of vol
%   lambda = Heston parameter: risk.
%   v0     = Heston parameter: initial variance.
%   rho    = Heston parameter: correlation
%   trap:  1 = "Little Trap" formulation
%          0 = Original Heston formulation
%   x = Gauss Legendre abscissas
%   w = Gauss Legendre weights
%   a = Lower Limit of integral (a=0)
%   b = Upper Limit of integral (b=500 or something large)
% OUTPUT -------------------------------------------------------
%   The Heston call or put price

% Numerical integration
for k=1:length(x);
	X = (a+b)/2 + (b-a)/2*x(k);
	int1(k) = w(k)*HestonProb(X,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
	int2(k) = w(k)*HestonProb(X,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
end

% Define P1 and P2
P1 = 1/2 + 1/pi*sum(int1)*(b-a)/2;
P2 = 1/2 + 1/pi*sum(int2)*(b-a)/2;

% The call price
HestonC = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;

% The put price by put-call parity
HestonP = HestonC - S*exp(-q*T) + K*exp(-r*T);

% Output the option price
if strcmp(PutCall,'C')
	y = HestonC;
else
	y = HestonP;
end

