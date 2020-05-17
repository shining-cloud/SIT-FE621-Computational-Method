function y = DoubleHestonPriceGaussLaguerre(PutCall,S,K,T,rf,q,param,x,w,trap)

% Double Heston (2009) call or put price by Gauss-Laguerre Quadrature
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   rf = Risk free rate.
%   q  = Dividend yield
%   param = Two sets of Double Heston parameters
%           [kappa1 theta1 sigma1 v01 rho1, 
%            kappa2 theta2 sigma2 v02 rho2] 
%   x = Gauss Laguerre abscissas
%   w = Gauss Laguerre weights
%   trap:  1 = "Little Trap" formulation (Gauthier and Possamai
%          0 = Original Heston formulation (Christoffersen et al)
% OUTPUT -------------------------------------------------------
%   The Double Heston call or put price

for k=1:length(x);
	u = x(k);
	f2(k) = DoubleHestonCF(u  ,param,T,S,rf,q,trap);
	f1(k) = DoubleHestonCF(u-i,param,T,S,rf,q,trap);
	int2(k) = w(k)*real(exp(-i*u*log(K))*f2(k)/i/u);
	int1(k) = w(k)*real(exp(-i*u*log(K))*f1(k)/i/u/S/exp((rf-q)*T));
end

% The ITM probabilities
P1 = 1/2 + 1/pi*sum(int1);
P2 = 1/2 + 1/pi*sum(int2);

% The Double Heston call price
HestonC = S*exp(-q*T)*P1 - K*exp(-rf*T)*P2;

% The put price by put-call parity
HestonP = HestonC - S*exp(-q*T) + K*exp(-rf*T);

% Output the option price
if strcmp(PutCall,'C')
	y = HestonC;
else
	y = HestonP;
end

