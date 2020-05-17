function y = HestonPriceGaussLaguerre2CF(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,x,w,CF)

% Heston (1993) call or put price by Gauss-Laguerre Quadrature
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
%   x = Gauss Laguerre abscissas
%   w = Gauss Laguerre weights
%   CF = 1 Original Heston (2 c.f.)
%        2 Single c.f.
% OUTPUT -------------------------------------------------------
%   The Heston call or put price

% Numerical integration
for k=1:length(x);
	phi    = x(k);
	% Characteristic function 2 and its integrand
	f2(k) = HestonCF(phi,kappa,theta,lambda,rho,sigma,T,S,r,q,v0,2,trap);
	% Characteristic function 1 and its integrand
	if CF==1
		f1(k) = HestonCF(phi,kappa,theta,lambda,rho,sigma,T,S,r,q,v0,1,trap);
	elseif CF==2
	    f1(k) = HestonCF(phi-i,kappa,theta,lambda,rho,sigma,T,S,r,q,v0,2,trap)/(S*exp((r-q)*T));
	end
	int2(k) = w(k) * real(exp(-i*phi*log(K))*f2(k)/i/phi);
	int1(k) = w(k) * real(exp(-i*phi*log(K))*f1(k)/i/phi);
end

% Define P1 and P2
P1 = 1/2 + 1/pi*sum(int1);
P2 = 1/2 + 1/pi*sum(int2);

% The call price
HestonC = S*P1 - K*exp(-r*T)*P2;

% The put price by put-call parity
HestonP = HestonC - S + K*exp(-r*T);

% Output the option price
if strcmp(PutCall,'C')
	y = HestonC;
else
	y = HestonP;
end

