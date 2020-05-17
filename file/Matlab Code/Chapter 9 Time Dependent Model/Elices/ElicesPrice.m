function y = ElicesPrice(PutCall,S,K,T,r,q,param,paramfixed,v0,trap,x,w)

% Heston (1993) price using Elices time dependent model
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T  = Maturities (shortest first)
%   r = Risk free rate.
%   param = parameters (kappa,theta,sigma,v0,rho) (oldest first)
%   trap:  1 = "Little Trap" formulation
%          0 = Original Heston formulation
%   x = Gauss Laguerre abscissas
%   w = Gauss Laguerre weights
% OUTPUT -------------------------------------------------------
%   The Heston call or put price

% Maturities are in a vector, shortest maturities first
%    T = (T(1),T(2),...., T(N-1),T(N))
% Fixed parameters are in a matrix, those for the shortest maturity first
%    paramfixed: (paramfixed(1),paramfixed(2), ...,paramfixed(N-1))
% Dimension of parameter rows is one less than dimension of maturities

Mat = T(end);

% Numerical integration
for k=1:length(x);
	phi    = x(k);
	weight = w(k);
    f2(k) = ElicesCF(phi  ,param,paramfixed,v0,T,S,r,q,trap);
	f1(k) = ElicesCF(phi-i,param,paramfixed,v0,T,S,r,q,trap) / (S*exp((r-q)*Mat));
	int2(k) = weight * real(exp(-i*phi*log(K))*f2(k)/i/phi);
	int1(k) = weight * real(exp(-i*phi*log(K))*f1(k)/i/phi);
end


% Define P1 and P2
P1 = 1/2 + 1/pi*sum(int1);
P2 = 1/2 + 1/pi*sum(int2);

% The call price
HestonC = S*exp(-q*Mat)*P1 - K*exp(-r*Mat)*P2;

% The put price by put-call parity
HestonP = HestonC - S*exp(-q*Mat) + K*exp(-r*Mat);

% Output the option price
if strcmp(PutCall,'C')
	y = HestonC;
else
	y = HestonP;
end



























