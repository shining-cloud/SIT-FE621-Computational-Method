function y = ElicesCF(phi1,param,paramfixed,v0,tau,S,r,q,trap)

% Returns the Bivariate Heston characteristic function
% Time dependent version
%    phi1 = integration variables
%    T  = Maturities (shortest first)
%    param = current parameters (kappa,theta,sigma,v0,rho)
%    paramfixed = past period parameters, increasing maturity order
% Option features.
%    S = spot price
%    r = risk free rate
%    trap  1 = Little trap formulation
%          0 = Original Heston formulation
% Fixed parameters are in order of increasing maturity
% [param(T1); param(T2);...;param(TN)]
% Dimension of parameter rows is one less than dimension of maturities

% Log of the stock price.
x0 = log(S);

% Maturity increments
N = length(tau);

% phi2(N) coefficients vector.  
% Last maturity at 2nd-to-last position in vector.
phi2(N) = 0 + 0i;
if N>=2
    phi2(N-1) = -i*C(phi1,phi2(N),param,tau(N),S,trap);
end
if N>=3
    for t=N-2:-1:1
        phi2(t) = -i*C(phi1,phi2(t+1),paramfixed(t+1,:),tau(t+1),S,trap);
    end
end
if N>=2
    phi20 = -i*C(phi1,phi2(1),paramfixed(1,:),tau(1),S,trap);
else
    phi20 = -i*C(phi1,phi2(1),param,tau(1),S,trap);
end    
    
% A coefficients.  
% Last maturity at last position in vector.
Ah(N) = A(phi1,0,param,tau(N),r,q,trap);  % Current params
if N>=2
    for t=N-1:-1:1
        Ah(t) = A(phi1,phi2(t),paramfixed(t,:),tau(t),r,q,trap);
    end
end

% C coefficient
if N>=2
    Ch = C(phi1,phi20,paramfixed(1,:),tau(1),S,trap);
else
    Ch = C(phi1,phi20,param,tau(1),S,trap);
end

y = exp(sum(Ah) + i*phi1*x0 + Ch*v0);

