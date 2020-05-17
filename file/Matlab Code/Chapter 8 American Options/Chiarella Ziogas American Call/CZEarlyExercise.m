function V = CZEarlyExercise(S0,tau,params,K,rf,q,xs,ws,xt,wt,Nt,b0,b1,a,b,c,d,DoubleType)

% Chiarella and Ziogas early exercise premium for Heston American calls
% params = Heston parameters
% Option features.
%   S0  = spot price
%   tau = maturity
%   params = Heston parameters
%   K   = strike price
%   rf  = risk free rate
%   q   = dividend yield
% Integration features
%   xs = abscissas for Gauss-Laguerre
%   ws = weights for Gauss-Laguerre
%   xt = abscissas for Gauss-Legendre
%   wt = weights for Gauss-Legendre
%   Nt = number of points for double trapezoidal rule for early exercise
%   b0 = Chiarella and Ziogas parameter
%   b1 = Chiarella and Ziogas parameter
%   integration grid (t,x) on (a,b) x (c,d)
%   DoubleType = type of double integration for early exercise

% The integrals
if strcmp(DoubleType,'GLe')
    Int1 = DoubleGaussLegendre(S0,tau,params,K,rf,q,b0,b1,xt,wt,xt,wt,a,b,c,d,1);
    Int2 = DoubleGaussLegendre(S0,tau,params,K,rf,q,b0,b1,xt,wt,xt,wt,a,b,c,d,2);
elseif strcmp(DoubleType,'Trapz')
    ht = (b-a)/Nt;
    hs = (d-c)/Nt;
    X = zeros(Nt,1);
    T = zeros(Nt,1);
    for j=1:Nt+1;
        T(j) = a + (j-1)*ht;
        X(j) = c + (j-1)*hs;
    end
    Int1 = DoubleTrapezoidal(params,S0,K,tau,rf,q,b0,b1,X,T,1);
    Int2 = DoubleTrapezoidal(params,S0,K,tau,rf,q,b0,b1,X,T,2);
end

% The early exercise premium
V1 = S0*(1-exp(-q*tau))/2 + (1/pi)*S0*q*exp(-q*tau)*Int1;
V2 = K*(1-exp(-rf*tau))/2 + (1/pi)*K*rf*exp(-rf*tau)*Int2;

V = V1 - V2;

