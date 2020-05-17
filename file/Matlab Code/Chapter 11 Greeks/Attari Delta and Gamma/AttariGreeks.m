function y = AttariGreeks(PutCall,kappa,theta,lambda,rho,sigma,T,K,S0,r,q,v0,trap,Greek,x,w)

% Returns the Greek for the Attari model
% PutCall = 'P'ut or 'C'all
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v      = initial variance
% Option features.
%    tau = maturity
%    K  = strike price
%    S0 = spot price
%    r  = risk free rate
%    q = dividend yield
%    Trap = 1 "Little Trap" formulation 
%           0  Original Heston formulation
% Greek = 'Delta' or 'Gamma'
% x,w = Integration abscissas and weights

for j=1:length(x);
    if strcmp(Greek,'Delta')
        A(j) = w(j) * AttariProbGreeks(x(j),kappa,theta,lambda,rho,sigma,T,K,S0,r,q,v0,trap,'Delta');
    elseif strcmp(Greek,'Gamma')
        A(j) = w(j) * AttariProbGreeks(x(j),kappa,theta,lambda,rho,sigma,T,K,S0,r,q,v0,trap,'Gamma');
    end
end
Integral = sum(A);

if strcmp(Greek,'Gamma')
    y = -exp(-r*T)*K/pi*Integral;
elseif strcmp(Greek,'Delta')
    if strcmp(PutCall,'C')
        y = 1 - exp(-r*T)*K/pi*Integral;
    else
        y = - exp(-r*T)*K/pi*Integral;
    end
end

