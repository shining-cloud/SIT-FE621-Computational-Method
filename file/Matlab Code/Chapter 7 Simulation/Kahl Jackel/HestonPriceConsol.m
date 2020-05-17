function y = HestonPriceConsol(PutCall,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,Lphi,Uphi,dphi)

% Heston (1993) price of a European option.
% Uses the consolidated form of the integrand
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
%    rho    = correlation between two Brownian motions
%    sigma  = volatility of variance
%    v0     = initial variance
% Option features.
%    PutCall = 'C'all or 'P'ut
%    K = strike price
%    S = spot price
%    r = risk free rate
%    q = dividend yield
%    T = maturity
% Integration features
%    L = lower limit
%    U = upper limit
%    dphi = integration increment

% Build the integration grid
phi = [Lphi:dphi:Uphi];
N = length(phi);

% Build the consolidated integrand
for k=1:N;
	inte(k) = HestonProbConsol(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
end

% The integral
I = trapz(inte)*dphi;

% The call price
HestonC = (1/2)*S*exp(-q*T) - (1/2)*K*exp(-r*T) + I/pi;

% The put price by put-call parity
HestonP = HestonC - S*exp(-q*T) + K*exp(-r*T);

% Output the option price
if strcmp(PutCall,'C')
	y = HestonC;
else
	y = HestonP;
end

