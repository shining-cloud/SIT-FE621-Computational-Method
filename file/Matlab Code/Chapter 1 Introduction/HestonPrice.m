function y = HestonPrice(PutCall,kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,trap,Lphi,Uphi,dphi)

% Heston (1993) price of a European option.
% Uses the original formulation by Heston
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

% Build the integrands for P1 and P2;
for k=1:N;
	int1(k) = HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
	int2(k) = HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
end

% The integrals
I1 = trapz(int1)*dphi;
I2 = trapz(int2)*dphi;

% The probabilities P1 and P2
P1 = 1/2 + 1/pi*I1;
P2 = 1/2 + 1/pi*I2;

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

