function y = HestonPriceZeroSigma(PutCall,kappa,theta,lambda,tau,K,S,r,q,Uphi,dphi,Lphi)

% Heston (1993) price of a European option when sigma = 0.
% Uses the original formulation by Heston
% Heston parameters:
%    kappa  = volatility mean reversion speed parameter
%    theta  = volatility mean reversion level parameter
%    lambda = risk parameter
% Option features.
%    PutCall = 'C'all or 'P'ut
%    K = strike price
%    S = spot price
%    r = risk free rate
%    q = dividend yield

% Get the integrands for the two Heston ITM probabilities p1 and p2
phi = [Lphi:dphi:Uphi];
for x=1:length(phi)
	P1_int(x) = HestonProbZeroSigma(phi(x),kappa,theta,lambda,tau,K,S,r,q,1);
	P2_int(x) = HestonProbZeroSigma(phi(x),kappa,theta,lambda,tau,K,S,r,q,2);
end

% Integrate to get p1 and p2
p1 = 0.5 + (1 / pi) * trapz(P1_int)*dphi;
p2 = 0.5 + (1 / pi) * trapz(P2_int)*dphi;

% Restrict p1 and p2 to the interval [0,1]
p1 = max(min(1,p1),0);
p2 = max(min(1,p2),0);

% Heston Call price directly.
HestonC = S*p1*exp(-q*tau) - K*exp(-r*tau)*p2;

% Heston put price by put-call parity;
HestonP = HestonC + K*exp(-r*tau) - S*exp(-q*tau);

if strcmp(PutCall,'C')
    y = HestonC;
else
	y = HestonP;
end

