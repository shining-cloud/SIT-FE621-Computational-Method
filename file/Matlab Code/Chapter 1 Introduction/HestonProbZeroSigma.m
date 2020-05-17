function y = HestonProbZeroSigma(phi,kappa,theta,lambda,tau,K,S,r,q,Pnum);

% Returns the integrand for the risk neutral probabilities P1 and P2
% when sigma = 0.
% phi = integration variable
% Pnum = 1 or 2 (for the probabilities)
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

% Log of the stock price.
x = log(S);

% Parameter "a" is the same for P1 and P2.
a = kappa*theta;

% Parameters "u" and "b" are different for P1 and P2.
if Pnum==1
	u = 0.5;
	b = kappa + lambda;
else
	u = -0.5;
	b = kappa + lambda;
end

% D and C when sigma = 0
D = (u*i*phi - phi^2/2)*(1-exp(-b*tau))/b;
C = (r-q)*i*phi*tau + a*(u*i*phi-0.5*phi^2)/b*(tau-(1-exp(-b*tau))/b);

% The characteristic function.
f = exp(C + D*theta + i*phi*x);

% Return the real part of the integrand.
y = real(exp(-i*phi*log(K))*f/i/phi);

