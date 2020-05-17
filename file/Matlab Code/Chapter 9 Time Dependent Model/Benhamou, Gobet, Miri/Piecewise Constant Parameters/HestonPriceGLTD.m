function y = HestonPriceGLTD(param,param0,tau,tau0,K,S,PutCall,rf,q,x,w)

% Option price for the Heston model with piece-wise constant parameters,
% from Nogel and Mikhailov.
% Uses 32-point Gauss Laguerre integration

% param  = Current vector of parameters fed into objective function
% param0 = matrix of old parameters (oldest row on top, newest row at the bottom)
% tau    = current maturity
% tau0   = vector of old maturities (oldest on top, newest at the bottom)
% K = strike price
% S = spot price
% PutCall = 'P'ut or 'C'all
% rf = risk free rate
% q = dividend yield
% x = vector of abscissas for the Gauss Laguerre integral
% w = vector of weights for the Gauss Laguerre integral

N = length(x);

% Numerical integration
for k=1:N;
	int1(k) = w(k)*HestonProbTD(x(k),param,param0,tau,tau0,K,S,rf,q,1);
	int2(k) = w(k)*HestonProbTD(x(k),param,param0,tau,tau0,K,S,rf,q,2);
end

% Define P1 and P2
P1 = 1/2 + 1/pi*sum(int1);
P2 = 1/2 + 1/pi*sum(int2);

% The Call price
Call = S*exp(-q*tau)*P1 - K*exp(-rf*tau)*P2;

% The option price
if strcmp(PutCall,'C')
	y = Call;
else
	y = Call - S*exp(-q*tau) + K*exp(-rf*tau);
end

