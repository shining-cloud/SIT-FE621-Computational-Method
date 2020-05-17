function y = CZEuroCall(S0,tau,params,K,rf,q,x,w)

% Chiarella and Ziogas form of the Heston call
% Uses a generalization of the characteristic functions f1 and f2

% Set xi=0 and psi=0 in f1 and f2
xi  = 0;
psi = 0;

% Create the integrands
for k=1:length(x);
	phi = x(k);
	Int1(k) = w(k) * real(exp(-i*phi*log(K)) * CZCharFun(S0,tau,xi,params,K,rf,q,phi,psi,1)/(i*phi));
    Int2(k) = w(k) * real(exp(-i*phi*log(K)) * CZCharFun(S0,tau,xi,params,K,rf,q,phi,psi,2)/(i*phi));
end

% Define the probabilities
P1 = 1/2 + (1/pi)*sum(Int1);
P2 = 1/2 + (1/pi)*sum(Int2);

% The call price
y = S0*exp(-q*tau)*P1 - K*exp(-rf*tau)*P2;
