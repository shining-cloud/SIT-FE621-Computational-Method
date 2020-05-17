function y = MNProb(phi,param,param0,tau,tau0,K,S,rf,q,Pnum)

% Integrand for the piece-wise constant parameter Heston model of Nogel and
% Mikhailov.

% param  = Current vector of parameters fed into objective function
% param0 = matrix of old parameters (oldest row on top, newest row at the bottom)
% tau    = current maturity
% tau0   = vector of old maturities (oldest on top, newest at the bottom)
% K = strike price
% S = Spot price
% rf = risk free rate
% q = dividend yield
% Pnum = 1 or 2, depending on which characteristic function is to be used

x = log(S);
N = length(tau0);

C = 0;

% Create the past Cj and Dj values corresponding to the old maturities
for t=1:N
	kappa  = param0(t,1);
	theta  = param0(t,2);
	sigma  = param0(t,3);
	v0     = param0(t,4);
	rho    = param0(t,5);
	T      = tau0(t);
	if (t==1)
		D0 = 0;
		C0 = 0;
	else
		D0 = D;
		C0 = C;
    end
	C = Ct(phi,kappa,theta,sigma,rho,rf,q,T,Pnum,C0,D0);
	D = Dt(phi,kappa,theta,sigma,rho,rf,T,Pnum,C0,D0);
	clear kappa theta sigma v0 rho T D0 C0
end

% Cj and Dj values for a single maturity
kappa  = param(1);
theta  = param(2);
sigma  = param(3);
v0     = param(4);
rho    = param(5);
T      = tau;
if N==0
	D0 = 0;
	C0 = 0;
else
	D0 = D;
	C0 = C;
end
C = Ct(phi,kappa,theta,sigma,rho,rf,q,T,Pnum,C0,D0);
D = Dt(phi,kappa,theta,sigma,rho,rf,T,Pnum,C0,D0);

% The characteristic function.
f = exp(C + D*v0 + i*phi*x);

% Return the real part of the integrand.
y = real(exp(-i*phi*log(K))*f/i/phi);


