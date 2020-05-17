function y = HestonPriceNewtonCoates(PutCall,S,K,T,r,q,kappa,theta,sigma,lambda,v0,rho,trap,method,a,b,N)

% Heston (1993) Call price by Newton-Coates Formulas
% Uses the original Heston formulation of the characteristic function,
% or the "Little Heston Trap" formulation of Albrecher et al.
% By Fabrice Douglas Rouah
% Returns the call price
% INPUTS -------------------------------------------------------
%   PutCall = 'C' Call or 'P' Put
%   S = Spot price.
%   K = Strike
%   T = Time to maturity.
%   r = Risk free rate.
%   kappa = Heston parameter: mean reversion speed.
%   theta = Heston parameter: mean reversion level.
%   sigma = Heston parameter: volatility of vol
%   lambda = Heston parameter: risk.
%   v0 = Heston parameter: initial variance.
%   rho = Heston parameter: correlation
%   trap:  1 = "Little Trap" formulation
%          0 = Original Heston formulation
%   method 1 = mid-point rule
%          2 = trapezoidal rule
%          3 = Simpson's rule
%          4 = Simpson's 3/8 rule
%   b = Upper limit for Newton-Cotes
%   h = Increment for Newton-Cotes
%   a = Lower limit for Newton-Cotes

% OUTPUT -------------------------------------------------------
%   The Heston call or put price

% Built the integration grid
% For Simpson's 3/8 rule, the code ensures that N-1 is divisible by 3
h = (b-a)/(N-1);
phi = [a:h:b];

% Integration methods
if method==1
	% Mid-Point rule --------------------------------------
	wt = h.*ones(1,N);
	for k=1:N-1;
		mid(k)  = (phi(k) + phi(k+1))/2;
		int1(k) = wt(k)*HestonProb(mid(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
		int2(k) = wt(k)*HestonProb(mid(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
	end
elseif method==2
	% Trapezoidal rule --------------------------------------------------
	% Weights
	wt = h.*[1/2 ones(1,N-2) 1/2];
	for k=1:N;
		int1(k) = wt(k)*HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
		int2(k) = wt(k)*HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
	end
elseif method==3
	% Simpson's Rule ----------------------------------------------------
	% Weights
	wt = (h/3).*[1 (3+(-1).^[2:N-1]) 1];
	for k=1:N;
		int1(k) = wt(k)*HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
		int2(k) = wt(k)*HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
	end
elseif method==4
	% Simpson's 3/8 rule --------------------------------------------
	% Ensure that N-1 is divisible by 3
	N = N-mod(N,3)+1;
	% Build the new grid
	h = (b-a)/(N-1);
	wt = (3*h/8).*[[1 3 3] repmat([2 3 3],1,(N-1)/3-1) 1];
	phi = [a:h:b];
	for k=1:N
		int1(k) = wt(k)*HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,1,trap);
		int2(k) = wt(k)*HestonProb(phi(k),kappa,theta,lambda,rho,sigma,T,K,S,r,q,v0,2,trap);
	end		
end

% The integrals
I1 = sum(int1);
I2 = sum(int2);

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

