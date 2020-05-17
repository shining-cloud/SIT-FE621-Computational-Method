function [CallFFT K lambdainc eta] = HestonCallFFTGreek(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,Trap,alpha,fast,rule,Greek)

% Fast Fourier Transform of the Heston Model
% INPUTS
%   N  = number of discretization points
%   uplimit = Upper limit of integration
%   S0 = spot price
%   r = risk free rate
%   q = dividend yield
%   tau = maturity
%   sigma = volatility
%   alpha = dampening factor
%   fast = fast versus slow algorithm.
%     fast = 1 fast version.  Uses vectorization.
%     fast = 0 slow version.  Uses loops
%   rule = integration rule
%     rule = 1 --> Trapezoidal
%     rule = 2 --> Simpson's Rule
% --------------------------------------
% Outputs
%   CallFFT = Black Scholes call prices using FFT
%   CallBS  = Black Scholes call prices using closed form
%         K = Strike prices
%       eta = increment for integration range
% lambdainc = increment for log-strike range
% --------------------------------------
% By Fabrice Douglas Rouah www.Volopta.com

% log spot price
s0 = log(S0);

% Specify the increments
eta = uplimit/N;
lambdainc = 2*pi/N/eta;

% Initialize and specify the weights
w = ones(N,1);
if strcmp(rule,'T')             % Trapezoidal rule
	w(1) = 1/2;
	w(N) = 1/2;
elseif strcmp(rule,'S')         % Simpson's rule
	w(1) = 1/3;
	w(N) = 1/3;
	for k=2:N-1
		if mod(k,2)==0
			w(k) = 4/3;
		else
			w(k) = 2/3;
		end
	end
end

% Specify the b parameter
b = N*lambdainc/2;

% Create the grid for the integration
v = eta.*[0:N-1]';

% Create the grid for the log-strikes
k = -b + lambdainc.*[0:N-1]' + s0;

% Create the strikes and identify ATM 
K = exp(k);

% Initialize the price vector;
CallFFT = zeros(N,1);

if fast==1
	% Implement the FFT - fast algorithm
	U = [0:N-1];
	J = [0:N-1];
	psi = HestonCFGreek(v-(alpha+1).*i,kappa,theta,lambda,rho,sigma,tau,S0,r,q,v0,Trap,Greek);
	phi = exp(-r*tau).*psi ./ (alpha.^2 + alpha - v.^2 + i.*v.*(2*alpha+1));
	x = exp(i.*(b-s0).*v).*phi.*w;
	e = exp(-i*2*pi/N.*(U'*J))*x;
	CallFFT = eta.*exp(-alpha.*k)./pi .* real(e);
elseif fast==0
	% Implement the FFT - slow algorithm
	for u=1:N
		for j=1:N
			psi(j) = HestonCFGreek(v(j)-(alpha+1)*i,kappa,theta,lambda,rho,sigma,tau,S0,r,q,v0,Trap,Greek);
			phi(j) = exp(-r*tau)*psi(j)/(alpha^2 + alpha - v(j)^2 + i*v(j)*(2*alpha+1));
			x(j) = exp(i*(b-s0)*v(j))*phi(j)*w(j);
			e(j) = exp(-i*2*pi/N*(j-1)*(u-1))*x(j);
		end
		CallFFT(u) = eta*exp(-alpha*k(u))/pi * real(sum(e));
	end
end
