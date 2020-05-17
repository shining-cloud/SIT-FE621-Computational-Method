function [CallFRFT K lambdainc eta] = HestonCallFRFT(N,uplimit,S0,r,q,tau,kappa,theta,lambda,rho,sigma,v0,Trap,alpha,rule,eta,lambdainc);

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
%   rule = integration rule
%     rule = 1 --> Trapezoidal
%     rule = 2 --> Simpson's Rule
%   eta = integration range increment
%   lambdainc = strike range increment
% --------------------------------------
% OUTPUTS
%   CallFFT = Heston call prices using FFT
%   CallBS  = Heston call prices using closed form
%         K = Strike prices
%       eta = increment for integration range
% lambdainc = increment for log-strike range
% --------------------------------------

% log spot price
s0 = log(S0);

% Initialize and specify the weights
w = zeros(N,1);
if strcmp(rule,'T')             % Trapezoidal rule
	w(1) = 1/2;
	w(N) = 1/2;
	w(2:N-1) = 1;
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
CallFRFT = zeros(N,1);

% Parameter for the FRFT
beta = lambdainc*eta/2/pi;

% Implement the FRFT - fast algorithm
U = [0:N-1];
J = [0:N-1];
psi = HestonCF(v-(alpha+1).*i,kappa,theta,lambda,rho,sigma,tau,S0,r,q,v0,Trap);
psi = conj(psi);
phi = exp(-r*tau).*psi./conj(alpha.^2 + alpha - v.^2 + i.*v.*(2*alpha+1));
x = conj(exp(i.*(b-s0).*v)).*phi.*w;
y = real(FRFT(x',beta));
CallFRFT = eta.*exp(-alpha.*k).*y'./pi;

