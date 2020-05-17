function [S V] = MMSimGreeks(params,S0,Mat,r,q,T,N,Zv,Zs)

% Andersen and Brotherton-Ratcliffe Moment Matching Scheme.
% INPUTS
%   negvar = 'R' Reflection or 'T' Truncation of negative variances
%   params = Heston parameters
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
%   Zv, Zs = (TxN) matrices of correlated N(0,1) variables
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances
%   F = Matrix of overriden negative variances

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);
lambda = params(6);

% Time increment
dt = Mat/T;

% Initialize the variance and stock processes
V = zeros(T,N);
S = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:) = S0;       % Spot price 
V(1,:) = v0;       % Heston v0 initial variance 

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		% Matched moment lognormal approximation
		num = 0.5*sigma^2*V(t-1,i)*(1-exp(-2*kappa*dt)) / kappa;
		den = (exp(-kappa*dt)*V(t-1,i) + (1-exp(-kappa*dt))*theta)^2;
		Gam = log(1 + num/den);
		V(t,i) = (exp(-kappa*dt)*V(t-1,i) + (1-exp(-kappa*dt))*theta)...
			   * exp(-0.5*Gam^2 + Gam*Zv(t,i));
		% Euler/Milstein discretization scheme for the log stock prices
		S(t,i) = S(t-1,i)*exp((r-q-V(t-1,i)/2)*dt + sqrt(V(t-1,i)*dt)*Zs(t,i));
	end
end

