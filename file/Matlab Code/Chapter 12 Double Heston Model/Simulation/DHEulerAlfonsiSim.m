function [S V1 V2 SimPrice] = DHEulerAlfonsiSim(scheme,params,S0,Strike,Mat,r,q,T,N,PutCall)

% Euler, Milstein, Implicit Milstein, and Weighted Implicit Milstein
% discretization of the Double Heston model.
% INPUTS
%   scheme = Scheme for the variance and stock price
%            'E' Euler, 'M' Milstein, 'IM' Implicit Milstein, or
%            'WM' Weighted Explicit-Implicit Scheme
%   negvar = 'R' Reflection or 'T' Truncation of negative variances
%   params = Heston parameters
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
%   alpha = Weight for explicit-implicit scheme
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances
%   F = Matrix of overriden negative variances

% Heston parameters
kappa1 = params(1);
theta1 = params(2);
sigma1 = params(3);
v01    = params(4);
rho1   = params(5);

kappa2 = params(6);
theta2 = params(7);
sigma2 = params(8);
v02    = params(9);
rho2   = params(10);

params1 = params(1:5);
params2 = params(6:10);

% Time increment
dt = Mat/T;

% Required quantities
K01 = -rho1*kappa1*theta1*dt/sigma1;
K11 = dt/2*(kappa1*rho1/sigma1 - 1/2) - rho1/sigma1;
K21 = dt/2*(kappa1*rho1/sigma1 - 1/2) + rho1/sigma1;
K31 = dt/2*(1-rho1^2);

K02 = -rho2*kappa2*theta2*dt/sigma2;
K12 = dt/2*(kappa2*rho2/sigma2 - 1/2) - rho2/sigma2;
K22 = dt/2*(kappa2*rho2/sigma2 - 1/2) + rho2/sigma2;
K32 = dt/2*(1-rho2^2);

% Initialize the variance and stock processes
V1 = zeros(T,N);
V2 = zeros(T,N);
S  = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:)  = S0;       % Spot price 
V1(1,:) = v01;       % Heston v0 initial variance 
V2(1,:) = v02;       % Heston v0 initial variance 

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		if strcmp(scheme, 'Euler')
			% Generate two in dependent N(0,1) variables
			G1 = randn(1);
			G2 = randn(1);
			% Euler discretization with full truncation for the variances
			V1(t,i) = V1(t-1,i) + kappa1*(theta1-V1(t-1,i))*dt + sigma1*sqrt(V1(t-1,i)*dt)*G1;
			V2(t,i) = V2(t-1,i) + kappa2*(theta2-V2(t-1,i))*dt + sigma2*sqrt(V2(t-1,i)*dt)*G2;
			V1(t,i) = max(0,V1(t,i));
			V2(t,i) = max(0,V2(t,i));
		elseif strcmp(scheme,'Alfonsi')
			% Alfonsi discretization
			V1(t,i) = AlfonsiV(params1,V1(t-1,i),dt);
			V2(t,i) = AlfonsiV(params2,V2(t-1,i),dt);
		end
		% Predictor-Corrector for the stock price
		B1 = randn(1);
		B2 = randn(1);
		logS = log(exp(-r*t*dt)*S(t-1,i)) ...
			 + K01 + K11*V1(t-1,i) + K21*V1(t,i) + sqrt(K31*(V1(t,i)+V1(t-1,i)))*B1 ...
			 + K02 + K12*V2(t-1,i) + K22*V2(t,i) + sqrt(K32*(V2(t,i)+V2(t-1,i)))*B2;
		S(t,i) = exp(logS)*exp(r*(t+1)*dt);
	end
end

% Terminal stock prices
ST = S(end,:);

% Payoff vectors
if strcmp(PutCall,'C')
	Payoff = max(ST - Strike,0);
elseif strcmp(PutCall,'P')
	Payoff = max(Strike - ST,0);
end
	
% Simulated price
SimPrice = exp(-r*Mat)*mean(Payoff);
