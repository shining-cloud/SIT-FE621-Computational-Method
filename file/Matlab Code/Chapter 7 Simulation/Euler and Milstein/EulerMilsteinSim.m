function [S V F] = EulerMilsteinSim(scheme,negvar,params,S0,Mat,r,q,T,N,alpha)

% Euler, Milstein, Implicit Milstein, and Weighted Implicit Milstein
% discretization of the Heston model.
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

% Flags for negative variances
F = 0;

% Starting values for the variance and stock processes
S(1,:) = S0;       % Spot price 
V(1,:) = v0;       % Heston v0 initial variance 

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		% Generate two dependent N(0,1) variables with correlation rho
		Zv = randn(1);
		Zs = rho*Zv + sqrt(1-rho^2)*randn(1); 
		
		if strcmp(scheme,'E')
			% Euler discretization for the variance
			V(t,i) = V(t-1,i) + kappa*(theta-V(t-1,i))*dt ...
			       + sigma*sqrt(V(t-1,i)*dt)*Zv;
		elseif strcmp(scheme,'M')
			% Milstein discretization for the variance.
			V(t,i) = V(t-1,i) + kappa*(theta-V(t-1,i))*dt ...
				   + sigma*sqrt(V(t-1,i)*dt)*Zv ...
				   + (1/4)*sigma^2*dt*(Zv^2-1);
		elseif strcmp(scheme,'IM')
	   		% Implicit Milstein for the variance.
			V(t,i) = (V(t-1,i) + kappa*theta*dt ...
				   + sigma*sqrt(V(t-1,i)*dt)*Zv...
				   + sigma^2*dt*(Zv^2-1)/4) / (1+kappa*dt);
		elseif strcmp(scheme,'WM')
			% Weighted Explicit-Implicit Milstein Scheme
			V(t,i) = (V(t-1,i) + kappa*(theta-alpha*V(t-1,i))*dt ...
				   + sigma*sqrt(V(t-1,i)*dt)*Zv...
				   + sigma^2*dt*(Zv^2-1)/4) / (1+(1-alpha)*kappa*dt);
		end
		
		% Apply the full truncation or reflection scheme to the variance
		if V(t,i) <= 0
			F = F+1;
			if strcmp(negvar,'R')          % Reflection: take -V
				V(t,i) = abs(V(t,i));
			elseif strcmp(negvar,'T')
				V(t,i) = max(0, V(t,i));   % Truncation: take max(0,V)
			end
		end
		
		% Discretize the log stock price
		S(t,i) = S(t-1,i)*exp((r-q-V(t-1,i)/2)*dt + sqrt(V(t-1,i)*dt)*Zs);
	end
end

