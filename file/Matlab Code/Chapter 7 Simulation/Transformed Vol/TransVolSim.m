function [S v] = TransVolSim(scheme,params,S0,Mat,r,q,T,N)

% Zhu (2010) Transformed Volatility discretization of the Heston model
% INPUTS
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

% Time increment
dt = Mat/T;

% Initialize the volatility and stock processes
v = zeros(T,N);
w = zeros(T,N);
S = zeros(T,N);
X = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:) = S0;        % Spot price 
X(1,:) = log(S0);   % Log spot price
v(1,:) = v0;        % Heston initial variance
w(1,:) = sqrt(v0);  % Heston initial volatility

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		% Generate two dependent N(0,1) variables with correlation rho
		Zv = randn(1);
		Zx = rho*Zv + sqrt(1-rho^2)*randn(1);

		if strcmp(scheme,'Euler')
			% Euler volatility scheme
			w(t,i) = w(t-1,i) + 0.5*kappa*((theta-sigma^2/4/kappa)/w(t-1,i) - w(t-1,i))*dt ...
				   + 0.5*sigma*sqrt(dt)*Zv;

		elseif strcmp(scheme,'TV')
			% Transformed Volatility scheme
			m1 = theta + (v(t-1,i) - theta)*exp(-kappa*dt);
			m2 = sigma^2/4/kappa*(1-exp(-kappa*dt));
            beta = sqrt(max(0,m1-m2));
			thetav = (beta - w(t-1,i)*exp(-kappa*dt/2))/(1-exp(-kappa*dt/2));
			w(t,i) = w(t-1,i) + 0.5*kappa*(thetav - w(t-1,i))*dt + 0.5*sigma*sqrt(dt)*Zv;
        end
        v(t,i) = w(t,i)^2;

        % Discretize the log stock price
		X(t,i) = X(t-1,i) + (r-q-w(t-1,i)^2/2)*dt + w(t-1,i)*sqrt(dt)*Zx;
		S(t,i) = exp(X(t,i));
	end
end

