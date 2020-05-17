function [S V] = AlfonsiSim(params,S0,Mat,r,q,T,N)

% Alfonsi discretization of the Double Heston model.
% INPUTS
%   params = Heston parameters
%   S0  = Spot Price
%   Mat = Maturity
%   r = riskf ree rate
%   q = dividend yield
%   T   = Number of time steps
%   N   = Number of stock price paths
% OUTPUTS
%   S = Vector of simulated stock prices
%   V = Vector of simulated variances

% Heston parameters
kappa = params(1);
theta = params(2);
sigma = params(3);
v0    = params(4);
rho   = params(5);

% Time increment
dt = Mat/T;

% Required quantities
K0 = -rho*kappa*theta*dt/sigma;
K1 = dt/2*(kappa*rho/sigma - 1/2) - rho/sigma;
K2 = dt/2*(kappa*rho/sigma - 1/2) + rho/sigma;
K3 = dt/2*(1-rho^2);

% Initialize the variance and stock processes
V = zeros(T,N);
S = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:) = S0;       % Spot price 
V(1,:) = v0;       % Heston v0 initial variance 

% Generate the stock and volatility paths
for i=1:N;
    for t=2:T;
        % Alfonsi discretization
        V(t,i) = AlfonsiV(params,V(t-1,i),dt);
        % Predictor-Corrector for the stock price
        B = randn(1);
        logS = log(exp(-r*t*dt)*S(t-1,i)) ...
            + K0 + K1*V(t-1,i) + K2*V(t,i) + sqrt(K3*(V(t,i)+V(t-1,i)))*B;
        S(t,i) = exp(logS)*exp(r*(t+1)*dt);
    end
end
