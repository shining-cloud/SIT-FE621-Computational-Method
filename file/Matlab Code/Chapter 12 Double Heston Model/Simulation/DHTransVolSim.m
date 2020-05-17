function [S w1 w2 SimPrice] = DHTransVolSim(scheme,params,S0,Strike,Mat,r,q,T,N,PutCall)

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
%   w = Vector of simulated volatilities
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

% Initialize the volatilities and stock processes
w1 = zeros(T,N);
w2 = zeros(T,N);
S  = zeros(T,N);

% Starting values for the variance and stock processes
S(1,:)  = S0;          % Spot price 
w1(1,:) = sqrt(v01);   % Heston initial volatility
w2(1,:) = sqrt(v02);   % Heston initial volatility

% Generate the stock and volatility paths
for i=1:N;
	for t=2:T;
		Zv1 = randn(1);
		Zv2 = randn(1);
		if strcmp(scheme,'ZhuEuler')
			% Euler volatility scheme
			w1(t,i) = w1(t-1,i) + 0.5*kappa1*((theta1 - sigma1^2/4/kappa1)/w1(t-1,i) - w1(t-1,i))*dt ...
				    + 0.5*sigma1*sqrt(dt)*Zv1;
			w2(t,i) = w2(t-1,i) + 0.5*kappa2*((theta2 - sigma2^2/4/kappa2)/w2(t-1,i) - w2(t-1,i))*dt ...
				    + 0.5*sigma2*sqrt(dt)*Zv2;
		elseif strcmp(scheme,'ZhuTV')
			% Zhu (2010) process for the transformed volatility
			m11 = theta1 + (w1(t-1,i)^2 - theta1)*exp(-kappa1*dt);
			m21 = theta2 + (w2(t-1,i)^2 - theta2)*exp(-kappa2*dt);
			m12 = sigma1^2/4/kappa1*(1-exp(-kappa1*dt));
			m22 = sigma2^2/4/kappa2*(1-exp(-kappa2*dt));
			beta1 = sqrt(max(0,m11-m12));
			beta2 = sqrt(max(0,m21-m22));
			thetav1 = (beta1 - w1(t-1,i)*exp(-kappa1*dt/2))/(1-exp(-kappa1*dt/2));
			thetav2 = (beta2 - w2(t-1,i)*exp(-kappa2*dt/2))/(1-exp(-kappa2*dt/2));
			w1(t,i) = w1(t-1,i) + 0.5*kappa1*(thetav1 - w1(t-1,i))*dt + 0.5*sigma1*sqrt(dt)*Zv1;
			w2(t,i) = w2(t-1,i) + 0.5*kappa2*(thetav2 - w2(t-1,i))*dt + 0.5*sigma2*sqrt(dt)*Zv2;
        end
		% Predictor-Corrector for the stock price
		B1 = randn(1);
		B2 = randn(1);
		logS = log(exp(-r*t*dt)*S(t-1,i)) ...
			 + K01 + K11*w1(t-1,i)^2 + K21*w1(t,i)^2 + sqrt(K31*(w1(t,i)^2+w1(t-1,i)^2))*B1 ...
			 + K02 + K12*w2(t-1,i)^2 + K22*w2(t,i)^2 + sqrt(K32*(w2(t,i)^2+w2(t-1,i)^2))*B2;
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

